/** @file
 * Experimental functions which are liable to
 * API breaks within QuEST minor version releases.
 * Some optional functions require compiling this
 * file against MPI, despite being outside of /comm/, 
 * and so require opt-in macros (QUEST_COMPILE_SUBCOMM)
 * 
 * @author Oliver Brown (custom QuESTEnv)
 * @author Ashmit JaiSarita Gupta (checkpointing)
 * @author Tyson Jones (structure, validation)
 */

#include "quest/include/config.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/modes.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <string>

#if QUEST_COMPILE_SUBCOMM && ! QUEST_COMPILE_MPI
    #error "Macro QUEST_COMPILE_SUBCOMM was true, but QUEST_COMPILE_MPI was illegally false."
#endif

#if QUEST_COMPILE_SUBCOMM
    #include <mpi.h>
#endif

#if QUEST_COMPILE_ADIOS2
    #include <adios2.h>

    #if QUEST_COMPILE_MPI
        #include <mpi.h>
    #endif
#endif



/*
 * EXTERNAL FUNCTIONS
 *
 * which we here regretfully 'extern' because we are either
 * unsure which header should expose them, or because they
 * contain deployment-specific types (like MPI_Comm) which
 * we do not wish to expose within internal headers 
 */


extern void validateAndInitCustomQuESTEnv(
    int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread, const char* caller);


extern Qureg validateAndCreateCustomQureg(
    int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, const char* caller);


#if QUEST_COMPILE_SUBCOMM // hide MPI_Comm
    extern bool comm_setMpiComm(MPI_Comm newComm, bool userOwnsMpi);
#endif


#if (QUEST_COMPILE_ADIOS2 && QUEST_COMPILE_MPI) // hide MPI_Comm
    extern MPI_Comm comm_getMpiComm();
#endif



/*
 * INTERNAL FUNCTIONS
 */


#if QUEST_COMPILE_ADIOS2
auto createAdios(bool useMpi) {

    // suppress unused warning when MPI not compiled (implies useMpi=false)
    (void) useMpi;

    // When the Qureg is distributed, ADIOS2 must be given QuEST's communicator so that each
    // node writes/reads its own slice of the shared file
    #if QUEST_COMPILE_MPI
        return useMpi?
            adios2::ADIOS(comm_getMpiComm()) :
            adios2::ADIOS();
    #else
        return adios2::ADIOS(); // implies useMpi=0
    #endif

    // caller need not call destructor; ADIOS2 uses RAII
}
#endif



/*
 * C API FUNCTIONS
 */


// enable invocation by both C and C++ binaries
extern "C" {


void initCustomMpiQuESTEnv(int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread) {
    validateAndInitCustomQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread, __func__);
}


#if QUEST_COMPILE_SUBCOMM // hide MPI_Comm
void initCustomMpiCommQuESTEnv(MPI_Comm userQuestComm, int useGpuAccel, int useMultithread) {

    // useDistrib and userOwnsMpi are implied by the user of this initialiser
    const int useDistrib = 1;
    const bool userOwnsMpi = true;

    // pre-validate that we are able to set the MPI communicator
    validate_mpiInitStatus(useDistrib, userOwnsMpi, __func__);
    validate_mpiSubCommIsNonNull(userQuestComm != MPI_COMM_NULL, __func__);

    // avoid re-setting the MPI comm (to avoid an internal error), which happens
    // if a user illegally re-calls this function, which will be subsequently
    // caught by the validation in validateAndInitCustomQuESTEnv() below
    if (!comm_isActive()) {
        bool success = comm_setMpiComm(userQuestComm, userOwnsMpi);
        validate_mpiSubCommSetSucceeded(success, __func__);
    }

    // perform remaining validation (some is harmlessly repeated) and init QuEST env
    validateAndInitCustomQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread, __func__);
}
#endif // QUEST_COMPILE_SUBCOMM


int getQuESTNumGpuThreadsPerBlock() {
    validate_envIsInit(__func__);
    
    return gpu_getNumThreadsPerBlock();
}


void setQuESTNumGpuThreadsPerBlock(int numTPB) {
    validate_envIsInit(__func__);

    // validation messages and queries depend upon GPU usage
    bool gpuIsActive = getQuESTEnv().isGpuAccelerated;
    validate_numGpuThreadsPerBlock(numTPB, gpuIsActive, __func__);

    gpu_setNumThreadsPerBlock(numTPB);
}


void saveQuregToFile(Qureg qureg, const char* fn) {
    validate_adios2IsCompiled(__func__);
    validate_quregFields(qureg, __func__);
    
    (void) fn; // suppress unused warning

#if QUEST_COMPILE_ADIOS2

    // When the QuEST env is distributed, but the given Qureg is not (and is instead
    // duplicated upon every node), we should permit only a single process (the root)
    // to use ADIOS2 to write to the (assumably, shared) filesystem. Note that non-root
    // nodes must not exit; they need to participate in validation syncs
    bool shouldSkipAdios = (! qureg.isDistributed) && (comm_getRank() > ROOT_RANK);

    // pedantic but safe - don't let ADIOS2 start reading amps prematurely
    if (qureg.isDistributed)
        comm_sync();
    
    // TODO:
    // We can optimise in GPU settings by giving ADIOS2 the device memory
    // pointers; but for now, we simply stage into CPU memory first
    if (qureg.isGpuAccelerated)
        gpu_copyGpuToCpu(qureg);

    // gratuitously re-create ADIOS2 at every call, for simplicity (occluded by IO)
    adios2::ADIOS adios = createAdios(qureg.isDistributed);
    adios2::IO io = adios.DeclareIO("QuESTQuregSave");

    // use BP5 specifically to avoid non-root-hangs upon rank exceptions
    io.SetEngine("BP5");

    // attempt to open the file
    adios2::Engine engine; // default ctor
    bool success = false;
    try {
        if (!shouldSkipAdios)
            engine = io.Open(fn, adios2::Mode::Write);
        success = true;
    } catch (...) {}
    validate_adiosCanOpenFileOnAllNodes(success, fn, __func__);

    // global single-value metadata; we deliberately record only the dimension
    // and precision, never incidental deployment fields (the loader chooses its
    // own deployment) nor derivable fields (like numAmps)
    adios2::Variable<int> vNumQubits  = io.DefineVariable<int>("numQubits");
    adios2::Variable<int> vNumNodes   = io.DefineVariable<int>("numNodes");
    adios2::Variable<int> vIsDensMatr = io.DefineVariable<int>("isDensityMatrix");
    adios2::Variable<size_t> vQrealBytes = io.DefineVariable<size_t>("qrealBytes"); // also encodes precision

    // amplitudes are stored as interleaved (real, imag) reals to stay agnostic
    // to precision and to ADIOS2's complex-type support; each node writes only
    // its local slice into the global array, avoiding excessive memory use
    // (these scalars are guaranteed not to overflow by createQureg validation)
    qindex globalReals = 2 * qureg.numAmps;
    qindex localReals  = 2 * qureg.numAmpsPerNode;
    qindex startReal   = localReals * qureg.rank;
    adios2::Variable<qreal> vAmpComponents = io.DefineVariable<qreal>(
        "ampComponents",
        { (size_t) globalReals },
        { (size_t) startReal },
        { (size_t) localReals });

    // attempt to write to file
    success = false;
    try {
        if (!shouldSkipAdios) {
            engine.Put(vNumQubits,  qureg.numQubits);
            engine.Put(vNumNodes,   qureg.numNodes);
            engine.Put(vIsDensMatr, qureg.isDensityMatrix);
            engine.Put(vQrealBytes, sizeof(qreal));
            engine.Put(vAmpComponents, reinterpret_cast<qreal*>(qureg.cpuAmps));
            engine.Close();
        }
        success = true;
    } catch (...) {}
    validate_adiosCanWriteToFileOnAllNodes(success, fn, __func__);

    // prevent any process from continuing until ADIOS2 is fully finished
    if (qureg.isDistributed)
        comm_sync();

#endif
}


Qureg createQuregFromFile(const char* fn) {
    validate_adios2IsCompiled(__func__);

#if QUEST_COMPILE_ADIOS2

    // pedantic but safe - don't let ADIOS2 start reading while other processes are working
    if (comm_isActive())
        comm_sync();

    // make ADIOS2 MPI-aware even when the subsequently-loaded Qureg is
    // auto-deployed to be non-distributed; every process will safely
    // parse the file and independently update its Qureg copy
    bool giveAdiosMpi = comm_isActive();

    // gratuitously re-create ADIOS2 at every call, for simplicity (occluded by IO)
    adios2::ADIOS adios = createAdios(giveAdiosMpi);
    adios2::IO io = adios.DeclareIO("QuESTQuregLoad");

    // use BP5 specifically to avoid non-root-hangs upon rank exceptions
    io.SetEngine("BP5");

    // attempt to open the file, and prepare to parse
    adios2::Engine engine; // default ctor
    bool success = false;
    try {
        engine = io.Open(fn, adios2::Mode::ReadRandomAccess);
        success = true;
    } catch (...) {}
    validate_adiosCanOpenFileOnAllNodes(success, fn, __func__);

    // check that the file contains the expected variables
    auto vNumQubits  = io.InquireVariable<int>("numQubits");
    auto vNumNodes   = io.InquireVariable<int>("numNodes");
    auto vIsDensMatr = io.InquireVariable<int>("isDensityMatrix");
    auto vQrealBytes = io.InquireVariable<size_t>("qrealBytes");
    auto vAmpComponents = io.InquireVariable<qreal>("ampComponents");
    bool areAllVarsPresent = vNumQubits && vNumNodes && vIsDensMatr && vQrealBytes && vAmpComponents;
    validate_adiosFileContainsFieldsOnAllNodes(areAllVarsPresent, __func__);

    // read dimension + precision metadata first, so we can size the new Qureg
    int numQubits = 0;
    int numNodes = 0;
    int isDensMatr = 0;
    size_t fileQrealBytes = 0;
    success = false;
    try {
        engine.Get(vNumQubits,  numQubits);
        engine.Get(vNumNodes,   numNodes);
        engine.Get(vIsDensMatr, isDensMatr);
        engine.Get(vQrealBytes, fileQrealBytes);
        engine.PerformGets();
        success = true;
    } catch (...) {}
    validate_adiosCanReadFileOnAllNodes(success, fn, __func__);

    // check the amps are of the expected precision, and so are parsable
    validate_newQuregFileMatchesPrecision(fileQrealBytes, __func__);

    // attempt to create a matching-dimension Qureg with automatically chosen deployments
    Qureg qureg = validateAndCreateCustomQureg(numQubits, isDensMatr, 
        modeflag::USE_AUTO, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);

    // auto-distribution MUST match checkpointed distribution (pre-free to avoid leak)
    if (qureg.numNodes != numNodes)
        destroyQureg(qureg);
    validate_newQuregNumNodesMatchesSavedFile(numNodes, qureg.numNodes, comm_getNumNodes(), numQubits, isDensMatr, __func__);

    // read only this node's slice of the global amplitude array into its buffer
    // (guaranteed not to overflow by above validateAndCreateCustomQureg validation)
    qindex localReals = 2 * qureg.numAmpsPerNode;
    qindex startReal  = 2 * ((qindex) qureg.rank) * qureg.numAmpsPerNode;
    vAmpComponents.SetSelection({ { (size_t) startReal }, { (size_t) localReals } });
    success = false;
    try {
        engine.Get(vAmpComponents, reinterpret_cast<qreal*>(qureg.cpuAmps)); // immediate; PerformGets redundant
        success = true;
    } catch (...) {}
    validate_adiosCanReadFileOnAllNodes(success, fn, __func__);

    // complete ADIOS2 work
    success = false;
    try {
        engine.Close();
        success = true;
    } catch (...) {}
    validate_adiosCanReadFileOnAllNodes(success, fn, __func__);

    // propagate the restored CPU amplitudes to the GPU, if deployed
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);

    return qureg;
#else
    // unreachable: the validation above always throws in non-checkpointing builds
    return Qureg{};
#endif
}


// end de-mangler
}


/*
 * C++ API FUNCTIONS
 */

void saveQuregToFile(Qureg qureg, std::string fn) {

    saveQuregToFile(qureg, fn.c_str());
}

Qureg createQuregFromFile(std::string fn) {

    return createQuregFromFile(fn.c_str());
}
