/** @file
 * Experimental functions which are liable to
 * API breaks within QuEST minor version releases.
 * Some optional functions require compiling this
 * file against MPI, despite being outside of /comm/, 
 * and so require opt-in macros (QUEST_COMPILE_SUBCOMM)
 * 
 * @author Oliver Brown
 */

#include "quest/include/config.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/modes.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#if QUEST_COMPILE_SUBCOMM && ! QUEST_COMPILE_MPI
    #error "Macro QUEST_COMPILE_SUBCOMM was true, but QUEST_COMPILE_MPI was illegally false."
#endif

#if QUEST_COMPILE_SUBCOMM
    #include <mpi.h>
#endif

#ifdef QUEST_COMPILE_ADIOS2
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
auto createAdios() {

    // In distributed builds, ADIOS2 must be given QuEST's communicator so that each
    // node's call collectively writes/reads its own slice of the shared file. Without
    // it, ADIOS2 runs serially per rank and the per-node slices never form one file.
    #if QUEST_COMPILE_MPI
        return adios2::ADIOS(comm_getMpiComm());
    #else
        return adios2::ADIOS();
    #endif
}
#endif



/*
 * API FUNCTIONS
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



    // TODO:
    // - make comment about size_t overflow risk
    // - fix Qureg{} return warning issue
    // - check restoration to a DISTRIBUTED qureg is correct


void saveQuregToFile(Qureg qureg, const char* fn) {
    validate_adios2IsCompiled(__func__);
    validate_quregFields(qureg, __func__);

#ifdef QUEST_COMPILE_ADIOS2

    // Pedantic but safe - don't let ADIOS2 start reading amps prematurely
    if (qureg.isDistributed)
        comm_sync();
    
    // TODO:
    // We can optimise in GPU settings by giving ADIOS2 the device memory
    // pointers; but for now, we simply stage into CPU memory first
    if (qureg.isGpuAccelerated)
        gpu_copyGpuToCpu(qureg);

    // gratuitously re-create ADIOS2 at every call, for simplicity (occluded by IO)
    adios2::ADIOS adios = createAdios();
    adios2::IO io = adios.DeclareIO("QuESTQuregSave");

    // attempt to open the file
    adios2::Engine engine; // default ctor
    try {
        engine = io.Open(fn, adios2::Mode::Write);
    } catch (...) {
        validate_adiosCanOpenFile(false, fn, __func__);
    }

    // global single-value metadata; we deliberately record only the dimension
    // and precision, never incidental deployment fields (the loader chooses its
    // own deployment) nor derivable fields (like numAmps)
    adios2::Variable<int> vNumQubits  = io.DefineVariable<int>("numQubits");
    adios2::Variable<int> vIsDensMatr = io.DefineVariable<int>("isDensityMatrix");
    adios2::Variable<size_t> vQrealBytes = io.DefineVariable<size_t>("qrealBytes"); // also encodes precision

    // amplitudes are stored as interleaved (real, imag) reals to stay agnostic
    // to precision and to ADIOS2's complex-type support; each node writes only
    // its local slice into the global array, avoiding excessive memory use
    // (these scalars are guaranteed not to overflow by createQureg validation)
    qindex globalReals = 2 * qureg.numAmps;
    qindex localReals  = 2 * qureg.numAmpsPerNode;
    qindex startReal   = 2 * ((qindex) qureg.rank) * qureg.numAmpsPerNode;
    adios2::Variable<qreal> vAmpComponents = io.DefineVariable<qreal>(
        "ampComponents",
        { (size_t) globalReals },
        { (size_t) startReal },
        { (size_t) localReals });

    // attempt to write to file
    try {
        engine.BeginStep();
        engine.Put(vNumQubits,  qureg.numQubits);
        engine.Put(vIsDensMatr, qureg.isDensityMatrix);
        engine.Put(vQrealBytes, sizeof(qreal));
        engine.Put(vAmpComponents, reinterpret_cast<qreal*>(qureg.cpuAmps));
        engine.EndStep();
        engine.Close();
    } catch (...) {
        // no need for a finally; RAII frees engine
        validate_adiosCanWriteToFile(false, fn, __func__);
    }

#endif
}


Qureg createQuregFromFile(const char* fn) {
    validate_adios2IsCompiled(__func__);

#ifdef QUEST_COMPILE_ADIOS2

    // gratuitously re-create ADIOS2 at every call, for simplicity (occluded by IO)
    adios2::ADIOS adios = createAdios();
    adios2::IO io = adios.DeclareIO("QuESTQuregLoad");

    // attempt to open the file, and prepare to parse
    adios2::Engine engine; // default ctor
    try {
        engine = io.Open(fn, adios2::Mode::Read);
        engine.BeginStep();
    } catch (...) {
        validate_adiosCanOpenFile(false, fn, __func__);
    }

    // check that the file contains the expected variables
    auto vNumQubits  = io.InquireVariable<int>("numQubits");
    auto vIsDensMatr = io.InquireVariable<int>("isDensityMatrix");
    auto vQrealBytes = io.InquireVariable<size_t>("qrealBytes");
    auto vAmpComponents = io.InquireVariable<qreal>("ampComponents");
    bool areAllVarsPresent = vNumQubits && vIsDensMatr && vQrealBytes && vAmpComponents;
    validate_adiosFileContainsFields(areAllVarsPresent, __func__);

    // read dimension + precision metadata first, so we can size the new Qureg
    int numQubits = 0;
    int isDensMatr = 0;
    size_t fileQrealBytes = 0;
    try {
        engine.Get(vNumQubits,  numQubits);
        engine.Get(vIsDensMatr, isDensMatr);
        engine.Get(vQrealBytes, fileQrealBytes);
        engine.PerformGets();
    } catch(...) {
        validate_adiosCanReadFile(false, fn, __func__);
    }

    // check the amps are of the expected precision, and so are parsable
    validate_newQuregFileMatchesPrecision(fileQrealBytes, __func__);

    // attempt to create a matching-dimension Qureg with automatically chosen deployments
    Qureg qureg = validateAndCreateCustomQureg(numQubits, isDensMatr, 
        modeflag::USE_AUTO, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);

    // read only this node's slice of the global amplitude array into its buffer
    // (guaranteed not to overflow by above validateAndCreateCustomQureg validation)
    qindex localReals = 2 * qureg.numAmpsPerNode;
    qindex startReal  = 2 * ((qindex) qureg.rank) * qureg.numAmpsPerNode;
    vAmpComponents.SetSelection({ { (size_t) startReal }, { (size_t) localReals } });
    try {
        engine.Get(vAmpComponents, reinterpret_cast<qreal*>(qureg.cpuAmps)); // immediate; PerformGets redundant
    } catch(...) {
        validate_adiosCanReadFile(false, fn, __func__);
    }

    // complete ADIOS2 work
    try {
        engine.EndStep();
        engine.Close();
    } catch(...) {
        validate_adiosCanReadFile(false, fn, __func__);
    }

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
