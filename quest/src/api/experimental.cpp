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


#if QUEST_COMPILE_SUBCOMM // hide MPI_Comm
    extern bool comm_setMpiComm(MPI_Comm newComm, bool userOwnsMpi);
#endif



/*
 * INTERNAL FUNCTIONS
 */


// TODO:
// below is broken; we must not give COMM_WORLD, but instead the QuEST
// subcommunicator. Must get this from comm somehow, though this requires
// exposing an MPI type across QuEST translation units. Hmm!!!


#if QUEST_COMPILE_ADIOS2
// In distributed builds, ADIOS2 must be given QuEST's communicator so that each
// node's call collectively writes/reads its own slice of the shared file. Without
// it, ADIOS2 runs serially per rank and the per-node slices never form one file.
static adios2::ADIOS makeAdios() {
#if QUEST_COMPILE_MPI
    return adios2::ADIOS(MPI_COMM_WORLD);
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


void saveQuregToFile(Qureg qureg, const char* fn) {
    validate_quregCheckpointingIsCompiled(__func__);

#ifdef QUEST_COMPILE_ADIOS2
    validate_quregFields(qureg, __func__);

    // ensure the CPU amplitudes reflect any GPU-resident state before writing
    syncQuregFromGpu(qureg);

    adios2::ADIOS adios = makeAdios();
    adios2::IO io = adios.DeclareIO("QuESTQuregSave");
    adios2::Engine engine = io.Open(fn, adios2::Mode::Write);

    // global single-value metadata; we deliberately record only the dimension
    // and precision, never incidental deployment fields (the loader chooses its
    // own deployment) nor derivable fields (like numAmps)
    adios2::Variable<int> vNumQubits  = io.DefineVariable<int>("numQubits");
    adios2::Variable<int> vIsDensMatr = io.DefineVariable<int>("isDensityMatrix");
    adios2::Variable<int> vQrealBytes = io.DefineVariable<int>("qrealBytes");

    // amplitudes are stored as interleaved (real, imag) reals to stay agnostic
    // to precision and to ADIOS2's complex-type support; each node writes only
    // its local slice into the global array, avoiding excessive memory use
    qindex globalReals = 2 * qureg.numAmps;
    qindex localReals  = 2 * qureg.numAmpsPerNode;
    qindex startReal   = 2 * ((qindex) qureg.rank) * qureg.numAmpsPerNode;
    adios2::Variable<qreal> vAmps = io.DefineVariable<qreal>(
        "amps",
        { (size_t) globalReals },
        { (size_t) startReal },
        { (size_t) localReals });

    int qrealBytes = (int) sizeof(qreal);

    engine.BeginStep();
    engine.Put(vNumQubits,  qureg.numQubits);
    engine.Put(vIsDensMatr, qureg.isDensityMatrix);
    engine.Put(vQrealBytes, qrealBytes);
    engine.Put(vAmps, reinterpret_cast<qreal*>(qureg.cpuAmps));
    engine.EndStep();
    engine.Close();
#endif
}


Qureg createQuregFromFile(const char* fn) {
    validate_quregCheckpointingIsCompiled(__func__);

#ifdef QUEST_COMPILE_ADIOS2
    adios2::ADIOS adios = makeAdios();
    adios2::IO io = adios.DeclareIO("QuESTQuregLoad");
    adios2::Engine engine = io.Open(fn, adios2::Mode::Read);

    engine.BeginStep();

    // read dimension + precision metadata first, so we can size the new Qureg
    int numQubits = 0;
    int isDensMatr = 0;
    int fileQrealBytes = 0;
    engine.Get(io.InquireVariable<int>("numQubits"),       numQubits);
    engine.Get(io.InquireVariable<int>("isDensityMatrix"), isDensMatr);
    engine.Get(io.InquireVariable<int>("qrealBytes"),      fileQrealBytes);
    engine.PerformGets();

    validate_quregFileMatchesPrecision(fileQrealBytes, __func__);

    // create a matching-dimension Qureg with automatically chosen deployments,
    // independent of those used when the file was saved
    Qureg qureg = (isDensMatr)?
        createDensityQureg(numQubits) :
        createQureg(numQubits);

    // read only this node's slice of the global amplitude array into its buffer
    qindex localReals = 2 * qureg.numAmpsPerNode;
    qindex startReal  = 2 * ((qindex) qureg.rank) * qureg.numAmpsPerNode;
    adios2::Variable<qreal> vAmps = io.InquireVariable<qreal>("amps");
    vAmps.SetSelection({ { (size_t) startReal }, { (size_t) localReals } });
    engine.Get(vAmps, reinterpret_cast<qreal*>(qureg.cpuAmps));

    engine.EndStep();
    engine.Close();

    // propagate the restored CPU amplitudes to the GPU, if deployed
    syncQuregToGpu(qureg);

    return qureg;
#else
    // unreachable: the validation above always throws in non-checkpointing builds
    return Qureg{};
#endif

    // TODO: fix above!!! Will warn non-init?
}


// end de-mangler
}
