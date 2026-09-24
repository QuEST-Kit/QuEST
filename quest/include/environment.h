/** @file
 * API signatures for managing the QuEST
 * execution environment.
 * 
 * @author Tyson Jones
 * @author Richard Meister (aided in design)
 * 
 * @defgroup environment Environment
 * @ingroup api
 * @brief Data structures for managing the QuEST execution environment.
 * @{
 */

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <stdbool.h>

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/*
 * QuESTEnv is a struct of which there will be a single, immutable
 * main instance, statically instantiated inside environment.cpp,
 * accessible anywhere via a getter, and which is consulted for
 * determining the deployment configuration. Users can obtain a
 * local copy of this struct with getQuESTEnv().
 */

/// @notyetdoced
typedef struct {

    // deployment modes which can be runtime disabled
    bool isMultithreaded;
    bool isGpuAccelerated;
    bool isDistributed;
    bool isMpiUserOwned;

    // deployment modes which cannot be directly changed after compilation
    bool isCuQuantumEnabled;

    // deployment configurations which can be changed via environment variables
    int isGpuSharingEnabled;
    int isMpiGpuAware;

    // distributed configuration
    int rank;
    int numNodes;

} QuESTEnv;


/** Initialises the QuEST execution environment.
 * 
 * This must be called before any other QuEST function, and performs tasks
 * like validating the environment, reading environment variables,
 * initialising external libraries like MPI or cuQuantum (when available),
 * and seeding random number generators.
 * 
 * This function prepares usage of all of QuEST's parallelisation facilities, such
 * as multithreading, GPU-acceleration and distribution, provided they are compiled
 * and appropriate hardware is available. The used facilities can be controlled with
 * initCustomQuESTEnv().
 * 
 * > [!TIP]
 * > The utilised facilities can be conveniently viewed with reportQuESTEnv().
 * 
 * When distributed execution is initialised with this function, QuEST takes control
 * of MPI, including its initialisation and finalization. User-owned MPI is possible
 * through initCustomMpiQuESTEnv().
 * 
 * Note that when cuQuantum was compiled, and a GPU is available at runtime, the
 * cuQuantum backend is always used over the custom GPU backend (which, infact, was
 * not compiled!). This means the GPU _must_ be compatible with cuQuantum.
 * 
 * > [!NOTE]
 * > Before exiting, the initialised QuEST environment should be finalized with
 * > finalizeQuESTEnv(). This is especially important in a distributed environment
 * > to avoid MPI errors.
 * 
 * @myexample
 * 
 * ```cpp
    int main() {
        initQuESTEnv();
        reportQuESTEnv();
        finalizeQuESTEnv();
        return 0;
    }
 * ```
 * 
 * @throws @validationerror
 * - if the QuEST environment was already initialised, or has already been finalised.
 * - if any environment variable has an invalid value.
 * - if distribution is enabled but MPI was already initialised.
 * - if distribution is enabled but QuEST is launched with a non-power-of-2 number of MPI processes.
 * - if distribution and GPU are enabled, and a GPU is used by more than one process, unless
 *   explicitly enabled through environment variable QUEST_PERMIT_NODES_TO_SHARE_GPU.
 * - if GPU is enabled and cuQuantum was compiled, but the GPU is not compatible with cuQuantum.
 * @see
 * - initCustomQuESTEnv()
 * - initCustomMpiQuESTEnv()
 * - initCustomMpiCommQuESTEnv()
 * - finalizeQuESTEnv()
 * - reportQuESTEnv()
 * @author Tyson Jones
 */
void initQuESTEnv();


/** Initialises the QuEST execution environment with the specified deployments.
 *
 * Each deployment flag may be @c 1 to force the deployment, @c 0 to disable it,
 * or @c -1 to let QuEST choose automatically. The environment must be initialised
 * exactly once, and cannot be re-initialised after finalizeQuESTEnv().
 *
 * @param[in] useDistrib      whether to force (@c =1), disable (@c =0), or automate (@c =-1) distribution.
 * @param[in] useGpuAccel     whether to force (@c =1), disable (@c =0), or automate (@c =-1) GPU acceleration.
 * @param[in] useMultithread  whether to force (@c =1), disable (@c =0), or automate (@c =-1) multithreading.
 * @throws @validationerror
 * - if any deployment flag is not @c 0, @c 1 or @c -1.
 * - if the QuEST environment was already initialised, or has already been finalised.
 * - if any environment variable has an invalid value.
 * - if distribution is enabled but MPI was already initialised.
 * - if distribution is enabled but QuEST is launched with a non-power-of-2 number of MPI processes.
 * - if distribution and GPU are enabled, and a GPU is used by more than one process, unless
 *   explicitly enabled through environment variable QUEST_PERMIT_NODES_TO_SHARE_GPU.
 * - if GPU is enabled and cuQuantum was compiled, but the GPU is not compatible with cuQuantum.
 * @see
 * - initQuESTEnv()
 * - initCustomMpiQuESTEnv()
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_environments.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_environments.cpp) examples
 * @author Tyson Jones
 */
void initCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread);


/** Finalises the active QuEST execution environment.
 *
 * This synchronises outstanding GPU/MPI work, clears QuEST's GPU cache, finalises
 * cuQuantum if active, and finalises MPI if QuEST initialised it. It does not
 * destroy any existing QuEST structs, such as Qureg or CompMatr, which should be
 * prior destroyed to avoid a leak.
 * 
 * @throws @validationerror
 * - if the QuEST environment is not initialised.
 * @see
 * - initQuESTEnv()
 * @author Tyson Jones
 */
void finalizeQuESTEnv();


/** Synchronises QuEST across all processes and machines, waiting for outstanding work to complete.
 * 
 * - When GPU acceleration is active, this function blocks until all outstanding GPU work is complete.
 * - When distribution is active, this function blocks until all MPI ranks are synchronised.
 * 
 * @throws @validationerror
 * - if the QuEST environment is not initialised.
 * @author Tyson Jones
 */
void syncQuESTEnv();


/** Prints a summary of the active QuEST execution environment.
 *
 * The report includes precision, compilation, deployment, CPU, GPU, distribution,
 * Qureg size-limit and automatic-deployment information.
 * 
 * @myexample
 * 
 * An example output:
 * 
 * ```text
    QuEST execution environment:
    [precision]
        qreal.................double (8 bytes)
        qcomp.................std::__1::complex<double> (16 bytes)
        qindex................long long int (8 bytes)
        validationEpsilon.....1e-12
    [compilation]
        isOmpCompiled...............1
        isMpiCompiled...............1
        isMpiSubCommCompiled........0
        isGpuCompiled...............0
        isHipCompiled...............0
        isCuQuantumCompiled.........0
        isCheckpointingCompiled.....0
    [deployment]
        isOmpEnabled...........1
        isMpiEnabled...........1
        isGpuEnabled...........0
        isCuQuantumEnabled.....0
    [cpu]
        numCpuCores.......14 per machine
        numOmpProcs.......14 per machine
        numOmpThrds.......14 per node
        cpuMemory.........36 GiB per machine
        cpuMemoryFree.....unknown
    [gpu]
        numGpus................N/A
        gpuDirect..............N/A
        gpuMemPools............N/A
        gpuMemory..............N/A
        gpuMemoryFree..........N/A
        gpuCache...............N/A
        numThreadsPerBlock.....N/A
    [distribution]
        isMpiUserOwned..........0
        isMpiGpuAware...........0
        isGpuSharingEnabled.....N/A
        numMpiNodes.............16
    [statevector limits]
        minQubitsForMpi.............4
        maxQubitsForCpu.............31
        maxQubitsForGpu.............N/A
        maxQubitsForMpiCpu..........34
        maxQubitsForMpiGpu..........N/A
        maxQubitsForMemOverflow.....58
        maxQubitsForIndOverflow.....63
    [density matrix limits]
        minQubitsForMpi.............4
        maxQubitsForCpu.............15
        maxQubitsForGpu.............N/A
        maxQubitsForMpiCpu..........19
        maxQubitsForMpiGpu..........N/A
        maxQubitsForMemOverflow.....28
        maxQubitsForIndOverflow.....31
    [statevector autodeployment]
        8 qubits......[omp] 
        30 qubits.....[omp] [mpi] 
    [density matrix autodeployment]
        4 qubits......[omp] 
        15 qubits.....[omp] [mpi]
 * ```
 * 
 * @throws @validationerror
 * - if the QuEST environment is not initialised.
 * @see
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_environments.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_environments.cpp) examples
 * @author Tyson Jones
 */
void reportQuESTEnv();

/** Indicates whether the QuEST execution environment is currently initialised.
 * 
 * Unlike other QuEST functions, this can be called at any time, including before
 * QuEST initialisation, and after finalisation.
 *
 * @returns @c 1 if the environment is initialised, otherwise @c 0.
 * @author Tyson Jones
 */
int isQuESTEnvInit();

/** Returns a copy of the active QuEST execution environment.
 *
 * The returned QuESTEnv describes the active deployment and MPI rank information.
 * This can be useful for making programmatical decisions based on the environment.
 * 
 * @myexample
 * 
 * ```cpp
    QuESTEnv env = getQuESTEnv();

    if (env.isDistributed && env.isGpuAccelerated && ! env.isMpiGpuAware)
        printf("What a waste!\n");
 * ```
 *
 * @returns A copy of the active QuESTEnv.
 * @throws @validationerror
 * - if the QuEST environment is not initialised.
 * @author Tyson Jones
 */
QuESTEnv getQuESTEnv();



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // ENVIRONMENT_H

/** @} */ // (end file-wide doxygen defgroup)
