/** @file
 * API definitions for managing QuESTEnv instances, which
 * themselves control and query the deployment environment.
 * 
 * @author Tyson Jones 
 */

#include "quest/include/environment.h"
#include "quest/include/precision.h"
#include "quest/include/modes.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/core/randomiser.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <typeinfo>
#include <cstring>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>
#include <tuple>

using std::string;



/*
 * PRIVATE QUESTENV SINGLETON
 *
 * Global to this file, accessible to other files only through 
 * getQuESTEnv() which returns a copy, which also has const fields.
 * The use of static ensures we never accidentally expose the "true"
 * runtime single instance to other files. We allocate the env
 * in heap memory (hence the pointer) so that we can defer 
 * initialisation of the const fields. The address being nullptr
 * indicates the QuESTEnv is not currently initialised; perhaps never,
 * or it was but has since been finalized.
 */


static QuESTEnv* globalEnvPtr = nullptr;



/*
 * PRIVATE QUESTENV INITIALISATION HISTORY
 *
 * indicating whether QuEST has ever been finalized. This is important, since 
 * the QuEST environment can only ever be initialised once, and can never
 * be re-initialised after finalisation, due to re-initialisation of MPI 
 * being undefined behaviour.
 */


static bool hasEnvBeenFinalized = false;



/*
 * PRIVATE QUESTENV INITIALISATION INNER FUNCTIONS
 */


void validateAndInitCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread, const char* caller) {

    // ensure that we are never re-initialising QuEST (even after finalize) because
    // this leads to undefined behaviour in distributed mode, as per the MPI
    validate_envNeverInit(globalEnvPtr != nullptr, hasEnvBeenFinalized, caller);

    // ensure the chosen deployment is compiled and supported by hardware.
    // note that these error messages will be printed by every node because
    // validation occurs before comm_init() below, so all processes spawned
    // by mpirun believe they are each the main rank. This seems unavoidable.
    validate_newEnvDeploymentMode(useDistrib, useGpuAccel, useMultithread, caller);

    // overwrite deployments left as modeflag::USE_AUTO
    autodep_chooseQuESTEnvDeployment(useDistrib, useGpuAccel, useMultithread);

    // optionally initialise MPI; necessary before completing validation,
    // and before any GPU initialisation and validation, since we will
    // perform that specifically upon the MPI-process-bound GPU(s). Further,
    // we can make sure validation errors are reported only by the root node.
    if (useDistrib)
        comm_init();

    validate_newEnvDistributedBetweenPower2Nodes(caller);

    /// @todo
    /// consider immediately disabling MPI here if comm_numNodes() == 1
    /// (also overwriting useDistrib = 0)

    // bind MPI nodes to unique GPUs; even when not distributed,
    // and before we have validated local GPUs are compatible
    if (useGpuAccel)
        gpu_bindLocalGPUsToNodes();

    // each MPI process must use a unique GPU. This is critical when
    // initializing cuQuantum, so we don't re-init cuStateVec on any
    // paticular GPU (causing runtime error), but still ensures we 
    // keep good performance in our custom backend GPU code; there is
    // no reason to use multi-nodes-per-GPU except for dev/debugging.
    if (useGpuAccel && useDistrib && ! PERMIT_NODES_TO_SHARE_GPU)
        validate_newEnvNodesEachHaveUniqueGpu(caller);

    /// @todo
    /// should we warn here if each machine contains
    /// more GPUs than deployed MPI-processes (some GPUs idle)?

    // use cuQuantum if compiled
    if (useGpuAccel && gpu_isCuQuantumCompiled()) {
        validate_gpuIsCuQuantumCompatible(caller); // assesses above bound GPU
        gpu_initCuQuantum();
    }

    // initialise RNG, used by measurements and random-state generation
    rand_setSeedsToDefault();

    // allocate space for the global QuESTEnv singleton (overwriting nullptr, unless malloc fails)
    globalEnvPtr = (QuESTEnv*) malloc(sizeof(QuESTEnv));

    // pedantically check that teeny tiny malloc just succeeded
    if (globalEnvPtr == nullptr)
        error_allocOfQuESTEnvFailed();

    /// @todo the below memcpy is naughty (QuESTEnv has no trivial copy-assignment) and causes compiler warning. Fix!

    // initialise it to a local env
    QuESTEnv env = {

        // bind deployment info
        .isMultithreaded  = useMultithread,
        .isGpuAccelerated = useGpuAccel,
        .isDistributed    = useDistrib,

        // set distributed info
        .rank     = (useDistrib)? comm_getRank()     : 0,
        .numNodes = (useDistrib)? comm_getNumNodes() : 1,
    };
    memcpy(globalEnvPtr, &env, sizeof(QuESTEnv));
}



/*
 * PRIVATE QUESTENV REPORTING INNER FUNCTIONS
 */


void printPrecisionInfo() {

    /// @todo
    /// - report MPI qcomp type?
    /// - report CUDA qcomp type?
    /// - report CUDA kernel qcomp type?

    print_table(
        "precision", {
        {"qreal",  printer_getQrealType()  + " (" + printer_getMemoryWithUnitStr(sizeof(qreal)) + ")"},

        /// @todo this is showing the backend C++ qcomp type, rather than that actually wieldable
        /// by the user which could the C-type. No idea how to solve this however!
        {"qcomp",  printer_getQcompType()  + " (" + printer_getMemoryWithUnitStr(sizeof(qcomp)) + ")"},

        {"qindex", printer_getQindexType() + " (" + printer_getMemoryWithUnitStr(sizeof(qindex)) + ")"},

        /// @todo this currently prints 0 when epsilon is inf (encoded by zero), i.e. disabled
        {"validationEpsilon", printer_toStr(validateconfig_getEpsilon())},
    });
}


void printCompilationInfo() {

    print_table(
        "compilation", {
        {"isMpiCompiled",      comm_isMpiCompiled()},
        {"isGpuCompiled",       gpu_isGpuCompiled()},
        {"isOmpCompiled",       cpu_isOpenmpCompiled()},
        {"isCuQuantumCompiled", gpu_isCuQuantumCompiled()},
    });
}


void printDeploymentInfo() {

    print_table(
        "deployment", {
        {"isMpiEnabled", globalEnvPtr->isDistributed},
        {"isGpuEnabled", globalEnvPtr->isGpuAccelerated},
        {"isOmpEnabled", globalEnvPtr->isMultithreaded},
    });
}


void printCpuInfo() {

    using namespace printer_substrings;

    // assume RAM is unknown unless it can be queried
    string ram = un;
    try { 
        ram = printer_getMemoryWithUnitStr(mem_tryGetLocalRamCapacityInBytes()) + pm; 
    } catch(mem::COULD_NOT_QUERY_RAM e){};

    /// @todo
    /// - CPU info e.g. speeds/caches?

    print_table(
        "cpu", {
        {"numCpuCores",   printer_toStr(std::thread::hardware_concurrency()) + pm},
        {"numOmpProcs",   (cpu_isOpenmpCompiled())? printer_toStr(cpu_getNumOpenmpProcessors()) + pm : na},
        {"numOmpThrds",   (cpu_isOpenmpCompiled())? printer_toStr(cpu_getCurrentNumThreads()) + pn : na},
        {"cpuMemory",     ram},
        {"cpuMemoryFree", un},
    });
}


void printGpuInfo() {

    using namespace printer_substrings;

    /// @todo below:
    /// - GPU compute capability
    /// - GPU #SVMs etc

    // must not query any GPU facilities unless confirmed compiled and available
    bool isComp = gpu_isGpuCompiled();
    bool isGpu = isComp && gpu_isGpuAvailable();

    print_table(
        "gpu", {
        {"numGpus",       isComp? printer_toStr(gpu_getNumberOfLocalGpus())    : na},
        {"gpuDirect",     isGpu?  printer_toStr(gpu_isDirectGpuCommPossible()) : na},
        {"gpuMemPools",   isGpu?  printer_toStr(gpu_doesGpuSupportMemPools())  : na},
        {"gpuMemory",     isGpu?  printer_getMemoryWithUnitStr(gpu_getTotalMemoryInBytes())            + pg : na},
        {"gpuMemoryFree", isGpu?  printer_getMemoryWithUnitStr(gpu_getCurrentAvailableMemoryInBytes()) + pg : na},
        {"gpuCache",      isGpu?  printer_getMemoryWithUnitStr(gpu_getCacheMemoryInBytes())            + pg : na},
    });
}


void printDistributionInfo() {

    using namespace printer_substrings;

    print_table(
        "distribution", {
        {"isMpiGpuAware", (comm_isMpiCompiled())? printer_toStr(comm_isMpiGpuAware()) : na},
        {"numMpiNodes",   printer_toStr(globalEnvPtr->numNodes)},
    });
}


void printQuregSizeLimits(bool isDensMatr) {

    using namespace printer_substrings;

    // for brevity
    int numNodes = globalEnvPtr->numNodes;

    // by default, CPU limits are unknown (because memory query might fail)
    string maxQbForCpu = un;
    string maxQbForMpiCpu = un;

    // max CPU registers are only determinable if RAM query succeeds
    try {
        qindex cpuMem = mem_tryGetLocalRamCapacityInBytes();
        maxQbForCpu = printer_toStr(mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, 1, cpuMem));

        // and the max MPI sizes are only relevant when env is distributed
        if (globalEnvPtr->isDistributed)
            maxQbForMpiCpu = printer_toStr(mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, numNodes, cpuMem));

        // when MPI irrelevant, change their status from "unknown" to "N/A"
        else
            maxQbForMpiCpu = na;

    // no problem if we can't query RAM; we simply don't report relevant limits
    } catch(mem::COULD_NOT_QUERY_RAM e) {};

    // GPU limits are default N/A because they're always determinable when relevant
    string maxQbForGpu = na;
    string maxQbForMpiGpu = na;

    // max GPU registers only relevant if env is GPU-accelerated
    if (globalEnvPtr->isGpuAccelerated) {
        qindex gpuMem = gpu_getCurrentAvailableMemoryInBytes();
        maxQbForGpu = printer_toStr(mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, 1, gpuMem));

        // and the max MPI sizes are further only relevant when env is distributed 
        if (globalEnvPtr->isDistributed)
            maxQbForMpiGpu = printer_toStr(mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, numNodes, gpuMem));
    }

    // tailor table title to type of Qureg
    string prefix = (isDensMatr)? "density matrix" : "statevector";
    string title = prefix + " limits";

    print_table(
        title, {
        {"minQubitsForMpi",     (numNodes>1)? printer_toStr(mem_getMinNumQubitsForDistribution(numNodes)) : na},
        {"maxQubitsForCpu",     maxQbForCpu},
        {"maxQubitsForGpu",     maxQbForGpu},
        {"maxQubitsForMpiCpu",  maxQbForMpiCpu},
        {"maxQubitsForMpiGpu",  maxQbForMpiGpu},
        {"maxQubitsForMemOverflow", printer_toStr(mem_getMaxNumQuregQubitsBeforeGlobalMemSizeofOverflow(isDensMatr, numNodes))},
        {"maxQubitsForIndOverflow", printer_toStr(mem_getMaxNumQuregQubitsBeforeIndexOverflow(isDensMatr))},
    });
}


void printQuregAutoDeployments(bool isDensMatr) {

    // build all table rows dynamically before print
    std::vector<std::tuple<string, string>> rows;

    // we will get auto-deployment for every possible number of qubits; silly but cheap and robust!
    int useDistrib,  useGpuAccel,  useMulti;
    int prevDistrib, prevGpuAccel, prevMulti;

    // assume all deployments disabled for 1 qubit
    prevDistrib  = 0;
    prevGpuAccel = 0;
    prevMulti    = 0;

    // test to theoretically max #qubits, surpassing max that can fit in RAM and GPUs, because
    // auto-deploy will still try to deploy there to (then subsequent validation will fail)
    int maxQubits = mem_getMaxNumQuregQubitsBeforeGlobalMemSizeofOverflow(isDensMatr, globalEnvPtr->numNodes);

    for (int numQubits=1; numQubits<maxQubits; numQubits++) {

        // re-choose auto deployment
        useDistrib  = modeflag::USE_AUTO;
        useGpuAccel = modeflag::USE_AUTO;
        useMulti    = modeflag::USE_AUTO;;
        autodep_chooseQuregDeployment(numQubits, isDensMatr, useDistrib, useGpuAccel, useMulti, *globalEnvPtr);

        // skip if deployments are unchanged
        if (useDistrib  == prevDistrib  &&
            useGpuAccel == prevGpuAccel &&
            useMulti    == prevMulti)
            continue; 

        // else prepare string summarising the new deployments (trailing space is fine)
        string value = "";
        if (useMulti)
            value += "[omp] "; // ordered by #qubits to attempt consistent printed columns
        if (useGpuAccel)
            value += "[gpu] ";
        if (useDistrib)
            value += "[mpi] ";

        // log the #qubits of the deployment change
        rows.push_back({printer_toStr(numQubits) + " qubits", value});

        // skip subsequent qubits with the same deployments
        prevDistrib  = useDistrib;
        prevGpuAccel = useGpuAccel;
        prevMulti    = useMulti;
    }

    // tailor table title to type of Qureg
    string prefix = (isDensMatr)? "density matrix" : "statevector";
    string title = prefix + " autodeployment";
    rows.empty()?
        print_table(title, "(no parallelisations available)"):
        print_table(title, rows);
}



/*
 * API FUNCTIONS
 */


// enable invocation by both C and C++ binaries
extern "C" {


void initCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread) {

    validateAndInitCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread, __func__);
}


void initQuESTEnv() {

    validateAndInitCustomQuESTEnv(modeflag::USE_AUTO, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);
}


int isQuESTEnvInit() {

    return (int) (globalEnvPtr != nullptr);
}


QuESTEnv getQuESTEnv() {
    validate_envIsInit(__func__);

    // returns a copy, so cheeky users calling memcpy() upon const struct still won't mutate
    return *globalEnvPtr;
}


void finalizeQuESTEnv() {
    validate_envIsInit(__func__);

    // NOTE:
    // calling this will not automatically
    // free the memory of existing Quregs

    if (globalEnvPtr->isGpuAccelerated)
        gpu_clearCache(); // syncs first

    if (globalEnvPtr->isGpuAccelerated && gpu_isCuQuantumCompiled())
        gpu_finalizeCuQuantum();

    if (globalEnvPtr->isDistributed) {
        comm_sync();
        comm_end();
    }

    // free global env's heap memory and flag it as unallocated
    free(globalEnvPtr);
    globalEnvPtr = nullptr;

    // flag that the environment was finalised, to ensure it is never re-initialised
    hasEnvBeenFinalized = true;
}


void syncQuESTEnv() {
    validate_envIsInit(__func__);

    if (globalEnvPtr->isGpuAccelerated)
        gpu_sync();

    if (globalEnvPtr->isDistributed)
        comm_sync();
}


void reportQuESTEnv() {
    validate_envIsInit(__func__);
    validate_numReportedNewlinesAboveZero(__func__); // because trailing newline mandatory

    /// @todo add function to write this output to file (useful for HPC debugging)

    print_label("QuEST execution environment");

    bool statevec = false;
    bool densmatr = true;

    // we attempt to report properties of available hardware facilities
    // (e.g. number of CPU cores, number of GPUs) even if the environment is not
    // making use of them, to inform the user how they might change deployment.
    printPrecisionInfo();
    printCompilationInfo();
    printDeploymentInfo();
    printCpuInfo();
    printGpuInfo();
    printDistributionInfo();
    printQuregSizeLimits(statevec);
    printQuregSizeLimits(densmatr);
    printQuregAutoDeployments(statevec);
    printQuregAutoDeployments(densmatr);

    // exclude mandatory newline above
    print_oneFewerNewlines();
}


void getEnvironmentString(char str[200]) {
    validate_envIsInit(__func__);

    QuESTEnv env = getQuESTEnv();

    int numThreads = cpu_isOpenmpCompiled()? cpu_getCurrentNumThreads() : 1;
    int cuQuantum = env.isGpuAccelerated && gpu_isCuQuantumCompiled();
    int gpuDirect = env.isGpuAccelerated && gpu_isDirectGpuCommPossible();

    snprintf(str, 200, "CUDA=%d OpenMP=%d MPI=%d threads=%d ranks=%d cuQuantum=%d gpuDirect=%d",
        env.isGpuAccelerated,
        env.isMultithreaded,
        env.isDistributed,
        numThreads,
        env.numNodes,
        cuQuantum,
        gpuDirect);
}


// end de-mangler
}
