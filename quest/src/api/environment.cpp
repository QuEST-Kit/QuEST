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
#include "quest/src/core/parser.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/envvars.hpp"
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


static QuESTEnv* global_envPtr = nullptr;



/*
 * PRIVATE QUESTENV INITIALISATION HISTORY
 *
 * indicating whether QuEST has ever been finalized. This is important, since 
 * the QuEST environment can only ever be initialised once, and can never
 * be re-initialised after finalisation, due to re-initialisation of MPI 
 * being undefined behaviour.
 */


static bool global_hasEnvBeenFinalized = false;



/*
 * PRIVATE QUESTENV INITIALISATION INNER FUNCTIONS
 */


void validateAndInitCustomQuESTEnv(int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread, const char* caller) {

    // ensure that we are never re-initialising QuEST (even after finalize) because
    // this leads to undefined behaviour in distributed mode, as per the MPI std,
    // regardless of whether the user owns MPI
    validate_envNeverInit(global_envPtr != nullptr, global_hasEnvBeenFinalized, caller);

    // load env-vars before validating deployment mode, because some env vars can
    // affect validation (such as QUEST_PERMIT_NODES_TO_SHARE_GPU). note that
    // some env-vars (like QUEST_DEFAULT_NUM_GPU_THREADS_PER_BLOCK) will be here
    // validated to have a correct format (like an int), but the validity of its
    // actual value will be checked later (since it requires deciding GPU-accel).
    envvars_validateAndLoadEnvVars(caller);
    validateconfig_setEpsilonToDefault();

    // ensure the chosen deployment is compiled and supported by hardware.
    // note that these error messages will be printed by every node because
    // validation occurs before comm_init() below, so all processes spawned
    // by mpirun believe they are each the main rank. This seems unavoidable.
    validate_newEnvDeploymentMode(useDistrib, useGpuAccel, useMultithread, caller);

    // overwrite deployments (left as modeflag::USE_AUTO=-1) with 0,1 (a bool),
    // which crucially, resolves useDistrib, permitting its consultation below
    autodep_chooseQuESTEnvDeployment(useDistrib, useGpuAccel, useMultithread);

    // ensure that current state of MPI is valid
    validate_mpiInitStatus(useDistrib, userOwnsMpi, caller);

    // optionally initialise MPI; necessary before completing validation,
    // and before any GPU initialisation and validation, since we will
    // perform that specifically upon the MPI-process-bound GPU(s). Further,
    // we can make sure validation errors are reported only by the root node.
    if (useDistrib)
        comm_init(userOwnsMpi);

    validate_newEnvDistributedBetweenPower2Nodes(caller);

    /// @todo
    /// consider immediately disabling MPI here if comm_numNodes() == 1
    /// (also overwriting useDistrib = 0)

    // bind MPI nodes to unique GPUs; even when not distributed,
    // and before we have validated local GPUs are compatible
    if (useGpuAccel)
        gpu_bindLocalGPUsToNodes();

    // consult environment variable to decide whether to allow GPU sharing 
    // (default = false) which informs whether below validation is triggered
    bool permitGpuSharing = envvars_getWhetherGpuSharingIsPermitted();

    // each MPI process should ordinarily use a unique GPU. This is 
    // critical when initializing cuQuantum so that we don't re-init 
    // cuStateVec on any paticular GPU (which can apparently cause a
    // so-far-unwitnessed runtime error), but is otherwise essential
    // for good performance. GPU sharing is useful for unit testing
    // however permitting a single GPU to test CUDA+MPI deployment
    if (useGpuAccel && useDistrib && ! permitGpuSharing)
        validate_newEnvNodesEachHaveUniqueGpu(caller);

    /// @todo
    /// should we warn here if each machine contains
    /// more GPUs than deployed MPI-processes (some GPUs idle)?

    // validate the initial numTPB env-var (if specified) is valid
    int initNumThreadsPerBlock = envvars_getDefaultNumGpuThreadsPerBlock();
    validate_numGpuThreadsPerBlock(initNumThreadsPerBlock, useGpuAccel, caller);
    gpu_setNumThreadsPerBlock(initNumThreadsPerBlock);

    // cuQuantum is always used in GPU-accelerated envs when available
    bool useCuQuantum = useGpuAccel && gpu_isCuQuantumCompiled();
    if (useCuQuantum) {
        validate_gpuIsCuQuantumCompatible(caller); // assesses above bound GPU
        gpu_initCuQuantum();
    }

    // MPI GPU-awareness detection is platform specific; sometimes it is
    // known at compile-time, other times according to env-vars
    bool isMpiGpuAware = comm_isMpiGpuAware();

    // initialise RNG, used by measurements and random-state generation
    rand_setSeedsToDefault();

    // allocate space for the global QuESTEnv singleton (overwriting nullptr, unless malloc fails)
    global_envPtr = (QuESTEnv*) malloc(sizeof(QuESTEnv));

    // pedantically check that teeny tiny malloc just succeeded
    if (global_envPtr == nullptr)
        error_allocOfQuESTEnvFailed();

    // bind deployment info to global instance (autocasting int to bool)
    global_envPtr->isMultithreaded     = useMultithread;
    global_envPtr->isGpuAccelerated    = useGpuAccel;
    global_envPtr->isDistributed       = useDistrib;
    global_envPtr->isMpiUserOwned      = userOwnsMpi;
    global_envPtr->isMpiGpuAware       = isMpiGpuAware;
    global_envPtr->isCuQuantumEnabled  = useCuQuantum;
    global_envPtr->isGpuSharingEnabled = permitGpuSharing;

    // bind distributed info
    global_envPtr->rank     = (useDistrib)? comm_getRank()     : 0;
    global_envPtr->numNodes = (useDistrib)? comm_getNumNodes() : 1;
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
        {"isOmpCompiled",         cpu_isOpenmpCompiled()},
        {"isMpiCompiled",         comm_isMpiCompiled()},
        {"isMpiSubCommCompiled",  comm_isMpiSubCommCompiled()},
        {"isGpuCompiled",         gpu_isGpuCompiled()},
        {"isHipCompiled",         gpu_isHipCompiled()},
        {"isCuQuantumCompiled",   gpu_isCuQuantumCompiled()},
    });
}


void printDeploymentInfo() {

    print_table(
        "deployment", {
        {"isOmpEnabled",        global_envPtr->isMultithreaded},
        {"isMpiEnabled",        global_envPtr->isDistributed},
        {"isGpuEnabled",        global_envPtr->isGpuAccelerated},
        {"isCuQuantumEnabled",  global_envPtr->isCuQuantumEnabled},
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
        {"numOmpThrds",   (cpu_isOpenmpCompiled())? printer_toStr(cpu_getAvailableNumThreads()) + pn : na},
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
        {"numThreadsPerBlock", isGpu? printer_toStr(gpu_getNumThreadsPerBlock()) : na},
    });
}


void printDistributionInfo() {

    using namespace printer_substrings;

    bool comm = global_envPtr->isDistributed;
    bool gpu  = global_envPtr->isGpuAccelerated;
    bool both = comm && gpu;

    print_table(
        "distribution", {
        {"isMpiUserOwned",      comm? printer_toStr(global_envPtr->isMpiUserOwned) : na},
        {"isMpiGpuAware",       comm? printer_toStr(global_envPtr->isMpiGpuAware ) : na},
        {"isGpuSharingEnabled", both? printer_toStr(global_envPtr->isGpuSharingEnabled) : na},
        {"numMpiNodes",         printer_toStr(global_envPtr->numNodes)},
    });
}


void printQuregSizeLimits(bool isDensMatr) {

    using namespace printer_substrings;

    // for brevity
    int numNodes = global_envPtr->numNodes;

    // by default, CPU limits are unknown (because memory query might fail)
    string maxQbForCpu = un;
    string maxQbForMpiCpu = un;

    // max CPU registers are only determinable if RAM query succeeds
    try {
        qindex cpuMem = mem_tryGetLocalRamCapacityInBytes();
        maxQbForCpu = printer_toStr(mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, 1, cpuMem));

        // and the max MPI sizes are only relevant when env is distributed
        if (global_envPtr->isDistributed)
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
    if (global_envPtr->isGpuAccelerated) {
        qindex gpuMem = gpu_getCurrentAvailableMemoryInBytes();
        maxQbForGpu = printer_toStr(mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, 1, gpuMem));

        // and the max MPI sizes are further only relevant when env is distributed 
        if (global_envPtr->isDistributed)
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
    int maxQubits = mem_getMaxNumQuregQubitsBeforeGlobalMemSizeofOverflow(isDensMatr, global_envPtr->numNodes);

    for (int numQubits=1; numQubits<maxQubits; numQubits++) {

        // re-choose auto deployment
        useDistrib  = modeflag::USE_AUTO;
        useGpuAccel = modeflag::USE_AUTO;
        useMulti    = modeflag::USE_AUTO;;
        autodep_chooseQuregDeployment(numQubits, isDensMatr, useDistrib, useGpuAccel, useMulti, *global_envPtr);

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

    const bool userOwnsMpi = false;
    validateAndInitCustomQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread, __func__);
}


void initQuESTEnv() {

    const bool userOwnsMpi = false;
    validateAndInitCustomQuESTEnv(modeflag::USE_AUTO, userOwnsMpi, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);
}


int isQuESTEnvInit() {

    return (int) (global_envPtr != nullptr);
}


QuESTEnv getQuESTEnv() {
    validate_envIsInit(__func__);

    // returns a copy, so cheeky users calling memcpy() upon const struct still won't mutate
    return *global_envPtr;
}


void finalizeQuESTEnv() {
    validate_envIsInit(__func__);

    // NOTE:
    // calling this will not automatically
    // free the memory of existing Quregs

    if (global_envPtr->isGpuAccelerated)
        gpu_clearCache(); // syncs first

    if (global_envPtr->isGpuAccelerated && gpu_isCuQuantumCompiled())
        gpu_finalizeCuQuantum();

    if (global_envPtr->isDistributed) {
        comm_sync();
        comm_end();
    }

    // free global env's heap memory and flag it as unallocated
    free(global_envPtr);
    global_envPtr = nullptr;

    // flag that the environment was finalised, to ensure it is never re-initialised
    global_hasEnvBeenFinalized = true;
}


void syncQuESTEnv() {
    validate_envIsInit(__func__);

    if (global_envPtr->isGpuAccelerated)
        gpu_sync();

    if (global_envPtr->isDistributed)
        comm_sync();
}


void reportQuESTEnv() {
    validate_envIsInit(__func__);
    validate_numReportedNewlinesAboveZero(__func__); // because trailing newline mandatory

    /// @todo add function to write this output to file (useful for HPC debugging)

    printer_sync();

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

    printer_sync();
}


void getQuESTEnvironmentString(char str[200]) {
    validate_envIsInit(__func__);

    int numThreads = cpu_isOpenmpCompiled()? cpu_getAvailableNumThreads() : 1;
    int cuQuantum = global_envPtr->isGpuAccelerated && gpu_isCuQuantumCompiled();
    int gpuDirect = global_envPtr->isGpuAccelerated && gpu_isDirectGpuCommPossible();

    snprintf(str, 200, "CUDA=%d OpenMP=%d MPI=%d userOwnsMPI=%d threads=%d ranks=%d cuQuantum=%d gpuDirect=%d",
        global_envPtr->isGpuAccelerated,
        global_envPtr->isMultithreaded,
        global_envPtr->isDistributed,
        global_envPtr->isMpiUserOwned,
        numThreads,
        global_envPtr->numNodes,
        cuQuantum,
        gpuDirect);
}


// end de-mangler
}
