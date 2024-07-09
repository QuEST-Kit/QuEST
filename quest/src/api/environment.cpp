/** @file
 * API definitions for managing QuESTEnv instances, which
 * themselves control and query the deployment environment. 
 */

#include "environment.h"
#include "precision.h"
#include "modes.h"

#include "../core/memory.hpp"
#include "../core/formatter.hpp"
#include "../core/autodeployer.hpp"
#include "../core/validation.hpp"
#include "../comm/comm_config.hpp"
#include "../cpu/cpu_config.hpp"

#include <iostream>
#include <typeinfo>
#include <cstring>
#include <string>
#include <thread>
#include <vector>
#include <tuple>

// provides substrings (by, na, pm, etc) used by reportQuESTEnv
using namespace form_substrings;



/*
 * PRIVATE QUESTENV SINGLETON
 *
 * Global to this file, accessible to other files only through 
 * getQuESTEnv() which returns a copy, which also has const fields.
 * The use of static ensures we never accidentally expose the "true"
 * runtime single instance to other files. We allocate the env
 * in heap memory (hence the pointer) so that we can defer 
 * initialisation of the const fields. The address being NULL
 * indicates the QuESTEnv is not currently initialised; perhaps never,
 * or it was but has since been finalized.
 */


static QuESTEnv* globalEnvPtr = NULL;



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
    validate_envNeverInit(globalEnvPtr != NULL, hasEnvBeenFinalized, caller);

    // ensure the chosen deployment is compiled and supported by hardware.
    // note that these error messages will be printed by every node because
    // validation occurs before comm_init() below, so all processes spawned
    // by mpirun believe they are each the main rank. This seems unavoidable.
    validate_newEnvDeploymentMode(useDistrib, useGpuAccel, useMultithread, caller);

    // overwrite deployments left as modeflag::USE_AUTO
    autodep_chooseQuESTEnvDeployment(useDistrib, useGpuAccel, useMultithread);

    // use of cuQuantum requires a modern GPU
    if (useGpuAccel && gpu_isCuQuantumCompiled())
        validate_gpuIsCuQuantumCompatible(caller);

    // optionally initialise MPI; necessary before completing validation
    if (useDistrib)
        comm_init();

    validate_newEnvDistributedBetweenPower2Nodes(caller);

    if (useGpuAccel && useDistrib) {

        // TODO:
        // validate 2^N local GPUs
    }

    // make a new, local env
    QuESTEnv env = {

        // bind deployment info
        .isMultithreaded  = useMultithread,
        .isGpuAccelerated = useGpuAccel,
        .isDistributed    = useDistrib,

        // set distributed info
        .rank     = (useDistrib)? comm_getRank()     : 0,
        .numNodes = (useDistrib)? comm_getNumNodes() : 1,
    };

    // in multi-GPU settings, bind each MPI process to one GPU
    if (useGpuAccel && useDistrib)
        gpu_bindLocalGPUsToNodes(env.rank);

    // in GPU settings, if cuQuantum is being used, initialise it
    if (useGpuAccel && gpu_isCuQuantumCompiled())
        gpu_initCuQuantum();

    // TODO: setup RNG

    // allocate space for the global QuESTEnv singleton (overwriting NULL, unless malloc fails)
    globalEnvPtr = (QuESTEnv*) malloc(sizeof(QuESTEnv));

    // pedantically check that teeny tiny malloc just succeeded
    if (globalEnvPtr == NULL)
        error_allocOfQuESTEnvFailed();

    // initialise it to our local env
    memcpy(globalEnvPtr, &env, sizeof(QuESTEnv));
}



/*
 * PRIVATE QUESTENV REPORTING INNER FUNCTIONS
 */


void printPrecisionInfo() {

    // TODO
    // - report MPI qcomp type?
    // - report CUDA qcomp type?
    // - report CUDA kernel qcomp type?

    form_printTable(
        "precision", {
        {"qreal",  form_getQrealType()  + " (" + form_str(sizeof(qreal))  + by + ")"},
        {"qcomp",  form_getQcompType()  + " (" + form_str(sizeof(qcomp))  + by + ")"},
        {"qindex", form_getQindexType() + " (" + form_str(sizeof(qindex)) + by + ")"},
        {"validationEpsilon", form_str(VALIDATION_EPSILON)},
    });
}


void printCompilationInfo() {

    form_printTable(
        "compilation", {
        {"isMpiCompiled",      comm_isMpiCompiled()},
        {"isGpuCompiled",       gpu_isGpuCompiled()},
        {"isOmpCompiled",       cpu_isOpenmpCompiled()},
        {"isCuQuantumCompiled", gpu_isCuQuantumCompiled()},
    });
}


void printDeploymentInfo() {

    form_printTable(
        "deployment", {
        {"isMpiEnabled", globalEnvPtr->isDistributed},
        {"isGpuEnabled", globalEnvPtr->isGpuAccelerated},
        {"isOmpEnabled", globalEnvPtr->isMultithreaded},
    });
}


void printCpuInfo() {

    // assume RAM is unknown unless it can be queried
    std::string ram = un;
    try { 
        ram = form_str(mem_tryGetLocalRamCapacityInBytes()) + by + pm; 
    } catch(mem::COULD_NOT_QUERY_RAM e){};

    // TODO
    // - CPU info e.g. speeds/caches?

    form_printTable(
        "cpu", {
        {"numCpuCores",   form_str(std::thread::hardware_concurrency()) + pm},
        {"numOmpProcs",   (cpu_isOpenmpCompiled())? form_str(cpu_getNumOpenmpProcessors()) + pm : na},
        {"numOmpThrds",   (cpu_isOpenmpCompiled())? form_str(cpu_getCurrentNumThreads()) + pn : na},
        {"cpuMemory",     ram},
        {"cpuMemoryFree", un},
    });
}


void printGpuInfo() {

    // TODO below:
    // - GPU compute capability
    // - GPU #SVMs etc

    form_printTable(
        "gpu", {
        {"numGpus",       (gpu_isGpuCompiled())? form_str(gpu_getNumberOfLocalGpus()) : un},
        {"gpuDirect",     (gpu_isGpuCompiled())? form_str(gpu_isDirectGpuCommPossible()) : na},
        {"gpuMemPools",   (gpu_isGpuCompiled())? form_str(gpu_doesGpuSupportMemPools()) : na},
        {"gpuMemory",     (gpu_isGpuCompiled())? form_str(gpu_getTotalMemoryInBytes()) + by + pg : na},
        {"gpuMemoryFree", (gpu_isGpuCompiled())? form_str(gpu_getTotalMemoryInBytes()) + by + pg : na},
    });
}


void printDistributionInfo() {

    form_printTable(
        "distribution", {
        {"isMpiGpuAware", (comm_isMpiCompiled())? form_str(comm_isMpiGpuAware()) : na},
        {"numMpiNodes",   form_str(globalEnvPtr->numNodes)},
    });
}


void printQuregSizeLimits(bool isDensMatr) {

    // for brevity
    int numNodes = globalEnvPtr->numNodes;

    // by default, CPU limits are unknown (because memory query might fail)
    std::string maxQbForCpu = un;
    std::string maxQbForMpiCpu = un;

    // max CPU registers are only determinable if RAM query succeeds
    try {
        qindex cpuMem = mem_tryGetLocalRamCapacityInBytes();
        maxQbForCpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, 1, cpuMem));

        // and the max MPI sizes are only relevant when env is distributed
        if (globalEnvPtr->isDistributed)
            maxQbForMpiCpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, numNodes, cpuMem));

        // when MPI irrelevant, change their status from "unknown" to "N/A"
        else
            maxQbForMpiCpu = na;

    // no problem if we can't query RAM; we simply don't report relevant limits
    } catch(mem::COULD_NOT_QUERY_RAM e) {};

    // GPU limits are default N/A because they're always determinable when relevant
    std::string maxQbForGpu = na;
    std::string maxQbForMpiGpu = na;

    // max GPU registers only relevant if env is GPU-accelerated
    if (globalEnvPtr->isGpuAccelerated) {
        qindex gpuMem = gpu_getCurrentAvailableMemoryInBytes();
        maxQbForGpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, 1, gpuMem));

        // and the max MPI sizes are further only relevant when env is distributed 
        if (globalEnvPtr->isDistributed)
            maxQbForMpiGpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, numNodes, gpuMem));
    }

    // tailor table title to type of Qureg
    std::string prefix = (isDensMatr)? "density matrix" : "statevector";
    std::string title = prefix + " limits";

    form_printTable(
        title, {
        {"minQubitsForMpi",     (numNodes>1)? form_str(mem_getMinNumQubitsForDistribution(numNodes)) : na},
        {"maxQubitsForCpu",     maxQbForCpu},
        {"maxQubitsForGpu",     maxQbForGpu},
        {"maxQubitsForMpiCpu",  maxQbForMpiCpu},
        {"maxQubitsForMpiGpu",  maxQbForMpiGpu},
        {"maxQubitsForMemOverflow", form_str(mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, numNodes))},
        {"maxQubitsForIndOverflow", form_str(mem_getMaxNumQubitsBeforeIndexOverflow(isDensMatr))},
    });
}


void printQuregAutoDeployments(bool isDensMatr) {

    // build all table rows dynamically before print
    std::vector<std::tuple<std::string, std::string>> rows;

    // we will get auto-deployment for every possible number of qubits; silly but cheap and robust!
    int useDistrib,  useGpuAccel,  useMulti;
    int prevDistrib, prevGpuAccel, prevMulti;

    // assume all deployments disabled for 1 qubit
    prevDistrib  = 0;
    prevGpuAccel = 0;
    prevMulti    = 0;

    // test to theoretically max #qubits, surpassing max that can fit in RAM and GPUs, because
    // auto-deploy will still try to deploy there to (then subsequent validation will fail)
    int maxQubits = mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, globalEnvPtr->numNodes);

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
        std::string value = "";
        if (useDistrib)
            value += "[mpi] ";
        if (useGpuAccel)
            value += "[gpu] ";
        if (useMulti)
            value += "[omp] ";

        // log the #qubits of the deployment change
        rows.push_back({form_str(numQubits) + " qubits", value});

        // skip subsequent qubits with the same deployments
        prevDistrib  = useDistrib;
        prevGpuAccel = useGpuAccel;
        prevMulti    = useMulti;
    }

    // tailor table title to type of Qureg
    std::string prefix = (isDensMatr)? "density matrix" : "statevector";
    std::string title = prefix + " autodeployment";
    form_printTable(title, rows);
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

    return (int) (globalEnvPtr != NULL);
}


QuESTEnv getQuESTEnv() {
    validate_envInit(__func__);

    // returns a copy, so cheeky users calling memcpy() upon const struct still won't mutate
    return *globalEnvPtr;
}


void finalizeQuESTEnv() {
    validate_envInit(__func__);

    if (globalEnvPtr->isGpuAccelerated && gpu_isCuQuantumCompiled())
        gpu_finalizeCuQuantum();

    if (globalEnvPtr->isDistributed)
        comm_end();

    // free global env's heap memory and flag as not active
    free(globalEnvPtr);
    globalEnvPtr = NULL;

    // flag that the environment was finalised, to ensure it is never re-initialised
    hasEnvBeenFinalized = true;
}


void reportQuESTEnv() {
    validate_envInit(__func__);

    // TODO: add function to write this output to file (useful for HPC debugging)

    // only root node reports (but no synch necesary)
    if (globalEnvPtr->rank != 0)
        return;

    std::cout << "QuEST execution environment:" << std::endl;

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
}


// end de-mangler
}
