/** @file
 * API definitions for managing QuESTEnv instances, which
 * themselves control and query the deployment environment. 
 */

#include "environment.h"
#include "precision.h"
#include "modes.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <typeinfo>
#include <cstring>
#include <string>
#include <thread>
#include <vector>
#include <tuple>

// provides substrings (by, na, pm, etc) used by reportQuESTEnv
using namespace printer_substrings;

using std::string;



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

    print_table(
        "precision", {
        {"qreal",  printer_getQrealType()  + " (" + printer_toStr(sizeof(qreal))  + by + ")"},
        {"qcomp",  printer_getQcompType()  + " (" + printer_toStr(sizeof(qcomp))  + by + ")"},
        {"qindex", printer_getQindexType() + " (" + printer_toStr(sizeof(qindex)) + by + ")"},
        {"validationEpsilon", printer_toStr(VALIDATION_EPSILON)},
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

    // assume RAM is unknown unless it can be queried
    string ram = un;
    try { 
        ram = printer_toStr(mem_tryGetLocalRamCapacityInBytes()) + by + pm; 
    } catch(mem::COULD_NOT_QUERY_RAM e){};

    // TODO
    // - CPU info e.g. speeds/caches?

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

    // TODO below:
    // - GPU compute capability
    // - GPU #SVMs etc

    print_table(
        "gpu", {
        {"numGpus",       (gpu_isGpuCompiled())? printer_toStr(gpu_getNumberOfLocalGpus()) : un},
        {"gpuDirect",     (gpu_isGpuCompiled())? printer_toStr(gpu_isDirectGpuCommPossible()) : na},
        {"gpuMemPools",   (gpu_isGpuCompiled())? printer_toStr(gpu_doesGpuSupportMemPools()) : na},
        {"gpuMemory",     (gpu_isGpuCompiled())? printer_toStr(gpu_getTotalMemoryInBytes()) + by + pg : na},
        {"gpuMemoryFree", (gpu_isGpuCompiled())? printer_toStr(gpu_getTotalMemoryInBytes()) + by + pg : na},
    });
}


void printDistributionInfo() {

    print_table(
        "distribution", {
        {"isMpiGpuAware", (comm_isMpiCompiled())? printer_toStr(comm_isMpiGpuAware()) : na},
        {"numMpiNodes",   printer_toStr(globalEnvPtr->numNodes)},
    });
}


void printQuregSizeLimits(bool isDensMatr) {

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
        {"maxQubitsForMemOverflow", printer_toStr(mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, numNodes))},
        {"maxQubitsForIndOverflow", printer_toStr(mem_getMaxNumQubitsBeforeIndexOverflow(isDensMatr))},
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

    return (int) (globalEnvPtr != NULL);
}


QuESTEnv getQuESTEnv() {
    validate_envIsInit(__func__);

    // returns a copy, so cheeky users calling memcpy() upon const struct still won't mutate
    return *globalEnvPtr;
}


void finalizeQuESTEnv() {
    validate_envIsInit(__func__);

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
    validate_envIsInit(__func__);

    // TODO: add function to write this output to file (useful for HPC debugging)

    print( "QuEST execution environment:");

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
