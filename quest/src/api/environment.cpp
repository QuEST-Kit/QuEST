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
#include "../gpu/gpu_config.hpp"

#include <iostream>
#include <string>
#include <typeinfo>
#include <thread>
#include <vector>
#include <tuple>

// provides substrings (by, na, pm, etc) used by reportQuESTEnv
using namespace form_substrings;



/*
 * PRIVATE QUESTENV SINGLETON
 *
 * Global to this file, accessible to other files only through getQuESTEnv().
 * The private bools indicate whether the QuEST environment is currently active,
 * and whether it has been finalised after being active respectively. This
 * difference is important, since the QuEST environment can only ever be initialised
 * once, even after finalisation.
 */

static QuESTEnv questEnv;
static bool isQuESTInit  = false;
static bool isQuESTFinal = false;



/*
 * PRIVATE QUESTENV INITIALISATION INNER FUNCTIONS
 */


void validateAndInitCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread, const char* caller) {

    // ensure that we are never re-initialising QuEST (even after finalize) because
    // this leads to undefined behaviour in distributed mode, as per the MPI
    validate_envNeverInit(isQuESTInit, isQuESTFinal, caller);

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

    // create a private, local env (so we don't modify global env in case of error)
    QuESTEnv env;

    // bind deployment info to QuESTEnv (may be overwritten still below)
    env.isDistributed = useDistrib;
    env.isGpuAccelerated = useGpuAccel;
    env.isMultithreaded = useMultithread;

    // assume no distribution, then revise below
    env.rank = 0;
    env.numNodes = 1;

    // initialise distribution (even for a single node), though we may subsequently error and finalize
    if (useDistrib) {
        comm_init();
        env.rank = comm_getRank();
        env.numNodes = comm_getNumNodes();
    }

    // validates numNodes=2^N and otherwise calls comm_end() before throwing error
    validate_newEnvDistributedBetweenPower2Nodes(env.numNodes, caller);

    // TODO:
    // validate 2^N local GPUs

    // in multi-GPU settings, bind each MPI process to one GPU
    if (useDistrib && useGpuAccel)
        gpu_bindLocalGPUsToNodes(env.rank);

    // in GPU settings, if cuQuantum is being used, initialise it
    if (useGpuAccel && gpu_isCuQuantumCompiled())
        gpu_initCuQuantum();

    // TODO: setup RNG

    // overwrite attributes of the global, static env
    questEnv = env;

    // declare successful initialisation
    isQuESTInit = true;
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


void printDeploymentInfo(QuESTEnv env) {

    form_printTable(
        "deployment", {
        {"isMpiEnabled", env.isDistributed},
        {"isGpuEnabled", env.isGpuAccelerated},
        {"isOmpEnabled", env.isMultithreaded},
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


void printDistributionInfo(QuESTEnv env) {

    form_printTable(
        "distribution", {
        {"isMpiGpuAware", (comm_isMpiCompiled())? form_str(comm_isMpiGpuAware()) : na},
        {"numMpiNodes",   form_str(env.numNodes)},
    });
}


void printQuregSizeLimits(bool isDensMatr, QuESTEnv env) {

    // by default, CPU limits are unknown (because memory query might fail)
    std::string maxQbForCpu = un;
    std::string maxQbForMpiCpu = un;

    // max CPU registers are only determinable if RAM query succeeds
    try {
        qindex cpuMem = mem_tryGetLocalRamCapacityInBytes();
        maxQbForCpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, 1, cpuMem));

        // and the max MPI sizes are only relevant when env is distributed
        if (env.isDistributed)
            maxQbForMpiCpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, env.numNodes, cpuMem));

        // when MPI irrelevant, change their status from "unknown" to "N/A"
        else
            maxQbForMpiCpu = na;

    // no problem if we can't query RAM; we simply don't report relevant limits
    } catch(mem::COULD_NOT_QUERY_RAM e) {};

    // GPU limits are default N/A because they're always determinable when relevant
    std::string maxQbForGpu = na;
    std::string maxQbForMpiGpu = na;

    // max GPU registers only relevant if env is GPU-accelerated
    if (env.isGpuAccelerated) {
        qindex gpuMem = gpu_getCurrentAvailableMemoryInBytes();
        maxQbForGpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, 1, gpuMem));

        // and the max MPI sizes are further only relevant when env is distributed 
        if (env.isDistributed)
            maxQbForMpiGpu = form_str(mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, env.numNodes, gpuMem));
    }

    // tailor table title to type of Qureg
    std::string prefix = (isDensMatr)? "density matrix" : "statevector";
    std::string title = prefix + " limits";

    form_printTable(
        title, {
        {"minQubitsForMpi",     (env.numNodes>1)? form_str(mem_getMinNumQubitsForDistribution(env.numNodes)) : na},
        {"maxQubitsForCpu",     maxQbForCpu},
        {"maxQubitsForGpu",     maxQbForGpu},
        {"maxQubitsForMpiCpu",  maxQbForMpiCpu},
        {"maxQubitsForMpiGpu",  maxQbForMpiGpu},
        {"maxQubitsForMemOverflow", form_str(mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, env.numNodes))},
        {"maxQubitsForIndOverflow", form_str(mem_getMaxNumQubitsBeforeIndexOverflow(isDensMatr))},
    });
}


void printQuregAutoDeployments(bool isDensMatr, QuESTEnv env) {

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
    int maxQubits = mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, env.numNodes);

    for (int numQubits=1; numQubits<maxQubits; numQubits++) {

        // re-choose auto deployment
        useDistrib  = modeflag::USE_AUTO;
        useGpuAccel = modeflag::USE_AUTO;
        useMulti    = modeflag::USE_AUTO;;
        autodep_chooseQuregDeployment(numQubits, isDensMatr, useDistrib, useGpuAccel, useMulti, env);

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

    return (int) isQuESTInit;
}


QuESTEnv getQuESTEnv() {
    validate_envInit(__func__);

    return questEnv;
}


void finalizeQuESTEnv() {
    validate_envInit(__func__);

    if (questEnv.isGpuAccelerated && gpu_isCuQuantumCompiled())
        gpu_finalizeCuQuantum();

    if (questEnv.isDistributed)
        comm_end();

    isQuESTInit = false;
    isQuESTFinal = true;
}


void reportQuESTEnv() {
    validate_envInit(__func__);

    // TODO: add function to write this output to file (useful for HPC debugging)

    // only root node reports (but no synch necesary)
    if (questEnv.rank != 0)
        return;

    std::cout << "QuEST execution environment:" << std::endl;

    bool statevec = false;
    bool densmatr = true;

    // we attempt to report properties of available hardware facilities
    // (e.g. number of CPU cores, number of GPUs) even if the environment is not
    // making use of them, to inform the user how they might change deployment.
    printPrecisionInfo();
    printCompilationInfo();
    printDeploymentInfo(questEnv);
    printCpuInfo();
    printGpuInfo();
    printDistributionInfo(questEnv);
    printQuregSizeLimits(statevec, questEnv);
    printQuregSizeLimits(densmatr, questEnv);
    printQuregAutoDeployments(statevec, questEnv);
    printQuregAutoDeployments(densmatr, questEnv);
}


// end de-mangler
}
