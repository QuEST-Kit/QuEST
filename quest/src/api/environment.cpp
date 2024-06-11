/** @file
 * API definitions for managing QuESTEnv instances, which
 * themselves control and query the deployment environment. 
 */

#include "quest/include/environment.h"
#include "quest/include/precision.h"
#include "quest/include/modes.h"

#include "quest/src/core/memory.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/comm/communication.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <string>
#include <typeinfo>
#include <thread>
#include <vector>
#include <tuple>



/*
 * PRIVATE QUESTENV CREATION INNER FUNCTIONS
 */


QuESTEnv validateAndCreateCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread, const char* caller) {

    // ensure the chosen deployment is compiled and supported by hardware.
    // note that these error messages will be printed by every node because
    // validation occurs before comm_init() below, so all processes spawned
    // by mpirun believe they are each the main rank. This seems unavoidable.
    validate_envNotYetInit(caller);
    validate_envDeploymentMode(useDistrib, useGpuAccel, useMultithread, caller);

    // overwrite deployments left as modeflag::USE_AUTO
    autodep_chooseQuESTEnvDeployment(useDistrib, useGpuAccel, useMultithread);

    // bind deployment info to QuESTEnv (may be overwritten still below)
    QuESTEnv env;
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
    validate_envDistributedBetweenPower2Nodes(env.numNodes, caller);

    // TODO:
    // validate 2^N local GPUs

    // in multi-GPU settings, bind each MPI process to one GPU
    if (useDistrib && useGpuAccel)
        gpu_bindLocalGPUsToNodes(env.rank);

    // TODO: setup RNG

    return env;
}



/*
 * PRIVATE QUESTENV REPORTING INNER FUNCTIONS
 */


// infixes used by QuESTEnv property printing
std::string by = " bytes";
std::string pn = " per node";
std::string pg = " per gpu";
std::string pm = " per machine";
std::string na = "N/A";
std::string un = "unknown";


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
        {"isMpiCompiled", comm_isMpiCompiled()},
        {"isGpuCompiled", gpu_isGpuCompiled()},
        {"isOmpCompiled", cpu_isOpenmpCompiled()},
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
        {"numOmpProcs",   (cpu_isOpenmpCompiled())? form_str(cpu_getNumOpenmpProcessors()) : na},
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
 * PUBLIC FUNCTIONS
 */


// enable invocation by both C and C++ binaries
extern "C" {


QuESTEnv createCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread) {

    return validateAndCreateCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread, __func__);
}


QuESTEnv createQuESTEnv() {

    return validateAndCreateCustomQuESTEnv(modeflag::USE_AUTO, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);
}


void destroyQuESTEnv(QuESTEnv env) {
    validate_existingEnv(env, __func__);

    if (env.isDistributed)
        comm_end();
}


void reportQuESTEnv(QuESTEnv env) {
    validate_existingEnv(env, __func__);

    // we attempt to report properties of available hardware facilities
    // (e.g. number of CPU cores, number of GPUs) even if the environment is not
    // making use of them, to inform the user how they might change deployment.

    // TODO: add function to write this output to file (useful for HPC debugging)

    // only root node reports, other nodes wait for synch
    if (env.rank != 0) {
        comm_synch();
        return;
    }

    std::cout << "QuEST execution environment:" << std::endl;

    bool statevec = false;
    bool densmatr = true;

    printPrecisionInfo();
    printCompilationInfo();
    printDeploymentInfo(env);
    printCpuInfo();
    printGpuInfo();
    printDistributionInfo(env);
    printQuregSizeLimits(statevec, env);
    printQuregSizeLimits(densmatr, env);
    printQuregAutoDeployments(statevec, env);
    printQuregAutoDeployments(densmatr, env);

    // free non-root nodes from synch barrier
    if (env.rank == 0)
        comm_synch();
}


// end de-mangler
}
