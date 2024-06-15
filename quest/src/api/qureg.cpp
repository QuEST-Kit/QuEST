/** @file
 * API definitions for creating and managing Quregs, and automatically
 * choosing their deployment modes.
 */

#include "quest/include/qureg.h"
#include "quest/include/environment.h"
#include "quest/include/initialisations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <string>
#include <iostream>

// provides substrings (by, na, pm, etc) used by reportQureg
using namespace form_substrings;



/*
 * PRIVATE INNER FUNCTIONS (C++)
 */

Qureg validateAndCreateCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, QuESTEnv env, const char* caller) {

    // ensure deployment is compatible with environment, considering available hardware and their memory capacities
    validate_newQuregParams(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env, caller);

    // automatically overwrite distrib, GPU, and multithread fields which were left as modeflag::USE_AUTO
    autodep_chooseQuregDeployment(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env);

    // throw error if the user had forced multithreading but GPU accel was auto-chosen
    validate_newQuregNotBothMultithreadedAndGpuAccel(useGpuAccel, useMultithread, caller);

    // bind deployment configuration
    Qureg qureg;
    qureg.isGpuAccelerated = useGpuAccel;
    qureg.isDistributed = useDistrib;
    qureg.isMultithreaded = useMultithread;

    // if distributed, inherit config from env
    if (useDistrib) {
        qureg.rank = env.rank;
        qureg.numNodes = env.numNodes;
        qureg.logNumNodes = logBase2(env.numNodes);

    // otherwise set config to single node
    } else {
        qureg.numNodes = 1;
        qureg.logNumNodes = 0;

        // but still retain the env's potentially unique rank
        // because non-distributed quregs are still duplicated
        // between every node, and have duplicate processes.
        qureg.rank = env.rank;
    }

    // set dimension
    qureg.isDensityMatrix = isDensMatr;
    qureg.numQubits = numQubits;
    qureg.numAmps = (isDensMatr)? 
        powerOf2(2*numQubits) : 
        powerOf2(  numQubits);

    // set dimension per node (even if not distributed)
    qureg.logNumAmpsPerNode = (isDensMatr)? 
        (2*numQubits - qureg.logNumNodes) :
        (  numQubits - qureg.logNumNodes);
    qureg.numAmpsPerNode = powerOf2(qureg.logNumAmpsPerNode);

    // pre-set all pointers to NULL so post-alloc validation can safely free or ignore
    qureg.cpuAmps = NULL;
    qureg.gpuAmps = NULL;
    qureg.cpuCommBuffer = NULL;
    qureg.gpuCommBuffer = NULL;

    // always allocate CPU memory
    qureg.cpuAmps = cpu_allocAmps(qureg.numAmpsPerNode);

    // conditionally allocate CPU communication buffer (even if numNodes == 1)
    if (useDistrib)
        qureg.cpuCommBuffer = cpu_allocAmps(qureg.numAmpsPerNode);

    // conditionally allocate GPU memory
    if (useGpuAccel)
        qureg.gpuAmps = gpu_allocAmps(qureg.numAmpsPerNode);

    // conditionally allocate GPU communication buffer
    if (useGpuAccel && useDistrib)
        qureg.gpuCommBuffer = gpu_allocAmps(qureg.numAmpsPerNode);

    // check all allocations succeeded (if any failed, validation frees all non-failures before throwing error)
    bool isNewQureg = true;
    validate_newOrExistingQuregAllocs(qureg, isNewQureg, caller);

    // initialise state to |0> or |0><0|
    initZeroState(qureg); 

    return qureg;
}



/*
 * PRIVATE QUREG REPORTING INNER FUNCTIONS
 */


void printDeploymentInfo(Qureg qureg) {

    form_printTable(
        "deployment", {
        {"isMpiEnabled", qureg.isDistributed},
        {"isGpuEnabled", qureg.isGpuAccelerated},
        {"isOmpEnabled", qureg.isMultithreaded},
    });
}

void printDimensionInfo(Qureg qureg) {

    // 2^N = M or 2^N x 2^N = M
    std::string str;
    str  = bt + form_str(qureg.numQubits);
    str += (qureg.isDensityMatrix)? mu + str : "";
    str += eq + form_str(qureg.numAmps);

    form_printTable(
        "dimension", {
        {"isDensMatr", form_str(qureg.isDensityMatrix)},
        {"numQubits",  form_str(qureg.numQubits)},
        {"numAmps",    str},
    });
}


void printDistributionInfo(Qureg qureg) {

    // not applicable when not distributed
    std::string nodesStr = na;
    std::string ampsStr  = na;

    // 2^N = M per node
    if (qureg.isDistributed) {
        nodesStr = bt + form_str(qureg.logNumNodes)       + eq + form_str(qureg.numNodes);
        ampsStr  = bt + form_str(qureg.logNumAmpsPerNode) + eq + form_str(qureg.numAmpsPerNode);
        ampsStr += pn;
    }

    form_printTable(
        "distribution", {
        {"numNodes", nodesStr},
        {"numAmps",  ampsStr},
    });
}


void printMemoryInfo(Qureg qureg) {

    size_t localArrayMem = mem_getLocalMemoryRequired(qureg.numAmpsPerNode);
    std::string localMemStr = form_str(localArrayMem) + by + ((qureg.isDistributed)? pn : "");

    qindex globalTotalMem = mem_getTotalGlobalMemoryUsed(qureg);
    std::string globalMemStr = (globalTotalMem == 0)? "overflowed" : (form_str(globalTotalMem) + by);

    form_printTable(
        "memory", {
        {"cpuAmps",       (qureg.cpuAmps       == NULL)? na : localMemStr},
        {"gpuAmps",       (qureg.gpuAmps       == NULL)? na : localMemStr},
        {"cpuCommBuffer", (qureg.cpuCommBuffer == NULL)? na : localMemStr},
        {"gpuCommBuffer", (qureg.gpuCommBuffer == NULL)? na : localMemStr},
        {"globalTotal",   globalMemStr},
    });
}



/*
 * PUBLIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {


Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, QuESTEnv env) {
    validate_envInit(env, __func__);

    return validateAndCreateCustomQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env, __func__);
}


Qureg createQureg(int numQubits, QuESTEnv env) {
    validate_envInit(env, __func__);

    int isDensMatr = 0;
    int autoMode = modeflag::USE_AUTO;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, autoMode, autoMode, autoMode, env, __func__);
}


Qureg createDensityQureg(int numQubits, QuESTEnv env) {
    validate_envInit(env, __func__);

    int isDensMatr = 1;
    int autoMode = modeflag::USE_AUTO;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, autoMode, autoMode, autoMode, env, __func__);
}


void destroyQureg(Qureg qureg) {
    validate_quregInit(qureg, __func__);

    // free CPU memory
    cpu_deallocAmps(qureg.cpuAmps);

    // free CPU communication buffer
    if (qureg.isDistributed)
        cpu_deallocAmps(qureg.cpuCommBuffer);

    // free GPU memory
    if (qureg.isGpuAccelerated)
        gpu_deallocAmps(qureg.gpuAmps);

    // free GPU communication buffer
    if (qureg.isGpuAccelerated && qureg.isDistributed)
        gpu_deallocAmps(qureg.gpuCommBuffer);

    // cannot set freed fields to NULL because qureg
    // wasn't passed-by-reference, and isn't returned.
}


void reportQureg(Qureg qureg) {
    validate_quregInit(qureg, __func__);

    // TODO: add function to write this output to file (useful for HPC debugging)

    // precondition: no reportable fields are at risk of overflow as a qindex
    // type, EXCEPT aggregate total memory between distributed nodes (in bytes)

    // only root node reports (but no synch necesary)
    if (qureg.rank != 0)
        return;

    std::cout << "Qureg:" << std::endl;
    printDeploymentInfo(qureg);
    printDimensionInfo(qureg);
    printDistributionInfo(qureg);
    printMemoryInfo(qureg);
}


// end de-mangler
}
