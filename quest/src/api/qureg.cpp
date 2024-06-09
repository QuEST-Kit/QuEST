/** @file
 * API definitions for creating and managing Quregs, and automatically
 * choosing their deployment modes.
 */

#include "quest/include/qureg.h"
#include "quest/include/environment.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"



/*
 * PRIVATE INNER FUNCTIONS (C++)
 */


Qureg validateAndCreateCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, QuESTEnv env, const char* caller) {

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
    validate_quregAllocs(qureg, isNewQureg, caller);

    return qureg;
}



/*
 * PUBLIC FUNCTIONS
 */


// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, QuESTEnv env) {
    validate_existingEnv(env, __func__);

    return validateAndCreateCustomQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env, __func__);
}


void destroyQureg(Qureg qureg) {

    // check below arrays are correctly allocated (if not, we do NOT free anything).
    // note this cannot detect whether the Qureg was already destroyed; see final comment
    bool isNewQureg = false;
    validate_quregAllocs(qureg, isNewQureg, __func__);

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


// end de-mangler
#ifdef __cplusplus
}
#endif
