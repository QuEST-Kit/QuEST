/** @file
 * API definitions for creating and managing Quregs, and automatically
 * choosing their deployment modes.
 */

#include "quest/include/qureg.h"
#include "quest/include/environment.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"



/*
 * PRIVATE INNER FUNCTIONS (C++)
 */


void chooseWhetherToUseDistribution(int numQubits, int isDensMatr, int &useDistrib, int useGpuAccel, int numEnvNodes) {

    // if the flag is already set, don't change it
    if (useDistrib != modeflag::USE_AUTO)
        return;

    // if the state has too few amplitudes or columns to uniformly distribute, don't distribute
    qindex numAmpsOrCols = powerOf2(numQubits);
    if (numAmpsOrCols < numEnvNodes) {
        useDistrib = 0;
        return;
    }

    // determine how much memory we'd need to store the whole Qureg in one node
    size_t nonDistribMemReq = mem_getLocalMemoryRequired(numQubits, isDensMatr, 1);

    // attempt to query RAM capacity
    try {
        size_t localCpuMem = mem_tryGetLocalRamCapacityInBytes();

        // if the state cannot fit in RAM, we MUST distribute (regardless of GPU memory; CPU is always allocated)
        if (nonDistribMemReq >= localCpuMem) {
            useDistrib = 1;
            return;
        }
    // if we cannot query RAM, we carry on (we will likely still distribute if RAM would have been exceeded)
    } catch (mem::COULD_NOT_QUERY_RAM &e) {}

    // similarly, if GPU deployment is possible but we cannot currently fit the state in the GPU, we must distribute
    if (useGpuAccel == 1 || useGpuAccel == modeflag::USE_AUTO) {
        size_t localGpuMem = gpu_getCurrentAvailableMemoryInBytes();
        if (nonDistribMemReq >= localGpuMem) {
            useDistrib = 1;
            return;
        }
    }

    // by now, we know that Qureg can definitely fit into a single GPU, or principally fit into RAM,
    // but we still wish to distribute it so that multiple Quregs don't choke up memory. We only avoid
    // distributing if the number of amplitudes per node would be trivially small.

    int effectiveNumQubitsPerNode = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numEnvNodes);
    useDistrib = (effectiveNumQubitsPerNode >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_DISTRIBUTION);
}


void chooseWhetherToUseGpuAcceleration(int numQubits, int isDensMatr, int useDistrib, int &useGpuAccel, int numQuregNodes) {

    // if the flag is already set, don't change it
    if (useGpuAccel != modeflag::USE_AUTO)
        return;

    // determine the 'effective number of qubits' each GPU would have to simulate, if distributed
    int effectiveNumQubits = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);

    // choose to GPU accelerate only if that's not too few
    useGpuAccel = (effectiveNumQubits >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_GPU_ACCELERATION);
}


void chooseWhetherToUseMultithreading(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int &useMultithread, int numQuregNodes) {

    // if the flag is already set (user-given, or inferred from env), don't change it
    if (useMultithread != modeflag::USE_AUTO)
        return;

    // if GPU-aceleration was chosen, disable auto multithreading...
    if (useGpuAccel) {
        useMultithread = 0;
        return;
    }

    // otherwise, we're not GPU-accelerating, and should choose to multithread based on Qureg size
    int effectiveNumQubits = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);
    useMultithread = (effectiveNumQubits >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_MULTITHREADING);
}


void chooseAutomaticDeployments(int numQubits, int isDensMatr, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env, const char* caller) {

    // preconditions:
    //  - the given configuration is compatible with env (assured by validation)
    //  - this means no deployment is forced (=1) which is incompatible with env
    //  - it also means GPU-acceleration and multithreading are not simultaneously forced
    //    (although they may still be left automatic and need explicit revision)

    // disable any automatic deployments not permitted by env (it's gauranteed we never overwrite =1 to =0)
    if (!env.isDistributed)
        useDistrib = 0;

    if (!env.isGpuAccelerated)
        useGpuAccel = 0;

    if (!env.isMultithreaded)
        useMultithread = 0;

    // overwrite any auto options (== modeflag::USE_AUTO)
    chooseWhetherToUseDistribution(numQubits, isDensMatr, useDistrib, useGpuAccel, env.numNodes);

    int numQuregNodes = (useDistrib)? env.numNodes : 1;
    chooseWhetherToUseGpuAcceleration(numQubits, isDensMatr, useDistrib, useGpuAccel, numQuregNodes);
    chooseWhetherToUseMultithreading(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, numQuregNodes);

    // throw error if the user had forced multithreading but the GPU was chosen/forced
    validate_quregNotBothMultithreadedAndGpuAccel(useGpuAccel, useMultithread, caller);
}


Qureg validateAndCreateCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, QuESTEnv env, const char* caller) {

    // ensure deployment is compatible with environment, considering available hardware and their memory capacities.
    // the only thing we cannot gaurantee to be validated is whether the user passed legal useMultithread; they might
    // force multithreading while permitting GPU accel, which if later chosen by auto-deployer, will raise an error.
    validate_quregParams(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env, caller);

    // automatically overwrite distrib, GPU, and thread fields which were left as modeflag::USE_AUTO
    chooseAutomaticDeployments(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env, caller);

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


Qureg createQureg(int numQubits, QuESTEnv env) {
    validate_existingEnv(env, __func__);

    int isDensMatr = 0;
    int autoMode = modeflag::USE_AUTO;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, autoMode, autoMode, autoMode, env, __func__);
}


Qureg createDensityQureg(int numQubits, QuESTEnv env) {
    validate_existingEnv(env, __func__);

    int isDensMatr = 1;
    int autoMode = modeflag::USE_AUTO;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, autoMode, autoMode, autoMode, env, __func__);
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
