/** @file
 * Functions which automatically choose QuESTEnv and Qureg deployment
 * parameters, replacing flag modeflag::USE_AUTO with 0 or 1,
 * depending on the compiled facilities, available hardware, and
 * Qureg dimensions.
 */

#include "modes.h"
#include "environment.h"

#include "../core/memory.hpp"
#include "../core/autodeployer.hpp"
#include "../comm/comm_config.hpp"
#include "../cpu/cpu_config.hpp"
#include "../gpu/gpu_config.hpp"



/*
 * QUESTENV DEPLOYMENTS
 */


void autodep_chooseQuESTEnvDeployment(int &useDistrib, int &useGpuAccel, int &useMultithread) {

    // replace automatic flags with all available deployments...
    if (useDistrib == modeflag::USE_AUTO)
        useDistrib = comm_isMpiCompiled();
    
    // where we require GPU is compiled AND available
    if (useGpuAccel == modeflag::USE_AUTO)
        useGpuAccel = gpu_isGpuCompiled() && gpu_isGpuAvailable();

    // and we require more than 1 thread available at QuESTEnv creation
    if (useMultithread == modeflag::USE_AUTO)
        useMultithread = (cpu_isOpenmpCompiled())? (cpu_getCurrentNumThreads() > 1) : 0;
}



/*
 * QUREG DEPLOYMENTS 
 */


void chooseWhetherToDistributeQureg(int numQubits, int isDensMatr, int &useDistrib, int useGpuAccel, int numEnvNodes) {

    // precondition: if the environment doesn't support distribution, useDistrib will already be 0

    // if the flag is already set, don't change it
    if (useDistrib != modeflag::USE_AUTO)
        return;

    // if the state has too few amplitudes or columns to uniformly distribute, don't distribute
    if (numQubits < mem_getMinNumQubitsForDistribution(numEnvNodes)) {
        useDistrib = 0;
        return;
    }

    // force distribution if necessary local memory exceeds single node RAM (if it's queryable)
    try {
        size_t localCpuMem = mem_tryGetLocalRamCapacityInBytes();
        if (!mem_canQuregFitInMemory(numQubits, isDensMatr, 1, localCpuMem)) {
            useDistrib = 1;
            return;
        }
    // it's ok if we cannot query RAM; if we'd have exceeded it, it's likely we'll exceed auto-threshold and will still distribute
    } catch (mem::COULD_NOT_QUERY_RAM &e) {}

    // force distribution if GPU deployment is possible but we exceed local VRAM
    if (useGpuAccel == 1 || useGpuAccel == modeflag::USE_AUTO) {
        size_t localGpuMem = gpu_getCurrentAvailableMemoryInBytes();
        if (!mem_canQuregFitInMemory(numQubits, isDensMatr, 1, localGpuMem)) {
            useDistrib = 1;
            return;
        }
    }

    // by now, we know that Qureg can definitely fit into a single GPU, or principally fit into RAM,
    // but we may still wish to distribute it so that multiple Quregs don't choke up memory.
    int effectiveNumQubitsPerNode = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numEnvNodes);
    useDistrib = (effectiveNumQubitsPerNode >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_DISTRIBUTION);
}


void chooseWhetherToGpuAccelQureg(int numQubits, int isDensMatr, int useDistrib, int &useGpuAccel, int numQuregNodes) {

    // if the flag is already set, don't change it
    if (useGpuAccel != modeflag::USE_AUTO)
        return;

    // determine the 'effective number of qubits' each GPU would have to simulate, if distributed
    int effectiveNumQubits = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);

    // choose to GPU accelerate only if that's not too few
    useGpuAccel = (effectiveNumQubits >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_GPU_ACCELERATION);

    // notice there was no automatic disabling of GPU acceleration in the scenario that the local
    // partition exceeded GPU memory. This is because such a scenario would be catastrophically
    // slow and astonish users by leaving GPUs idle in intensive simulation. Instead, we auto-deploy
    // to GPU and subsequent validation will notice we exceeded GPU memory.
}


void chooseWhetherToMultithreadQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int &useMultithread, int numQuregNodes) {

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


void autodep_chooseQuregDeployment(int numQubits, int isDensMatr, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env) {

    // preconditions:
    //  - the given configuration is compatible with env (assured by prior validation)
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
    chooseWhetherToDistributeQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, env.numNodes);
    int numQuregNodes = (useDistrib)? env.numNodes : 1;
    chooseWhetherToGpuAccelQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, numQuregNodes);
    chooseWhetherToMultithreadQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, numQuregNodes);
}
