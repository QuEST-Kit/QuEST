/** @file
 * Functions which automatically choose QuESTEnv, Qureg and
 * FullStateDiagMatr deployment parameters, replacing flag 
 * modeflag::USE_AUTO with 0 or 1, depending on the compiled 
 * facilities, available hardware, and object dimensions.
 * 
 * @author Tyson Jones
 * @author Richard Meister (aided in design)
 */

#include "quest/include/modes.h"
#include "quest/include/environment.h"

#include "quest/src/core/memory.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"



/*
 * QUESTENV DEPLOYMENT
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
 * QUREG DEPLOYMENT
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

    // force distribution if GPU deployment is available but we exceed local VRAM; 
    // this is preferable over falling back to CPU-only which would be astonishingly slow 
    if (useGpuAccel == 1 || useGpuAccel == modeflag::USE_AUTO) {
        size_t localGpuMem = gpu_getCurrentAvailableMemoryInBytes();
        if (!mem_canQuregFitInMemory(numQubits, isDensMatr, 1, localGpuMem)) {
            useDistrib = 1;
            return;
        }
    }

    // to reach here, we know that Qureg can fit into the remaining memory of a single GPU, or principally 
    // fit into RAM, but we may still wish to distribute for improved parallelisation and to avoid memory saturation
    int effectiveNumQubitsPerNode = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numEnvNodes);
    useDistrib = (effectiveNumQubitsPerNode >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_DISTRIBUTION);
}


void chooseWhetherToGpuAccelQureg(int numQubits, int isDensMatr, int &useGpuAccel, int numQuregNodes) {

    // if the flag is already set, don't change it
    if (useGpuAccel != modeflag::USE_AUTO)
        return;

    // determine the 'effective number of qubits' each GPU would have to simulate, if distributed...
    int effectiveNumQubits = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);

    // and choose to GPU accelerate only if that's not too few
    useGpuAccel = (effectiveNumQubits >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_GPU_ACCELERATION);

    // notice there was no automatic disabling of GPU acceleration in the scenario that the local
    // partition exceeded GPU memory. This is because such a scenario would be catastrophically
    // slow and astonish users by leaving GPUs idle in intensive simulation. Instead, we auto-deploy
    // to GPU anyway and subsequent validation will notice we exceeded GPU memory and report an error.
}


void chooseWhetherToMultithreadQureg(int numQubits, int isDensMatr, int &useMultithread, int numQuregNodes) {

    // if the flag is already set (user-given, or inferred from env), don't change it
    if (useMultithread != modeflag::USE_AUTO)
        return;

    // otherwise, choose to multithread based on Qureg size
    int effectiveNumQubits = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);
    useMultithread = (effectiveNumQubits >= MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_MULTITHREADING);

    // note the qureg may be simultaneously GPU-accelerated and so never use its
    // multithreaded CPU routines, except in functions which accept multiple Quregs
}


void autodep_chooseQuregDeployment(int numQubits, int isDensMatr, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env) {

    // preconditions:
    //  - the given configuration is compatible with env (assured by prior validation)
    //  - this means no deployment is forced (=1) which is incompatible with env

    // disable any automatic deployments not permitted by env (it's gauranteed we never overwrite =1 to =0)
    if (!env.isDistributed)
        useDistrib = 0;
    if (!env.isGpuAccelerated)
        useGpuAccel = 0;
    if (!env.isMultithreaded)
        useMultithread = 0;

    // disable distribution if env is distributed over 1 node. this can occur
    // because env auto-deployer does not know the number of nodes in advance
    if (env.numNodes == 1)
        useDistrib = 0;

    // overwrite useDistrib
    chooseWhetherToDistributeQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, env.numNodes);
    int numQuregNodes = (useDistrib)? env.numNodes : 1;
    
    // overwrite useGpuAccel
    chooseWhetherToGpuAccelQureg(numQubits, isDensMatr, useGpuAccel, numQuregNodes);

    // overwrite useMultithread
    chooseWhetherToMultithreadQureg(numQubits, isDensMatr, useMultithread, numQuregNodes);
}



/*
 * FULL-STATE DIAGONAL MATRIX DEPLOYMENT
 */


void autodep_chooseFullStateDiagMatrDeployment(int numQubits, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env) {

    // we choose to deploy just like a equivalently-sized statevector Qureg would 
    // deploy, so that automatically-deployed matrices and Quregs are always 
    // compatible. This retains compatibility with density-matrix Quregs too;
    // a non-distributed DiagMatr is compatible with both distributed and not
    // Quregs (because its elements are available to every node), and the GPU
    // memory of the DiagMatr is quadratically smaller than of the density-
    // matrix Qureg, so we can always allocate/copy as needed. Ergo, we pretend
    // the FullStateDiagMatr is a statevector Qureg.
    int isDensMatr = 0;

    autodep_chooseQuregDeployment(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env);
}
