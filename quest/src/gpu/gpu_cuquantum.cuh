/** @file
 * Subroutines which invoke cuQuantum. This file is only ever included
 * when COMPILE_CUQUANTUM=1 and COMPILE_CUDA=1 so it can 
 * safely invoke CUDA signatures without guards.
 */


// because this file uses a global instance of CuQuantumConfig (not inlined),
// it must never be included by multiple translation units; it can only ever
// be included by gpu_subroutines.cpp
#ifdef GPU_CUQUANTUM_HPP
    #error "File gpu_cuquantum.hpp was erroneously included by multiple source files."
#endif

#ifndef GPU_CUQUANTUM_HPP
#define GPU_CUQUANTUM_HPP


// check preprocessors and compilers are valid before #includes to avoid compile errors

#if ! COMPILE_CUQUANTUM
    #error "A file being compiled somehow included gpu_cuquantum.hpp despite QuEST not being compiled in cuQuantum mode."
#endif

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_cuquantum.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#ifndef __NVCC__
    #error "A file which included gpu_cuquantum.hpp was attemptedly compiled with a non-CUDA compiler."
#endif


#include "quest/src/gpu/gpu_config.hpp"
#include "quest/src/gpu/gpu_types.cuh"

#include <custatevec.h>
#include <vector>

using std::vector;



/*
 * PRECISION FLAG
 *
 * which we will use for both the state precision AND gate matrix precision,
 * because QuEST uses only a single qcomp type for both in the API.
 */

#if (FLOAT_PRECISION == 1)
    #define CUQUANTUM_QCOMP CUDA_C_32F

#elif (FLOAT_PRECISION == 2)
    #define CUQUANTUM_QCOMP CUDA_C_64F

#else
    #error "Build bug; precision.h should have prevented non-float non-double qcomp precision on GPU (and cuQuantum)."
    
#endif



/*
 * ENVIRONMENT MANAGEMENT
 */


// sets the size at which the CUDA memory pool will 
// automatically deallocate temporary memory. Below this
// size, temporary memory structures (like a CompMatr)
// will persist in GPU memory to save time. This is
// only relevant to GPU-mode with cuQuantum enabled,
// and is effected at createQuESTEnv().
size_t CUQUANTUM_MEM_POOL_BYTES = 16*(1<<15); // 1 MiB ~ 8 qubit complex<double> matrix


struct CuQuantumConfig {
    cudaStream_t cuStream;
    custatevecHandle_t cuQuantumHandle;
};

// singleton handle to cuQuantum env needed for applying gates and finalizing env
CuQuantumConfig config;


cudaMemPool_t getExistingMemPool() {

    // validation gaurantees memory pool already exists
    cudaMemPool_t memPool;
    int deviceId;
    CUDA_CHECK( cudaGetDevice(&deviceId) );
    CUDA_CHECK( cudaDeviceGetMemPool(&memPool, deviceId) );
    return memPool;
}


void adjustMemPoolSize(cudaMemPool_t memPool) {

    // find existing memPool threshold (above which memory gets freed at every stream synch)
    size_t currMaxMem;
    CUDA_CHECK( cudaMemPoolGetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &currMaxMem) ); 

    // optionally increase memPool threshold
    if (currMaxMem < CUQUANTUM_MEM_POOL_BYTES)
        CUDA_CHECK( cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &CUQUANTUM_MEM_POOL_BYTES) ); 
}

int allocMemInPool(void* ctx, void** ptr, size_t size, cudaStream_t stream) {
    cudaMemPool_t& pool = *static_cast<cudaMemPool_t*>(ctx);
    return cudaMallocFromPoolAsync(ptr, size, pool, stream); 
}

int deallocMemInPool(void* ctx, void* ptr, size_t size, cudaStream_t stream) {
    return cudaFreeAsync(ptr, stream); 
}


void gpu_initCuQuantum() {

    // create new stream and cuQuantum handle, binding to global config
    CUDA_CHECK( custatevecCreate(&config.cuQuantumHandle) );
    CUDA_CHECK( cudaStreamCreate(&config.cuStream) );

    // get and configure existing mem-pool (for later automatic alloc/dealloc of gate matrices)
    cudaMemPool_t memPool = getExistingMemPool();
    adjustMemPoolSize(memPool);

    // create a temporary cuQuantum memory handler
    custatevecDeviceMemHandler_t memHandler;
    memHandler.ctx = &memPool;
    memHandler.device_alloc = allocMemInPool;
    memHandler.device_free = deallocMemInPool;
    strcpy(memHandler.name, "mempool");

    // bind memory handler and stream to cuQuantum handle
    CUDA_CHECK( custatevecSetDeviceMemHandler(config.cuQuantumHandle, &memHandler) );
    CUDA_CHECK( custatevecSetStream(config.cuQuantumHandle, config.cuStream) );
}


void gpu_finalizeCuQuantum() {

    CUDA_CHECK( cudaStreamDestroy(config.cuStream) );
    CUDA_CHECK( custatevecDestroy(config.cuQuantumHandle) );
}



/*
 * MATRICES
 */


void cuquantum_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    // our SWAP targets are bundled into pairs
    int2 targPairs[] = {{targ1, targ2}};;
    int numTargPairs = 1;

    CUDA_CHECK( custatevecSwapIndexBits(
        config.cuQuantumHandle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        targPairs, numTargPairs,

        // swap mask params seem to be in the reverse order to the remainder of the cuStateVec API
        ctrlStates.data(), ctrls.data(), ctrls.size()) );
}


// there is no bespoke cuquantum_statevec_anyCtrlSwap_subB()


void cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, cu_qcomp* flatMatrElems) {

    // this funciton is called 'subA' instead of just 'sub', because it is also called in 
    // the one-target case whereby it is strictly the embarrassingly parallel _subA scenario

    // do not adjoint matrix
    int adj = 0;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    CUDA_CHECK( custatevecApplyMatrix(
        config.cuQuantumHandle, 
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode, 
        flatMatrElems, CUQUANTUM_QCOMP, CUSTATEVEC_MATRIX_LAYOUT_ROW, adj, 
        targs.data(), targs.size(),
        ctrls.data(), ctrlStates.data(), ctrls.size(), 
        CUSTATEVEC_COMPUTE_DEFAULT,
        work, workSize) );
}


// there is no bespoke cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subB()


void cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, cu_qcomp* flatMatrElems) {

    // apply no permutation matrix
    custatevecIndex_t *perm = nullptr;

    // do not adjoint elems
    int adj = 0;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    CUDA_CHECK( custatevecApplyGeneralizedPermutationMatrix(
        config.cuQuantumHandle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        perm, flatMatrElems, CUQUANTUM_QCOMP, adj, 
        targs.data(), targs.size(), 
        ctrls.data(), ctrlStates.data(), ctrls.size(),
        work, workSize) );
}



/*
 * DECOHERENCE
 */


void cuquantum_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob) {

    // effect the superoperator as a two-qubit diagonal upon a statevector suffix state
    cu_qcomp a = {1,        0};
    cu_qcomp b = {1-2*prob, 0};
    cu_qcomp elems[] = {a, b, b, a};
    vector<int> targs {qubit, qubit + qureg.numQubits};

    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, {}, {}, targs, elems);
}


void cuquantum_densmatr_oneQubitDephasing_subB(Qureg qureg, int ketQubit, qreal prob) {

    // we need to merely scale every amp where ketQubit differs from braBit, which is
    // equivalent to a state-controlled global phase upon a statevector, which is 
    // itself a same-element one-qubit diagonal applied to any target
    int braBit = getBit(qureg.rank, ketQubit - qureg.logNumColsPerNode);
    cu_qcomp fac = {1 - 2*prob, 0};
    cu_qcomp elems[] = {fac, fac};

    // we choose to target the largest possible qubit, expecting best cuStateVec performance
    int targ = qureg.numQubits * 2 - 1; // leftmost bra qubit

    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, {ketQubit}, {!braBit}, {targ}, elems);
}


void cuquantum_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubitA, int qubitB, qreal prob) {

    // TODO:
    // only 75% of the amps are changed, each of which is multiplied by the same scalar,
    // but our below method multiplies all amps with 16 separate scalars - can we accel?
    // we are applying a prefactor b to all states except where the ket & bra qubits
    // are the same, i.e. we skip |00><00|, |01><01|, |10><10|, |11><11|

    // effect the superoperator as a four-qubit diagonal upon a statevector suffix state
    cu_qcomp a = {1,          0};
    cu_qcomp b = {1-4*prob/3, 0};
    cu_qcomp elems[] = {a,b,b,b, b,a,b,b, b,b,a,b, b,b,b,a};
    vector<int> targs {qubitA, qubitB, qubitA + qureg.numQubits, qubitB + qureg.numQubits};

    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, {}, {}, targs, elems);
}


// there is no bespoke cuquantum_densmatr_twoQubitDephasing_subB()



#endif // GPU_CUQUANTUM_HPP