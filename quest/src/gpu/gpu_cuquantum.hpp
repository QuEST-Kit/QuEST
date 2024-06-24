/** @file
 * Subroutines which invoke Thrust. This file is only ever included
 * when ENABLE_CUQUANTUM=1 and ENABLE_GPU_ACCELERATION=1 so it can 
 * safely invoke CUDA signatures without guards.
 */

// because this file uses a global instance of CuQuantumConfig (not inlined),
// it must never be included by multiple translation units; it can only ever
// be included by gpu_subroutines.cpp
#ifdef GPU_CUQUANTUM
    #error "File gpu_cuquantum.hpp was erroneously included by multiple source files."
#endif

#ifndef GPU_CUQUANTUM_HPP
#define GPU_CUQUANTUM_HPP


#if ! ENABLE_CUQUANTUM
    #error "A file being compiled somehow included gpu_cuquantum.hpp despite QuEST not being compiled in cuQuantum mode."
#endif

#if ! ENABLE_GPU_ACCELERATION
    #error "A file being compiled somehow included gpu_cuquantum.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif


#include "quest/include/debug.h"

#include "quest/src/gpu/gpu_types.hpp"

#include <custatevec.h>



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
    cudaGetDevice(&deviceId);
    cudaDeviceGetMemPool(&memPool, deviceId);
    return memPool;
}


void adjustMemPoolSize(cudaMemPool_t memPool) {

    // find existing memPool threshold (above which memory gets freed at every stream synch)
    size_t currMaxMem;
    cudaMemPoolGetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &currMaxMem); 

    // optionally increase memPool threshold
    if (currMaxMem < CUQUANTUM_MEM_POOL_BYTES)
        cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &CUQUANTUM_MEM_POOL_BYTES); 
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
    custatevecCreate(&config.cuQuantumHandle);
    cudaStreamCreate(&config.cuStream);

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
    custatevecSetDeviceMemHandler(config.cuQuantumHandle, &memHandler);
    custatevecSetStream(config.cuQuantumHandle, config.cuStream);
}


void gpu_finalizeCuQuantum() {

    cudaStreamDestroy(config.cuStream);
    custatevecDestroy(config.cuQuantumHandle);
}



/*
 * INTERNAL CUQUANTUM WRAPPERS (to reduce boilerplate)
 */


void applyMatrix(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, cu_qcomp* matrElems) {

    // do not adjoint matrix
    int matrAdj = 0;

    // condition all ctrls on =1 state
    int* ctrlVals = nullptr;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    custatevecApplyMatrix(
        config.cuQuantumHandle, 
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode, 
        matrElems, CUQUANTUM_QCOMP, CUSTATEVEC_MATRIX_LAYOUT_ROW, matrAdj, 
        targs, numTargs,
        ctrls, ctrlVals, numCtrls, 
        CUSTATEVEC_COMPUTE_DEFAULT,
        work, workSize);
}

void applyMatrix(Qureg qureg, std::vector<int> ctrls, std::vector<int> targs, std::vector<cu_qcomp> matr) {

    applyMatrix(qureg, ctrls.data(), ctrls.size(), targs.data(), targs.size(), matr.data());
}



/*
 * GATES
 */

void cuquantum_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matr) {

    applyMatrix(qureg, {}, {target}, unpackMatrixToCuQcomps(matr));
}



#endif // GPU_CUQUANTUM_HPP