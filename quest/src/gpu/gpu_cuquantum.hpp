/** @file
 * Subroutines which invoke cuQuantum. This file is only ever included
 * when COMPILE_CUQUANTUM=1 and COMPILE_CUDA=1 so it can 
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
#include "quest/src/gpu/gpu_types.hpp"

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
int CUQUANTUM_MEM_POOL_BYTES = 16*(1<<15); // 1 MiB ~ 8 qubit complex<double> matrix


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

    CUDA_CHECK( custatevecApplyMatrix(
        config.cuQuantumHandle, 
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode, 
        matrElems, CUQUANTUM_QCOMP, CUSTATEVEC_MATRIX_LAYOUT_ROW, matrAdj, 
        targs, numTargs,
        ctrls, ctrlVals, numCtrls, 
        CUSTATEVEC_COMPUTE_DEFAULT,
        work, workSize) );
}

void applyMatrix(Qureg qureg, vector<int> ctrls, vector<int> targs, vector<cu_qcomp> matr) {

    applyMatrix(qureg, ctrls.data(), ctrls.size(), targs.data(), targs.size(), matr.data());
}



/*
 * GATES
 */

void cuquantum_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matr) {

    applyMatrix(qureg, {}, {target}, unpackMatrixToCuQcomps(matr));
}



#endif // GPU_CUQUANTUM_HPP