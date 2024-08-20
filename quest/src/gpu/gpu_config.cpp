/** @file
 * Utility functions for querying GPU hardware.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/environment.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"


#if COMPILE_CUDA && ! (defined(__NVCC__) || defined(__HIPCC__))
    #error \
        "Attempted to compile gpu_config.cpp in GPU-accelerated mode with a non-GPU compiler. "\
        "Please compile this file with a CUDA (nvcc) or ROCm (hipcc) compiler."
#endif


#if COMPILE_CUDA
    #include <cuda.h>
    #include <cuda_runtime.h>
#endif



/*
 * CUDA ERROR HANDLING
 *
 * which is only defined when CUDA-compiling, since it is invoked only a macro (defined
 * in gpu_config.hpp) which wraps CUDA API calls
 */

#if COMPILE_CUDA

void assertCudaCallSucceeded(int result, const char* call, const char* caller, const char* file, int line) {

    // result (int) is actually type cudaError_t but we cannot use this CUDA-defined type
    // in gpu_config.hpp (since it's included by non-CUDA-compiled files), and we wish to keep
    // the signature consistent.
    cudaError_t code = (cudaError_t) result;

    if (result != cudaSuccess)
        error_cudaCallFailed(cudaGetErrorString(code), call, caller, file, line);
}

#endif



/*
 * CUQUANTUM MANAGEMENT
 *
 * these functions are defined in gpu_cuquantum.hpp when
 * COMPILE_CUQUANTUM is 1, but are otherwise defaulted to
 * the internal errors below. This slight inelegance
 * enables us to keep gpu_cuquantum.hpp as a single header
 * file, without exposing it to code beyond gpu/
 */


#if ! COMPILE_CUQUANTUM

void gpu_initCuQuantum() {
    error_cuQuantumInitOrFinalizedButNotCompiled();
}

void gpu_finalizeCuQuantum() {
    error_cuQuantumInitOrFinalizedButNotCompiled();
}

#endif



/*
 * HARDWARE AVAILABILITY
 */


bool gpu_isGpuCompiled() {
    return (bool) COMPILE_CUDA;
}


bool gpu_isCuQuantumCompiled() {
    return (bool) COMPILE_CUQUANTUM;
}


bool gpu_isGpuAvailable() {
#if COMPILE_CUDA

    // DEBUG: cudaGetDeviceProperties is (for some reason) being auto-suffixed with _v2
    // in Cuda 12, which is the only sm=90 compatible version we can use. But then the
    // function is not found in -lcuda and -lcudart, WTF

    // ask CUDA for the number of available "devices"
    int numDevices;
    cudaError_t successCode = cudaGetDeviceCount(&numDevices);

    // if the query failed, we can't use any devices anyway, so we abort
    if (successCode != cudaSuccess)
        return false;

    // so for each reported device...
    for (int deviceInd=0; deviceInd < numDevices; deviceInd++) {

        // query its properties
        struct cudaDeviceProp props;
        successCode = cudaGetDeviceProperties(&props, deviceInd);

        // if the query failed, device is anyway unusable
        if (successCode != cudaSuccess) 
            continue;

        // if the device is a real GPU, it's 'major' compute capability is != 9999 (meaning emulation)
        if (props.major != 9999)
            return true;
    }

    // no non-emulation devices were found
    return false;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif
}


bool gpu_isDirectGpuCommPossible() {
#if COMPILE_CUDA

    if (!comm_isMpiGpuAware())
        return false;

    if (!gpu_isGpuAvailable())
        return false;

    // TODO:
    // and are GPUs compatible?
    // (the above might need to call a GPU-compiled func)

    return true;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif
}


int gpu_getNumberOfLocalGpus() {
#if COMPILE_CUDA

    // TODO: this will over-report, since it may include virtual devices!
    // see gpu_isGpuAvailable()

    int num;
    CUDA_CHECK( cudaGetDeviceCount(&num) );
    return num;

#else
    error_gpuQueriedButGpuNotCompiled();
    return -1;
#endif
}


size_t gpu_getCurrentAvailableMemoryInBytes() {
#if COMPILE_CUDA

    // note that in distributed settings, all GPUs
    // are being simultaneously queried, and it is
    // possible their values differ per-node

    size_t free, total;
    CUDA_CHECK( cudaMemGetInfo(&free, &total) );
    return free;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}


size_t gpu_getTotalMemoryInBytes() {
#if COMPILE_CUDA

    size_t free, total;
    CUDA_CHECK( cudaMemGetInfo(&free, &total) );
    return total;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}


bool gpu_doesGpuSupportMemPools() {
#if COMPILE_CUDA

    // consult only the first device (garuanteed already to exist)
    int deviceId, supports;
    CUDA_CHECK( cudaGetDevice(&deviceId) );
    CUDA_CHECK( cudaDeviceGetAttribute(&supports, cudaDevAttrMemoryPoolsSupported, deviceId) );
    return (bool) supports;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif
}



/*
 * ENVIRONMENT MANAGEMENT
 */


void gpu_bindLocalGPUsToNodes(int rank) {
#if COMPILE_CUDA

    int numLocalGpus = gpu_getNumberOfLocalGpus();
    int localGpuInd = rank % numLocalGpus;
    CUDA_CHECK( cudaSetDevice(localGpuInd) );

#else
    error_gpuQueriedButGpuNotCompiled();
#endif 
}


void gpu_sync() {
#if COMPILE_CUDA

    CUDA_CHECK( cudaDeviceSynchronize() );

#else
    error_gpuSyncedButGpuNotCompiled();
#endif
}



/*
 * MEMORY MANAGEMENT
 */


// initial value of first element of freshly GPU-allocated memory which
// indicates memory has not yet been synced/overwritten since creation.
// This is used by validation to detect users forgetting to sync memory
// of API structures, so should be an arbitrary value users are unlikely 
// to set as the first element (which we will validate anyway). We don't
// do this for Quregs which are always overwritten at creation.
const qcomp UNSYNCED_GPU_MEM_FLAG = qcomp(3.141592653, 12345.67890);


qcomp* gpu_allocAmps(qindex numLocalAmps) {
#if COMPILE_CUDA

    size_t numBytes = mem_getLocalQuregMemoryRequired(numLocalAmps);

    // attempt to malloc
    qcomp* ptr;
    cudaError_t errCode = cudaMalloc(&ptr, numBytes);

    // intercept memory-alloc error and merely return NULL pointer (to be handled by validation)
    if (errCode == cudaErrorMemoryAllocation)
        return NULL;

    // pass all other unexpected errors to internal error handling
    CUDA_CHECK(errCode);

    // mark that the gpu memory is fresh and needs overwriting, by overwriting ptr[0] to unsyc flag
    CUDA_CHECK( cudaMemcpy(ptr, &UNSYNCED_GPU_MEM_FLAG, sizeof(qcomp), cudaMemcpyHostToDevice) );

    return ptr;

#else
    error_gpuAllocButGpuNotCompiled();
    return NULL;
#endif
}


void gpu_deallocAmps(qcomp* amps) {
#if COMPILE_CUDA

    // cudaFree on NULL is fine
    CUDA_CHECK( cudaFree(amps) );

#else
    error_gpuDeallocButGpuNotCompiled();
#endif
}


void gpu_copyCpuToGpu(Qureg qureg, qcomp* cpuArr, qcomp* gpuArr, qindex numElems) {
#if COMPILE_CUDA

    assert_quregIsGpuAccelerated(qureg);

    size_t numBytes = numElems * sizeof(qcomp);
    CUDA_CHECK( cudaMemcpy(gpuArr, cpuArr, numBytes, cudaMemcpyHostToDevice) );

#else
    error_gpuCopyButGpuNotCompiled();
#endif
}

void gpu_copyCpuToGpu(Qureg qureg) {
    gpu_copyCpuToGpu(qureg, qureg.cpuAmps, qureg.gpuAmps, qureg.numAmpsPerNode);
}


void gpu_copyGpuToCpu(Qureg qureg, qcomp* gpuArr, qcomp* cpuArr, qindex numElems) {
#if COMPILE_CUDA

    assert_quregIsGpuAccelerated(qureg);

    size_t numBytes = numElems * sizeof(qcomp);
    CUDA_CHECK( cudaMemcpy(cpuArr, gpuArr, numBytes, cudaMemcpyDeviceToHost) );

#else
    error_gpuCopyButGpuNotCompiled();
#endif
}

void gpu_copyGpuToCpu(Qureg qureg) {
    gpu_copyGpuToCpu(qureg, qureg.gpuAmps, qureg.cpuAmps, qureg.numAmpsPerNode);
}


void gpu_copyCpuToGpu(CompMatr matr) {
#if COMPILE_CUDA

    if (matr.gpuElems == NULL || ! getQuESTEnv().isGpuAccelerated)
        error_gpuCopyButMatrixNotGpuAccelerated();

    // copy each CPU row into flattened GPU memory. we make each memcpy asynch,
    // but it's unclear it helps, nor whether single-stream sync is necessary
    size_t numBytesPerRow = matr.numRows * sizeof(**matr.cpuElems);
    gpu_sync();

    for (qindex r=0; r<matr.numRows; r++) {
        qcomp* cpuRow = matr.cpuElems[r];
        qcomp* gpuSlice = &matr.gpuElems[r*matr.numRows];
        CUDA_CHECK( cudaMemcpyAsync(gpuSlice, cpuRow, numBytesPerRow, cudaMemcpyHostToDevice) );
    }

    gpu_sync();
    
#else
    error_gpuCopyButGpuNotCompiled();
#endif
}


void gpu_copyCpuToGpu(DiagMatr matr) {
#if COMPILE_CUDA

    if (matr.gpuElems == NULL || ! getQuESTEnv().isGpuAccelerated)
        error_gpuCopyButMatrixNotGpuAccelerated();

    size_t numBytes = matr.numElems * sizeof(qcomp);
    CUDA_CHECK( cudaMemcpy(matr.gpuElems, matr.cpuElems, numBytes, cudaMemcpyHostToDevice) );
    
#else
    error_gpuCopyButGpuNotCompiled();
#endif
}


void gpu_copyCpuToGpu(FullStateDiagMatr matr) {
#if COMPILE_CUDA

    if (matr.gpuElems == NULL || ! getQuESTEnv().isGpuAccelerated)
        error_gpuCopyButMatrixNotGpuAccelerated();

    size_t numBytes = matr.numElemsPerNode * sizeof(qcomp);
    CUDA_CHECK( cudaMemcpy(matr.gpuElems, matr.cpuElems, numBytes, cudaMemcpyHostToDevice) );
    
#else
    error_gpuCopyButGpuNotCompiled();
#endif
}


bool gpu_haveGpuAmpsBeenSynced(qcomp* gpuArr) {
#if COMPILE_CUDA

    if (gpuArr == NULL || ! getQuESTEnv().isGpuAccelerated)
        error_gpuCopyButMatrixNotGpuAccelerated();

    // obtain first element from device memory
    qcomp firstElem;
    CUDA_CHECK( cudaMemcpy(&firstElem, gpuArr, sizeof(qcomp), cudaMemcpyDeviceToHost) );

    // check whether it is still the unsync'd flag
    return firstElem != UNSYNCED_GPU_MEM_FLAG;

#else
    error_gpuCopyButGpuNotCompiled();
    return false;
#endif
}


bool gpu_doCpuAmpsHaveUnsyncMemFlag(qcomp firstCpuAmp) {

    // we permit the unsync flag to appear in CPU-only matrices, so we
    // should never be asking this question unless env is GPU-accelerated. 
    // Indeed permitting the flag when CPU-only could astonish users
    // if they later enabled GPU-accel and encounter a validation error
    // (not actually this internal error), but that's less astonishing then 
    // getting a GPU-related error message when running in CPU-mode!
    if (!getQuESTEnv().isGpuAccelerated)
        error_gpuMemSyncQueriedButEnvNotGpuAccelerated();

    return firstCpuAmp == UNSYNCED_GPU_MEM_FLAG;
}



/*
 * CACHE MANAGEMENT
 */


// persistent but variably-sized cache memory used by the any-targ dense
// matrix kernel as working memory, which is lazily runtime expanded when
// necessary, and only ever cleared when triggered by the user
qcomp* gpuCache = NULL;
qindex gpuCacheLen = 0;


qcomp* gpu_getCacheOfSize(qindex numElemsPerThread, qindex numThreads) {
#if COMPILE_CUDA

    qindex numNewElems = numElemsPerThread * numThreads;

    // return existing cache if it's already sufficiently big
    if (numNewElems <= gpuCacheLen)
        return gpuCache;

    // otherwise, resize the cache
    gpuCacheLen = numNewElems;
    CUDA_CHECK( cudaFree(gpuCache) );
    CUDA_CHECK( cudaMalloc(&gpuCache, gpuCacheLen * sizeof *gpuCache) );

    return gpuCache;

#else
    error_gpuCacheModifiedButGpuNotCompiled();
    return NULL;
#endif
}


void gpu_clearCache() {

#if COMPILE_CUDA

    // cudaFree on NULL is fine
    CUDA_CHECK( cudaFree(gpuCache) );

    gpuCache = NULL;
    gpuCacheLen = 0;

#else
    error_gpuCacheModifiedButGpuNotCompiled();
#endif
}


size_t gpu_getCacheMemoryInBytes() {

    // query permitted even when not GPU accelerated
    return gpuCacheLen * sizeof *gpuCache;
}
