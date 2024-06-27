/** @file
 * Utility functions for querying GPU hardware.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"


#if ENABLE_GPU_ACCELERATION && ! (defined(__NVCC__) || defined(__HIPCC__))
    #error \
        "Attempted to compile gpu_config.cpp in GPU-accelerated mode with a non-GPU compiler. "\
        "Please compile this file with a CUDA (nvcc) or ROCm (hipcc) compiler."
#endif


#if ENABLE_GPU_ACCELERATION
    #include <cuda.h>
    #include <cuda_runtime.h>
#endif



/*
 * CUDA ERROR HANDLING
 *
 * which is only defined when CUDA-compiling, since it is invoked only a macro (defined
 * in gpu_config.hpp) which wraps CUDA API calls
 */

#if ENABLE_GPU_ACCELERATION

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
 * ENABLE_CUQUANTUM is 1, but are otherwise defaulted to
 * the internal errors below. This slight inelegance
 * enables us to keep gpu_cuquantum.hpp as a single header
 * file, without exposing it to code beyond gpu/
 */


#if ! ENABLE_CUQUANTUM

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
    return (bool) ENABLE_GPU_ACCELERATION;
}


bool gpu_isCuQuantumCompiled() {
    return (bool) ENABLE_CUQUANTUM;
}


bool gpu_isGpuAvailable() {
#if ENABLE_GPU_ACCELERATION

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
#if ENABLE_GPU_ACCELERATION

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
#if ENABLE_GPU_ACCELERATION

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
#if ENABLE_GPU_ACCELERATION

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
#if ENABLE_GPU_ACCELERATION

    size_t free, total;
    CUDA_CHECK( cudaMemGetInfo(&free, &total) );
    return total;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}


bool gpu_doesGpuSupportMemPools() {
#if ENABLE_GPU_ACCELERATION

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
#if ENABLE_GPU_ACCELERATION

    int numLocalGpus = gpu_getNumberOfLocalGpus();
    int localGpuInd = rank % numLocalGpus;
    CUDA_CHECK( cudaSetDevice(localGpuInd) );

#else
    error_gpuQueriedButGpuNotCompiled();
#endif 
}


void gpu_synch() {
#if ENABLE_GPU_ACCELERATION

    CUDA_CHECK( cudaDeviceSynchronize() );

#endif

    // TODO: validation? eh
}




/*
 * MEMORY MANAGEMENT
 */


qcomp* gpu_allocAmps(qindex numLocalAmps) {
#if ENABLE_GPU_ACCELERATION

    size_t numBytes = mem_getLocalMemoryRequired(numLocalAmps);

    qcomp* ptr;
    CUDA_CHECK( cudaMalloc(&ptr, numBytes) );
    return ptr;

#else
    error_gpuAllocButGpuNotCompiled();
    return NULL;
#endif
}


void gpu_deallocAmps(qcomp* amps) {
#if ENABLE_GPU_ACCELERATION

    CUDA_CHECK( cudaFree(amps) );

#else
    error_gpuDeallocButGpuNotCompiled();
#endif
}


void gpu_copyCpuToGpu(Qureg qureg, qcomp* cpuArr, qcomp* gpuArr, qindex numElems) {
#if ENABLE_GPU_ACCELERATION

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
#if ENABLE_GPU_ACCELERATION

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
