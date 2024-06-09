/** @file
 * Utility functions for querying GPU hardware.
 */

#include "quest/include/modes.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/comm/communication.hpp"


#if ENABLE_GPU_ACCELERATION && ! (defined(__NVCC__) || defined(__HIPCC__))
    #error \
        "Attempted to compile config.cpp in GPU-accelerated mode with a non-GPU compiler. "\
        "Please compile this file with a CUDA (nvcc) or ROCm (hipcc) compiler."
#endif


#if ENABLE_GPU_ACCELERATION
    #include <cuda.h>
    #include <cuda_runtime.h>
#endif



bool gpu_isGpuCompiled() {
    return (bool) ENABLE_GPU_ACCELERATION;
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
    cudaGetDeviceCount(&num);
    return num;

#else
    error_gpuQueriedButGpuNotCompiled();
    return -1;
#endif
}


void gpu_bindLocalGPUsToNodes(int rank) {
#if ENABLE_GPU_ACCELERATION

    int numLocalGpus = gpu_getNumberOfLocalGpus();
    int localGpuInd = rank % numLocalGpus;
    cudaSetDevice(localGpuInd);

#else
    error_gpuQueriedButGpuNotCompiled();
#endif 
}


size_t gpu_getCurrentAvailableMemoryInBytes() {
#if ENABLE_GPU_ACCELERATION

    // note that in distributed settings, all GPUs
    // are being simultaneously queried, and it is
    // possible their values differ per-node

    size_t free, total;
    cudaMemGetInfo(&free, &total);
    return free;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}


size_t gpu_getTotalMemoryInBytes() {
#if ENABLE_GPU_ACCELERATION

    size_t free, total;
    cudaMemGetInfo(&free, &total);
    return total;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}