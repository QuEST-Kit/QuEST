#ifndef CUDA_TO_HIP_H
#define CUDA_TO_HIP_H

#include <hip/hip_runtime.h>

// Set of macros pointing CUDA functions to their analogous HIP function.

#define WARPSIZE 64
static constexpr int maxWarpsPerBlock = 1024/WARPSIZE;

#define CUFFT_D2Z HIPFFT_D2Z
#define CUFFT_FORWARD HIPFFT_FORWARD
#define CUFFT_INVERSE HIPFFT_BACKWARD
#define CUFFT_Z2D HIPFFT_Z2D
#define CUFFT_Z2Z HIPFFT_Z2Z

#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaError hipError_t
#define cudaError_t hipError_t
#define cudaErrorInsufficientDriver hipErrorInsufficientDriver
#define cudaErrorNoDevice hipErrorNoDevice
#define cudaEvent_t hipEvent_t
#define cudaEventCreate hipEventCreate
#define cudaEventElapsedTime hipEventElapsedTime
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaFree hipFree
#define cudaFreeHost hipHostFree
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError
#define cudaHostAlloc hipHostMalloc
#define cudaHostAllocDefault hipHostMallocDefault
#define cudaMalloc hipMalloc
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemGetInfo hipMemGetInfo
#define cudaMemset hipMemset
#define cudaReadModeElementType hipReadModeElementType
#define cudaSetDevice hipSetDevice
#define cudaSuccess hipSuccess
#define cudaDeviceProp hipDeviceProp_t
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaThreadSynchronize hipDeviceSynchronize


#define cufftDestroy hipfftDestroy
#define cufftDoubleComplex hipfftDoubleComplex
#define cufftDoubleReal hipfftDoubleReal
#define cufftExecD2Z hipfftExecD2Z
#define cufftExecZ2D hipfftExecZ2D
#define cufftExecZ2Z hipfftExecZ2Z
#define cufftHandle hipfftHandle
#define cufftPlan3d hipfftPlan3d
#define cufftPlanMany hipfftPlanMany

static void __attribute__((unused)) check(const hipError_t err, const char *const file, const int line)
{
  if (err == hipSuccess) return;
  fprintf(stderr,"HIP ERROR AT LINE %d OF FILE '%s': %s %s\n",line,file,hipGetErrorName(err),hipGetErrorString(err));
  fflush(stderr);
  exit(err);
}

#endif //CUDA_TO_HIP_H