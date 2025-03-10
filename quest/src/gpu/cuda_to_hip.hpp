#ifndef CUDA_TO_HIP_HPP
#define CUDA_TO_HIP_HPP

#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>

// Set of macros pointing CUDA functions to their analogous HIP function.

#define WARPSIZE 64
static constexpr int maxWarpsPerBlock = 1024/WARPSIZE;

#define CUFFT_D2Z HIPFFT_D2Z
#define CUFFT_FORWARD HIPFFT_FORWARD
#define CUFFT_INVERSE HIPFFT_BACKWARD
#define CUFFT_Z2D HIPFFT_Z2D
#define CUFFT_Z2Z HIPFFT_Z2Z

#define cudaUUID_t hipUUID
#define cudaError hipError_t
#define cudaError_t hipError_t
#define cudaEvent_t hipEvent_t
#define cudaDeviceProp hipDeviceProp_t

#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaErrorInsufficientDriver hipErrorInsufficientDriver
#define cudaErrorNoDevice hipErrorNoDevice
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
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaThreadSynchronize hipDeviceSynchronize
#define cudaErrorMemoryAllocation hipErrorMemoryAllocation
#define cudaDeviceGetAttribute hipDeviceGetAttribute
#define cudaDevAttrMemoryPoolsSupported hipDeviceAttributeMemoryPoolsSupported
#define cudaDevAttrMaxThreadsPerBlock hipDeviceAttributeMaxThreadsPerBlock
#define cudaDevAttrMultiProcessorCount hipDeviceAttributeMultiprocessorCount

#define cufftDestroy hipfftDestroy
#define cufftDoubleComplex hipfftDoubleComplex
#define cufftDoubleReal hipfftDoubleReal
#define cufftExecD2Z hipfftExecD2Z
#define cufftExecZ2D hipfftExecZ2D
#define cufftExecZ2Z hipfftExecZ2Z
#define cufftHandle hipfftHandle
#define cufftPlan3d hipfftPlan3d
#define cufftPlanMany hipfftPlanMany

#define cuFloatComplex hipFloatComplex
#define cuDoubleComplex hipDoubleComplex
#define make_cuFloatComplex make_hipFloatComplex
#define make_cuDoubleComplex make_hipDoubleComplex


static void __attribute__((unused)) check(const hipError_t err, const char *const file, const int line)
{
  if (err == hipSuccess) return;
  fprintf(stderr,"HIP ERROR AT LINE %d OF FILE '%s': %s %s\n",line,file,hipGetErrorName(err),hipGetErrorString(err));
  fflush(stderr);
  exit(err);
}

#endif //CUDA_TO_HIP_HPP