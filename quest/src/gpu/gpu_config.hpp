/** @file
 * Utility functions for querying GPU hardware. Note this header is included by
 * /core/ and parsed by non-CUDA compilers, so it must not contain any CUDA-specific
 * signatures.
 */

#ifndef GPU_CONFIG_HPP
#define GPU_CONFIG_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include <cstdlib>



/*
 * CUDA ERROR HANDLING
 */

#if ENABLE_GPU_ACCELERATION

#define CUDA_CHECK(cmd) do { assertCudaCallSucceeded((cmd), #cmd, __func__, __FILE__, __LINE__); } while (0)

void assertCudaCallSucceeded(int code, const char* call, const char* caller, const char* file, int line);

#endif 



/*
 * HARDWARE AVAILABILITY
 */

bool gpu_isGpuCompiled();

bool gpu_isCuQuantumCompiled();

bool gpu_isGpuAvailable();

bool gpu_isDirectGpuCommPossible();

int gpu_getNumberOfLocalGpus();

size_t gpu_getCurrentAvailableMemoryInBytes();

size_t gpu_getTotalMemoryInBytes();

bool gpu_doesGpuSupportMemPools();



/*
 * ENVIRONMENT MANAGEMENT
 */

void gpu_initCuQuantum();

void gpu_finalizeCuQuantum();

void gpu_bindLocalGPUsToNodes(int rank);

void gpu_synch();



/*
 * MEMORY MANAGEMENT
 */

qcomp* gpu_allocAmps(qindex numLocalAmps);

void gpu_deallocAmps(qcomp* amps);

void gpu_copyCpuToGpu(Qureg qureg, qcomp* cpuArr, qcomp* gpuArr, qindex numElems);

void gpu_copyCpuToGpu(Qureg qureg);

void gpu_copyGpuToCpu(Qureg qureg, qcomp* gpuArr, qcomp* cpuArr, qindex numElems);

void gpu_copyGpuToCpu(Qureg qureg);



#endif // GPU_CONFIG_HPP