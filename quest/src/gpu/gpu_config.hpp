/** @file
 * Utility signatures for querying GPU hardware,
 * and allocating and copying GPU VRAM data.
 * 
 * Note that this header is included by /core/ and
 * parsed by non-CUDA compilers, so must never contain
 * any CUDA-specific signatures
 * 
 * @author Tyson Jones
 */

#ifndef GPU_CONFIG_HPP
#define GPU_CONFIG_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"



/*
 * CUDA ERROR HANDLING
 */

#if COMPILE_CUDA

#define CUDA_CHECK(cmd) \
    assertCudaCallSucceeded((int) (cmd), #cmd, __func__, __FILE__, __LINE__)

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

int gpu_getComputeCapability();

size_t gpu_getCurrentAvailableMemoryInBytes();

size_t gpu_getTotalMemoryInBytes();

bool gpu_doesGpuSupportMemPools();

qindex gpu_getMaxNumConcurrentThreads();



/*
 * ENVIRONMENT MANAGEMENT
 */

void gpu_bindLocalGPUsToNodes();

bool gpu_areAnyNodesBoundToSameGpu();

void gpu_sync();

void gpu_initCuQuantum();

void gpu_finalizeCuQuantum();



/*
 * MEMORY MANAGEMENT
 */

qcomp* gpu_allocArray(qindex numLocalAmps);
void gpu_deallocArray(qcomp* amps);

void gpu_copyArray(qcomp* dest, qcomp* src, qindex dim);

void gpu_copyCpuToGpu(qcomp* cpuArr, qcomp* gpuArr, qindex numElems);
void gpu_copyGpuToCpu(qcomp* gpuArr, qcomp* cpuArr, qindex numElems);

void gpu_copyGpuToCpu(Qureg qureg, qcomp* gpuArr, qcomp* cpuArr, qindex numElems);
void gpu_copyGpuToCpu(Qureg qureg);

void gpu_copyCpuToGpu(Qureg qureg, qcomp* cpuArr, qcomp* gpuArr, qindex numElems);
void gpu_copyCpuToGpu(Qureg qureg);

void gpu_copyCpuToGpu(CompMatr matr);
void gpu_copyGpuToCpu(CompMatr matr);

void gpu_copyCpuToGpu(DiagMatr matr);
void gpu_copyGpuToCpu(DiagMatr matr);

void gpu_copyCpuToGpu(SuperOp op);
void gpu_copyGpuToCpu(SuperOp op);

// funnily, there is no need for GPU-to-CPU or FullStateDiagMatr;
// the invoking printer.cpp function uses localiser_fullstatediagamtr_...
// copying to handle the distributed nuisance, which spoofs a Qureg
void gpu_copyCpuToGpu(FullStateDiagMatr matr);


/*
 * CACHE MANAGEMENT
 */

qcomp* gpu_getCacheOfSize(qindex numElemsPerThread, qindex numThreads);

void gpu_clearCache();

size_t gpu_getCacheMemoryInBytes();



#endif // GPU_CONFIG_HPP