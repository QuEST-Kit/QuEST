/** @file
 * Utility functions for querying GPU hardware.
 */

#ifndef GPU_CONFIG_HPP
#define GPU_CONFIG_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include <cstdlib>



/*
 * HARDWARE AVAILABILITY
 */

bool gpu_isGpuCompiled();

bool gpu_isGpuAvailable();

bool gpu_isDirectGpuCommPossible();

int gpu_getNumberOfLocalGpus();

size_t gpu_getCurrentAvailableMemoryInBytes();

size_t gpu_getTotalMemoryInBytes();



/*
 * ENVIRONMENT MANAGEMENT
 */

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