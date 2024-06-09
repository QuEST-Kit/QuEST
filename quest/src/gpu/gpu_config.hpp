/** @file
 * Utility functions for querying GPU hardware.
 */

#ifndef GPU_CONFIG_HPP
#define GPU_CONFIG_HPP

#include "quest/include/types.h"

#include <cstdlib>



/*
 * HARDWARE AVAILABILITY
 */

bool gpu_isGpuCompiled();

bool gpu_isGpuAvailable();

bool gpu_isDirectGpuCommPossible();

int gpu_getNumberOfLocalGpus();

void gpu_bindLocalGPUsToNodes(int rank);

size_t gpu_getCurrentAvailableMemoryInBytes();

size_t gpu_getTotalMemoryInBytes();



/*
 * MEMORY MANAGEMENT
 */

qcomp* gpu_allocAmps(qindex numLocalAmps);

void gpu_deallocAmps(qcomp* amps);



#endif // GPU_CONFIG_HPP