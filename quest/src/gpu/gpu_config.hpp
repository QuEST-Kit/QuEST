/** @file
 * Utility functions for querying GPU hardware.
 */

#ifndef GPU_CONFIG_HPP
#define GPU_CONFIG_HPP

#include <cstdlib>



bool gpu_isGpuCompiled();

bool gpu_isGpuAvailable();

bool gpu_isDirectGpuCommPossible();

int gpu_getNumberOfLocalGpus();

void gpu_bindLocalGPUsToNodes(int rank);

size_t gpu_getCurrentAvailableMemoryInBytes();

size_t gpu_getTotalMemoryInBytes();



#endif // GPU_CONFIG_HPP