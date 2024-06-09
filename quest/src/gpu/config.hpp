/** @file
 * Utility functions for querying GPU hardware, used by gpu.cpp.
 */

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <cstdlib>



bool gpu_isGpuCompiled();

bool gpu_isGpuAvailable();

bool gpu_isDirectGpuCommPossible();

int gpu_getNumberOfLocalGpus();

void gpu_bindLocalGPUsToNodes(int rank);

size_t gpu_getCurrentAvailableMemoryInBytes();

size_t gpu_getTotalMemoryInBytes();



#endif // CONFIG_HPP