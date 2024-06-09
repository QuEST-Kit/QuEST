/** @file
 * Utility signatures for querying CPU hardware.
 */

#ifndef CPU_CONFIG_HPP
#define CPU_CONFIG_HPP



/*
 * OPENMP CONFIG
 */

bool cpu_isOpenmpCompiled();

int cpu_getCurrentNumThreads();

int cpu_getNumOpenmpProcessors();



#endif // CPU_CONFIG_HPP