/** @file
 * CPU signatures of the subroutines called by accelerator.cpp. 
 */

#ifndef OMP_SUBROUTINES_HPP
#define OMP_SUBROUTINES_HPP



/*
 * OPENMP CONFIG
 */

bool cpu_isOpenmpCompiled();

int cpu_getCurrentNumThreads();

int cpu_getNumOpenmpProcessors();



#endif // OMP_SUBROUTINES_HPP