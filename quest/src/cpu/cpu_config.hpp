/** @file
 * Utility signatures for querying CPU hardware.
 */

#ifndef CPU_CONFIG_HPP
#define CPU_CONFIG_HPP

#include "types.h"



/*
 * OPENMP CONFIG
 */

bool cpu_isOpenmpCompiled();

int cpu_getCurrentNumThreads();

int cpu_getNumOpenmpProcessors();



/*
 * MEMORY MANAGEMENT
 */

qcomp* cpu_allocAmps(qindex numLocalAmps);

void cpu_deallocAmps(qcomp* amps);



#endif // CPU_CONFIG_HPP