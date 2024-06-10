/** @file
 * Internal functions which query available CPU memory (in an
 * attemptedly OS-agnostic way), determine the maximum
 * number of qubits which can be simulated, and query
 * number types. Note GPU memory querying is performed by 
 * the dedicated GPU backend, though this file is always 
 * compiled (even in GPU mode) because GPU-acceleration still 
 * requires accompanying CPU memory arrays.
 */

#ifndef MEMORY_HPP
#define MEMORY_HPP

#include "quest/include/types.h"



/*
 * HARDWARE QUERYING
 */

namespace mem { typedef bool COULD_NOT_QUERY_RAM; }

qindex mem_tryGetLocalRamCapacityInBytes();



/*
 * MEMORY COST QUERYING
 */

bool mem_canQuregFitInMemory(int numQubits, bool isDensMatr, int numNodes, qindex memBytesPerNode);

int mem_getEffectiveNumStateVecQubitsPerNode(int numQubits, bool isDensMatr, int numNodes);

int mem_getMinNumQubitsForDistribution(int numNodes);

int mem_getMaxNumQubitsWhichCanFitInMemory(bool isDensMatr, int numNodes, qindex memBytesPerNode);

int mem_getMaxNumQubitsBeforeIndexOverflow(bool isDensMatr);

int mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(bool isDensMatr, int numNodes);

size_t mem_getLocalMemoryRequired(int numQubits, bool isDensMatr, int numNodes) ;



#endif // MEMORY_HPP