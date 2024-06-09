/** @file
 * Internal functions which query available CPU memory (in an
 * attemptedly OS-agnostic way) and determine the maximum
 * number of qubits which can be simualted. Note GPU memory
 * querying is performed by the dedicated GPU backend, 
 * though this file is always compiled (even in GPU mode) 
 * because GPU-acceleration still requires accompanying
 * CPU memory arrays.
 */

#ifndef MEMORY_HPP
#define MEMORY_HPP

#include "quest/include/types.h"



// exception thrown by mem_tryGetTotalRamCapacityInBytes()
namespace mem { typedef bool COULD_NOT_QUERY_RAM; }



qindex mem_tryGetLocalRamCapacityInBytes();

int mem_getEffectiveNumStateVecQubitsPerNode(int numQubits, bool isDensMatr, int numNodes);

bool mem_canQuregFitInMemory(int numQubits, bool isDensMatr, int numNodes, qindex memBytesPerNode);



#endif // MEMORY_HPP