/** @file
 * Internal signatures which query available CPU memory (in an
 * attemptedly OS-agnostic way), and provide needed memory
 * querents. Note GPU memory querying is performed by 
 * the dedicated GPU backend, though this file is always 
 * compiled (even in GPU mode) because GPU-acceleration still 
 * requires accompanying CPU memory arrays. This file does not
 * perform any allocation of memory; that is instead performed
 * by cpu_config.cpp, to be symmetric with the GPU-memory
 * allocators in gpu_config.cpp, and use NUMA strategies.
 * 
 * @author Tyson Jones
 */

#ifndef MEMORY_HPP
#define MEMORY_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"



/*
 * HARDWARE QUERYING
 */


namespace mem { typedef bool COULD_NOT_QUERY_RAM; }

qindex mem_tryGetLocalRamCapacityInBytes();



/*
 * MEMORY USAGE
 */


int mem_getEffectiveNumStateVecQubitsPerNode(int numQubits, bool isDensMatr, int numNodes);

qindex mem_getTotalGlobalMemoryUsed(Qureg qureg);



/*
 * MEMORY REQUIRED
 */

size_t mem_getLocalQuregMemoryRequired(int numQubits, bool isDensityMatr, int numNodes);
size_t mem_getLocalQuregMemoryRequired(qindex numAmpsPerNode);

size_t mem_getLocalMatrixMemoryRequired(int numQubits, bool isDenseMatrix, int numNodes);

size_t mem_getLocalSuperOpMemoryRequired(int numQubits);



/*
 * QUBIT BOUNDS
 */


int mem_getMaxNumQuregQubitsWhichCanFitInMemory(bool isDensityMatrix, int numNodes, qindex memBytesPerNode);

int mem_getMaxNumMatrixQubitsWhichCanFitInMemory(bool isDenseMatrix, int numNodes, qindex memBytesPerNode);

int mem_getMaxNumSuperOpQubitsWhichCanFitInMemory(qindex memBytesPerNode);


int mem_getMinNumQubitsForDistribution(int numNodes);


int mem_getMaxNumQuregQubitsBeforeIndexOverflow(bool isDensityMatrix);

int mem_getMaxNumMatrixQubitsBeforeIndexOverflow(bool isDenseMatrix);

int mem_getMaxNumSuperOpQubitsBeforeIndexOverflow();

qindex mem_getMaxNumKrausMapMatricesBeforeIndexOverflow(int numQubits);


int mem_getMaxNumQuregQubitsBeforeGlobalMemSizeofOverflow(bool isDensityMatrix, int numNodes);

int mem_getMaxNumMatrixQubitsBeforeGlobalMemSizeofOverflow(bool isDenseMatrix, int numNodes);

int mem_getMaxNumSuperOpQubitsBeforeGlobalMemSizeofOverflow();

qindex mem_getMaxNumKrausMapMatricesBeforeLocalMemSizeofOverflow(int numQubits);



/*
 * SUFFICIENT MEMORY QUERYING
 */


bool mem_canQuregFitInMemory(int numQubits, bool isDensMatr, int numNodes, qindex memBytesPerNode);

bool mem_canMatrixFitInMemory(int numQubits, bool isDense, int numNodes, qindex memBytesPerNode);

bool mem_canSuperOpFitInMemory(int numQubits, qindex numBytesPerNode);



/*
 * MEMORY ALLOCATION SUCCESS
 */


bool mem_isAllocated(int* heapflag);
bool mem_isAllocated(PauliStr* array);
bool mem_isAllocated(qcomp* array);
bool mem_isAllocated(qcomp** matrix, qindex numRows);
bool mem_isAllocated(qcomp*** matrixList, qindex numRows, int numMatrices);

bool mem_isOuterAllocated(qcomp*   ptr);
bool mem_isOuterAllocated(qcomp**  ptr);
bool mem_isOuterAllocated(qcomp*** ptr);



#endif // MEMORY_HPP