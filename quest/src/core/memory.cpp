/** @file
 * Internal functions which query available CPU memory (in an
 * attemptedly OS-agnostic way), and provided needed memory
 * querents. Note GPU memory querying is performed by 
 * the dedicated GPU backend, though this file is always 
 * compiled (even in GPU mode) because GPU-acceleration still 
 * requires accompanying CPU memory arrays. This file does not
 * perform any allocation of memory; that is instead performed
 * by cpu_config.cpp, to be symmetric with the GPU-memory
 * allocators in gpu_config.cpp.
 */

#include "quest/include/types.h"

#include "../core/memory.hpp"
#include "../core/bitwise.hpp"
#include "../core/errors.hpp"

#include <cstdlib>



/*
 * HARDWARE QUERYING
 */


qindex mem_tryGetLocalRamCapacityInBytes() {

    // TODO:
    //      attempt to find total Ram

    // if we're unable to find total RAM, throw an exception
    // (which the caller should catch and gracefully continue)
    throw (mem::COULD_NOT_QUERY_RAM) false;
}



/*
 * MEMORY USAGE
 */


int mem_getEffectiveNumStateVecQubitsPerNode(int numQubits, bool isDensMatr, int numNodes) {

    // compute logs directly to avoid overflows (even though validation should preclude them)
    qindex logNumAmpsTotal = ((isDensMatr)? 2 : 1) * numQubits;
    qindex logNumAmpsPerNode = logNumAmpsTotal - logBase2(numNodes);
    return logNumAmpsPerNode;
}


qindex mem_getTotalGlobalMemoryUsed(Qureg qureg) {

    // TODO:
    //  if sizeof(qcomp) is a power of 2 (which it almost always is, c'mon now),
    //  then we could instead return the LOG of the total memory and always
    //  avoid overflow, permitting reporters to display mem=2^exp.
    //  it would also make changing units (e.g. to GB) easier.

    // work out individual array costs
    qindex memLocalArray = (qindex) mem_getLocalQuregMemoryRequired(qureg.numAmpsPerNode); // never overflows
    int numLocalArrays = 
        (qureg.cpuAmps != nullptr) + (qureg.cpuCommBuffer != nullptr) + 
        (qureg.gpuAmps != nullptr) + (qureg.gpuCommBuffer != nullptr);  // but *4 might

    // if total local costs would overflow qindex, return 0
    qindex maxQindex = std::numeric_limits<qindex>::max();
    qindex maxLocalArrayMem = maxQindex / numLocalArrays; // floors
    if (memLocalArray > maxLocalArrayMem)
        return 0;

    // if qureg is non-distributed, compute local CPU+GPU+buffers costs and return
    qindex memLocalTotal = numLocalArrays * memLocalArray;
    if (!qureg.isDistributed)
        return memLocalTotal;

    // else if total global costs would overflow qindex, return 0
    qindex maxLocalTotalMem = maxQindex / qureg.numNodes; // floors
    if (memLocalTotal > maxLocalTotalMem)
        return 0;

    // else compute total costs between all nodes
    qindex memGlobalTotal = memLocalTotal * qureg.numNodes;
    return memGlobalTotal;
}



/*
 * MEMORY REQUIRED
 */


size_t getLocalMemoryRequired(int numQubits, int numNodes, bool isDenseMatrix, bool hasBuffers) {

    // assert no-overflow precondition
    if (numQubits > mem_getMaxNumQuregQubitsBeforeLocalMemSizeofOverflow(isDenseMatrix, numNodes))
        error_memSizeQueriedButWouldOverflow();

    // no risk of overflow; we have already validated numAmpsTotal fits in qindex
    qindex numAmpsTotal = (isDenseMatrix)? powerOf2(2*numQubits) : powerOf2(numQubits);
    qindex numAmpsPerNode = numAmpsTotal / numNodes; // divides evenly

    // communication buffers double costs
    if (hasBuffers && numNodes > 1)
        numAmpsPerNode *= 2;

    // return number of bytes to store local amps
    return numAmpsPerNode * sizeof(qcomp);
}


size_t mem_getLocalQuregMemoryRequired(int numQubits, bool isDensityMatr, int numNodes) {

    // Quregs may need buffers for inter-node communication, depending on numNodes > 1
    bool hasBuffers = true;
    return getLocalMemoryRequired(numQubits, numNodes, isDensityMatr, hasBuffers);
}


size_t mem_getLocalQuregMemoryRequired(qindex numAmpsPerNode) {

    // assert no-overflow precondition
    qindex maxNumAmpsPerNode = std::numeric_limits<size_t>::max() / sizeof(qcomp); // floors
    if (numAmpsPerNode > maxNumAmpsPerNode)
        error_memSizeQueriedButWouldOverflow();

    // return number of bytes to store local array, EXCLUDING communication buffer
    return numAmpsPerNode * sizeof(qcomp);
}


size_t mem_getLocalMatrixMemoryRequired(int numQubits, bool isDenseMatrix, int numNodes) {

    // matrix types don't store buffers - they'll use those of Quregs they're applied to
    bool hasBuffers = false;
    return getLocalMemoryRequired(numQubits, numNodes, isDenseMatrix, hasBuffers);
}



/*
 * QUBIT BOUNDS
 */


int getMaxNumQubitsWhichCanFitInMemory(bool isDensMatr, int numNodes, bool hasBuffer, bool isSuperOp, qindex memBytesPerNode) {

    // distribution requires communication buffers, doubling costs, halving fittable amps-per-qureg
    qindex maxLocalNumAmps = memBytesPerNode / sizeof(qcomp); // floors
    if (hasBuffer && numNodes > 1)
        maxLocalNumAmps /= 2; // floors

    // density matrices require square more memory, so halve (flooring) the number of qubits
    int maxLocalNumQubits = std::floor(std::log2(maxLocalNumAmps));
    if (isDensMatr)
        maxLocalNumQubits /= 2; // floors

    // superoperators require square more memory still, so halve (flooring) the number of qubits
    if (isSuperOp)
        maxLocalNumQubits /= 2; // floors

    // doubling nodes permits 1 additional qubit
    int maxGlobalNumQubits = maxLocalNumQubits + logBase2(numNodes);
    return maxGlobalNumQubits;
}


int mem_getMaxNumQuregQubitsWhichCanFitInMemory(bool isDensityMatrix, int numNodes, qindex memBytesPerNode) {

    bool hasBuffer = true;
    bool isSuperOp = false;
    return getMaxNumQubitsWhichCanFitInMemory(isDensityMatrix, numNodes, hasBuffer, isSuperOp, memBytesPerNode);
}

int mem_getMaxNumMatrixQubitsWhichCanFitInMemory(bool isDenseMatrix, int numNodes, qindex memBytesPerNode) {

    // matrix types don't store buffers - they'll use those of Quregs they're applied to
    bool hasBuffer = false;
    bool isSuperOp = false;
    return getMaxNumQubitsWhichCanFitInMemory(isDenseMatrix, numNodes, hasBuffer, isSuperOp, memBytesPerNode);
}


int mem_getMinNumQubitsForDistribution(int numNodes) {

    return logBase2(numNodes);
}


int mem_getMaxNumQuregQubitsBeforeIndexOverflow(bool isDensityMatrix) {

    // cannot store more amplitudes than can be counted by the qindex type (even when distributed)
    qindex maxNumAmps = std::numeric_limits<qindex>::max();
    int maxNumQubits = std::floor(std::log2(maxNumAmps) / (qreal) ((isDensityMatrix)? 2 : 1));
    return maxNumQubits;
}

int mem_getMaxNumMatrixQubitsBeforeIndexOverflow(bool isDenseMatrix) {

    // matrices have the same number of amplitudes as a same-dimension Qureg
    return mem_getMaxNumQuregQubitsBeforeIndexOverflow(isDenseMatrix);
}


int getMaxNumQubitsBeforeLocalMemSizeofOverflow(bool isDensMatr, int numNodes, bool hasBuffer, bool isSuperOp) {

    // we return largest N satisfying 2^(2N + [numNodes > 1]) * sizeof(qcomp) / numNodes <= max[sizeof]
    size_t maxSizeof = std::numeric_limits<size_t>::max();
    size_t maxLocalNumAmps = maxSizeof / sizeof(qcomp); // floors
    size_t maxLocalNumQubits = std::floor(std::log2(maxLocalNumAmps));
    size_t maxGlobalNumQubits = maxLocalNumQubits + logBase2(numNodes);

    // distributing Quregs requires communication buffers, doubling memory, decreasing qubits by 1
    if (hasBuffer && numNodes > 1)
        maxGlobalNumQubits -= 1;

    // density matrices have square-more amps, halving the number of qubtis (AFTER buffer subtraction)
    if (isDensMatr)
        maxGlobalNumQubits /= 2; // floors

    // superoperators are square-bigger than their constituent dense matrices
    if (isSuperOp)
        maxGlobalNumQubits /= 2; // floors

    return maxGlobalNumQubits;
}

int mem_getMaxNumQuregQubitsBeforeLocalMemSizeofOverflow(bool isDensityMatrix, int numNodes) {

    bool hasBuffer = true;
    bool isSuperOp = false;
    return getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensityMatrix, numNodes, hasBuffer, isSuperOp);
}

int mem_getMaxNumMatrixQubitsBeforeLocalMemSizeofOverflow(bool isDenseMatrix, int numNodes) {

    // matrix types don't store buffers - they'll use those of Quregs they're applied to
    bool hasBuffer = false;
    bool isSuperOp = false;
    return getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDenseMatrix, numNodes, hasBuffer, isSuperOp);
}



/*
 * SUFFICIENT MEMORY QUERYING
 */


bool mem_canQuregFitInMemory(int numQubits, bool isDensMatr, int numNodes, qindex memBytesPerNode) {

    return numQubits <= mem_getMaxNumQuregQubitsWhichCanFitInMemory(isDensMatr, numNodes, memBytesPerNode);
}


bool mem_canMatrixFitInMemory(int numQubits, bool isDense, int numNodes, qindex memBytesPerNode) {

    // this function's logic is similar to mem_canQuregFitInMemory(), where diagonal matrices are
    // like statevectors and dense matrices are like density-matrices, except that distributed
    // matrices (numNodes > 1) do not store (nor need to account for) communication buffers

    // distributing the matrix shrinks the local number of qubits stored, effectively
    int localNumQubits = numQubits - logBase2(numNodes);

    // work out the maximum "local" qubits that can fit in memory
    qindex maxLocalNumElems = memBytesPerNode / sizeof(qcomp); // floors
    int maxLocalNumQubits  = std::floor(std::log2(maxLocalNumElems));

    // dense matrices (as opposed to diagonals) require square more memory
    if (isDense)
        maxLocalNumQubits /= 2; // floors

    return localNumQubits <= maxLocalNumQubits;
}
