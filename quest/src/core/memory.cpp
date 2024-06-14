/** @file
 * Internal functions which query available CPU memory (in an
 * attemptedly OS-agnostic way), and provided needed memory
 * querents. Note GPU memory querying is performed by 
 * the dedicated GPU backend, though this file is always 
 * compiled (even in GPU mode) because GPU-acceleration still 
 * requires accompanying CPU memory arrays.
 */

#include "quest/include/types.h"
#include "quest/include/quest.h"

#include "quest/src/core/memory.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"

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
 * MEMORY COST QUERYING
 */


int mem_getEffectiveNumStateVecQubitsPerNode(int numQubits, bool isDensMatr, int numNodes) {

    // compute logs directly to avoid overflows (even though validation should preclude them)
    qindex logNumAmpsTotal = ((isDensMatr)? 2 : 1) * numQubits;
    qindex logNumAmpsPerNode = logNumAmpsTotal - logBase2(numNodes);
    return logNumAmpsPerNode;
}


int mem_getMinNumQubitsForDistribution(int numNodes) {

    return logBase2(numNodes);
}


int mem_getMaxNumQubitsWhichCanFitInMemory(bool isDensMatr, int numNodes, qindex memBytesPerNode) {

    // distribution requires communication buffers, doubling costs, halving fittable amps-per-qureg
    qindex maxLocalNumAmps = memBytesPerNode / sizeof(qcomp); // floors
    if (numNodes > 1)
        maxLocalNumAmps = maxLocalNumAmps / 2; // floors

    // density matrices require square more memory, so halve (flooring) the number of qubits
    int maxLocalNumQubits = std::floor(std::log2(maxLocalNumAmps));
    if (isDensMatr)
        maxLocalNumQubits /= 2; // floors

    // doubling nodes permits 1 additional qubit
    int maxGlobalNumQubits = maxLocalNumQubits + logBase2(numNodes);
    return maxGlobalNumQubits;
}


bool mem_canQuregFitInMemory(int numQubits, bool isDensMatr, int numNodes, qindex memBytesPerNode) {

    return numQubits <= mem_getMaxNumQubitsWhichCanFitInMemory(isDensMatr, numNodes, memBytesPerNode);
}


int mem_getMaxNumQubitsBeforeIndexOverflow(bool isDensMatr) {

    // cannot store more amplitudes than can be counted by the qindex type (even when distributed)
    qindex maxNumAmps = std::numeric_limits<qindex>::max();
    int maxNumQubits = std::floor(std::log2(maxNumAmps) / (qreal) ((isDensMatr)? 2 : 1));
    return maxNumQubits;
}


int mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(bool isDensMatr, int numNodes) {

    // we return largest N satisfying 2^(2N + [numNodes > 1]) * sizeof(qcomp) / numNodes <= max[sizeof]
    size_t maxSizeof = std::numeric_limits<size_t>::max();
    size_t maxLocalNumAmps = maxSizeof / sizeof(qcomp); // floors
    size_t maxLocalNumQubits = std::floor(std::log2(maxLocalNumAmps));
    size_t maxGlobalNumQubits = maxLocalNumQubits + logBase2(numNodes);

    // distribution requires communication buffers, doubling memory, decreasing qubits by 1
    if (numNodes > 1)
        maxGlobalNumQubits -= 1;

    // density matrices have square-more amps, halving the number of qubtis (AFTER buffer subtraction)
    if (isDensMatr)
        maxGlobalNumQubits = maxGlobalNumQubits / 2; // floors

    return maxGlobalNumQubits;
}


size_t mem_getLocalMemoryRequired(int numQubits, bool isDensMatr, int numNodes) {

    // assert no-overflow precondition
    if (numQubits > mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, numNodes))
        error_memSizeQueriedButWouldOverflow();

    // no risk of overflow; we have already validated numAmpsTotal fits in qindex
    qindex numAmpsTotal = (isDensMatr)? powerOf2(2*numQubits) : powerOf2(numQubits);
    qindex numAmpsPerNode = numAmpsTotal / numNodes; // divides evenly

    // distribution requires communication buffers, doubling costs
    if (numNodes > 1)
        numAmpsPerNode *= 2;

    // return number of bytes to store local amps
    return numAmpsPerNode * sizeof(qcomp);
}


size_t mem_getLocalMemoryRequired(qindex numAmpsPerNode) {

    // assert no-overflow precondition
    qindex maxNumAmpsPerNode = std::numeric_limits<size_t>::max() / sizeof(qcomp); // floors
    if (numAmpsPerNode > maxNumAmpsPerNode)
        error_memSizeQueriedButWouldOverflow();

    // return number of bytes to store local array, EXCLUDING communication buffer
    return numAmpsPerNode * sizeof(qcomp);
}


qindex mem_getTotalGlobalMemoryUsed(Qureg qureg) {

    // TODO:
    //  if sizeof(qcomp) is a power of 2 (which it almost always is, c'mon now),
    //  then we could instead return the LOG of the total memory and always
    //  avoid overflow, permitting reporters to display mem=2^exp.
    //  it would also make changing units (e.g. to GB) easier.

    // work out individual array costs
    size_t memLocalArray = mem_getLocalMemoryRequired(qureg.numAmpsPerNode);
    int numLocalArrays = 
        (qureg.cpuAmps != NULL) + (qureg.cpuCommBuffer != NULL) + 
        (qureg.gpuAmps != NULL) + (qureg.gpuCommBuffer != NULL);

    // if total local costs would overflow qindex, return 0
    qindex maxQindex = std::numeric_limits<qindex>::max();
    qindex maxLocalArrayMem = maxQindex / memLocalArray; // floors
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
