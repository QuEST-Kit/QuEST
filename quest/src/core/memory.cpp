/** @file
 * Internal functions which query available CPU memory (in an
 * attemptedly OS-agnostic way) and determine the maximum
 * number of qubits which can be simualted. Note GPU memory
 * querying is performed by the dedicated GPU backend, 
 * though this file is always compiled (even in GPU mode) 
 * because GPU-acceleration still requires accompanying
 * CPU memory arrays.
 */

#include "quest/include/types.h"

#include "quest/src/core/memory.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"

#include <cstdlib>



qindex mem_tryGetLocalRamCapacityInBytes() {

    // TODO:
    //      attempt to find total Ram

    // if we're unable to find total RAM, throw an exception
    // (which the caller should catch and gracefully continue)
    throw (mem::COULD_NOT_QUERY_RAM) false;
}


int mem_getEffectiveNumStateVecQubitsPerNode(int numQubits, bool isDensMatr, int numNodes) {

    // compute logs directly to avoid overflows (even though validation should preclude them)
    qindex logNumAmpsTotal = ((isDensMatr)? 2 : 1) * numQubits;
    qindex logNumAmpsPerNode = logNumAmpsTotal - logBase2(numNodes);
    return logNumAmpsPerNode;
}


bool mem_canQuregFitInMemory(int numQubits, bool isDensMatr, int numNodes, qindex memBytesPerNode) {

    // work in logs to avoid overflows
    qindex logNumAmpsPerNode = mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numNodes);

    // distribution requires communication buffers, doubling costs
    if (numNodes > 1)
        logNumAmpsPerNode += 1;

    // use floating-point division here because sizeof(qcomp) might not be a power-of-2
    qindex maxNumAmpsPerNode = std::floor(memBytesPerNode / (qreal) sizeof(qcomp));
    qindex maxLogNumAmpsPerNode = std::floor(std::log2(maxNumAmpsPerNode));

    // strict inequality because node memory must also fit code pages, etc
    return logNumAmpsPerNode < maxLogNumAmpsPerNode;
}


int mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(bool isDensMatr, int numNodes) {

    // we return largest N satisfying 2^(2N + [numNodes > 1]) * sizeof(qcomp) / numNodes <= max[sizeof]
    size_t maxSizeof = std::numeric_limits<size_t>::max();
    size_t maxLocalNumAmps = std::floor(maxSizeof / (qreal) sizeof(qcomp));
    size_t maxLocalNumQubits = std::floor(std::log2(maxLocalNumAmps));
    size_t maxGlobalNumQubits = maxLocalNumQubits + logBase2(numNodes);

    // distribution requires communication buffers, doubling memory, decreasing qubits by 1
    if (numNodes > 1)
        maxGlobalNumQubits -= 1;

    // density matrices have square-more amps, halving the number of qubtis (AFTER buffer subtraction)
    if (isDensMatr)
        maxGlobalNumQubits = std::floor(maxGlobalNumQubits / 2);

    return maxGlobalNumQubits;
}


size_t mem_getLocalMemoryRequired(int numQubits, bool isDensMatr, int numNodes) {

    // assert no-overflow precondition
    if (numQubits <= mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, numNodes))
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
