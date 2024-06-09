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

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"

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
