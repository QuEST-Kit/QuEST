/** @file
 * Internal functions which query available CPU memory (in an
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

#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/memory.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"

#include <cstdlib>



/*
 * MEMORY BOUNDS
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

int mem_getMaxNumSuperOpQubitsWhichCanFitInMemory(qindex memBytesPerNode) {

    // superoperators have square-bigger superoperators than dense matrices, and are never distributed
    int numNodes = 1;
    bool isDense = true;
    bool hasBuffer = false;
    bool isSuperOp = true;
    return getMaxNumQubitsWhichCanFitInMemory(isDense, numNodes, hasBuffer, isSuperOp, memBytesPerNode);
}


int mem_getMinNumQubitsForDistribution(int numNodes) {

    return logBase2(numNodes);
}


int mem_getMaxNumQuregQubitsBeforeIndexOverflow(bool isDensityMatrix) {

    // cannot store more amplitudes than can be counted by the qindex type (even when distributed)
    qindex maxNumAmps = std::numeric_limits<qindex>::max();
    int maxNumQubits = std::floor(std::log2(maxNumAmps) / static_cast<qreal>((isDensityMatrix)? 2 : 1));
    return maxNumQubits;
}

int mem_getMaxNumMatrixQubitsBeforeIndexOverflow(bool isDenseMatrix) {

    // matrices have the same number of amplitudes as a same-dimension Qureg
    return mem_getMaxNumQuregQubitsBeforeIndexOverflow(isDenseMatrix);
}

int mem_getMaxNumSuperOpQubitsBeforeIndexOverflow() {

    // the superoperator is square-bigger than a dense matrix
    bool isDense = true;
    int maxMatrixQubits = mem_getMaxNumMatrixQubitsBeforeIndexOverflow(isDense);
    int maxSuperOpQubits = maxMatrixQubits / 2; // floors
    return maxSuperOpQubits;
}

qindex mem_getMaxNumKrausMapMatricesBeforeIndexOverflow(int numQubits) {

    qindex numElemPerMatrix = powerOf2(2 * numQubits);
    qindex maxNumTotalElems = std::numeric_limits<qindex>::max();
    qindex maxNumMatrices = maxNumTotalElems / numElemPerMatrix; // floors
    return maxNumMatrices;
}

int getMaxNumQubitsBeforeGlobalMemSizeofOverflow(bool isDensMatr, int numNodes, bool hasBuffer, bool isSuperOp) {

    // this function assumes we must be able to store the total 
    // CPU memory (in bytes) used by a single data structure, 
    // aggregate across all nodes, in a single size_t primitive. 
    // This is a defensively-designed constraint; we do not ever
    // actually need to know the full memory, but assuring that
    // we could principally store it in a size_t will futureproof
    // reporter functions against future overflows etc. Note it
    // does not meaningfully restrict the maximum simulatable size
    // except on ~8 EiB supercomputers. Looking at you, Jupiter!

    size_t maxSizeof = std::numeric_limits<size_t>::max();
    size_t maxGlobalNumAmps = maxSizeof / sizeof(qcomp); // floors
    size_t maxGlobalNumQubits = std::floor(std::log2(maxGlobalNumAmps)); // floors

    // distributing Quregs requires communication buffers, doubling memory, decreasing qubits by 1
    if (hasBuffer && numNodes > 1)
        maxGlobalNumQubits -= 1;

    // density matrices have square-more amps, halving the number of qubtis (AFTER buffer subtraction)
    if (isDensMatr)
        maxGlobalNumQubits /= 2; // floors

    // superoperators are square-bigger than their constituent dense matrices
    if (isSuperOp)
        maxGlobalNumQubits /= 2; // floors

    /// @todo
    /// above sometimes overestimates by one; suggesting N can fit
    /// in fact only N-1 can fit without causing overflows. It's
    /// a chore to correct this precision-agnostically, and we cannot
    /// just try get the total-memory and provoke the overflow because
    /// a call to mem_getLocalQuregMemoryRequired() would recurse! As
    /// this function is only needed for ridiculously overzealous
    /// validation, and does not risk any logical error, we simply
    /// subtract one to avoid the overflowing edge-case. The returned
    /// max-size remains completely unreachable by users of course!
    maxGlobalNumQubits -= 1;

    return maxGlobalNumQubits;
}

int mem_getMaxNumQuregQubitsBeforeGlobalMemSizeofOverflow(bool isDensityMatrix, int numNodes) {

    bool hasBuffer = true;
    bool isSuperOp = false;
    return getMaxNumQubitsBeforeGlobalMemSizeofOverflow(isDensityMatrix, numNodes, hasBuffer, isSuperOp);
}

int mem_getMaxNumMatrixQubitsBeforeGlobalMemSizeofOverflow(bool isDenseMatrix, int numNodes) {

    // matrix types don't store buffers - they'll use those of Quregs they're applied to
    bool hasBuffer = false;
    bool isSuperOp = false;
    return getMaxNumQubitsBeforeGlobalMemSizeofOverflow(isDenseMatrix, numNodes, hasBuffer, isSuperOp);
}

int mem_getMaxNumSuperOpQubitsBeforeGlobalMemSizeofOverflow() {

    // superoperators have square-bigger superoperators than dense matrices, and are never distributed
    int numNodes = 1;
    bool isDense = true;
    bool hasBuffer = false;
    bool isSuperOp = true;
    return getMaxNumQubitsBeforeGlobalMemSizeofOverflow(isDense, numNodes, hasBuffer, isSuperOp);
}

qindex mem_getMaxNumKrausMapMatricesBeforeLocalMemSizeofOverflow(int numQubits) {

    qindex numMatrWithMaxTotalElems = mem_getMaxNumKrausMapMatricesBeforeIndexOverflow(numQubits);
    qindex numMatrWithMaxTotalMem = numMatrWithMaxTotalElems / sizeof(qcomp); // floors
    return numMatrWithMaxTotalMem;
}



/*
 * HARDWARE QUERYING
 */


qindex mem_tryGetLocalRamCapacityInBytes() {

    /// @todo attempt to find total Ram

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

    /// @todo
    ///  if sizeof(qcomp) is a power of 2 (which it almost always is, c'mon now),
    ///  then we could instead return the LOG of the total memory and always
    ///  avoid overflow, permitting reporters to display mem=2^exp.
    ///  it would also make changing units (e.g. to GB) easier.

    // work out individual array costs
    qindex memLocalArray = (qindex) mem_getLocalQuregMemoryRequired(qureg.numAmpsPerNode); // never overflows
    int numLocalArrays = 
        mem_isAllocated(qureg.cpuAmps) + mem_isAllocated(qureg.cpuCommBuffer) +
        mem_isAllocated(qureg.gpuAmps) + mem_isAllocated(qureg.gpuCommBuffer);  // but 4*memLocalArray might overflow

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
    if (numQubits > getMaxNumQubitsBeforeGlobalMemSizeofOverflow(isDenseMatrix, numNodes, hasBuffers, false)) // isSuperop=false
        error_memSizeQueriedButWouldOverflow();

    // no risk of overflow; we have already validated numAmpsTotal fits in qindex
    qindex numAmpsTotal = (isDenseMatrix)? powerOf2(2*numQubits) : powerOf2(numQubits);
    qindex numAmpsPerNode = numAmpsTotal / numNodes; // divides evenly

    // communication buffers double costs
    if (hasBuffers && numNodes > 1)
        numAmpsPerNode *= 2;

    // beware that we must cast to a size_t (which can be greater 
    // than qindex) BEFORE multiplying, to avoid overflows
    return static_cast<size_t>(numAmpsPerNode) *  sizeof(qcomp);
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


size_t mem_getLocalSuperOpMemoryRequired(int numQubits) {

    // superoperators have square-bigger superoperators than dense matrices, and are never distributed
    int numMatrixQubits = 2 * numQubits;
    bool isDense = true;
    int numNodes = 1;
    return mem_getLocalMatrixMemoryRequired(numMatrixQubits, isDense, numNodes);
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


bool mem_canSuperOpFitInMemory(int numQubits, qindex numBytesPerNode) {

    // superoperators are square-bigger than their constituent dense matrices, and are never distributed
    int numMatrixQubits = 2 * numQubits;
    int numNodes = 1;
    bool isDense = true;
    return mem_canMatrixFitInMemory(numMatrixQubits, isDense, numNodes, numBytesPerNode);
}



/*
 * MEMORY ALLOCATION SUCCESS
 *
 * which check that some or all nested pointers are
 * non-NULL, indicating that all allocations involved
 * in a multidimensional data structure were successful.
 * Some of these functions are trivial NULL checks, but
 * are still abstracted here so that data structures can
 * later be adjusted (e.g. CPU matrices may be flattened).
 */


template <typename T>
bool isNonNull(T ptr) {

    // note that (ptr == None) implies (ptr == nullptr)
    return ptr != nullptr;
}


// fast checks which can be used in validation of existing
// heap objects, which check only the outer alloc is valid.
// this is useful for checking whether a user has manually
// modified a heap pointer to be NULL because they are 
// tracking freed objects, but they cannot be used to check
// intiail memory mallocs were successful.

bool mem_isOuterAllocated(qcomp*   ptr) { return isNonNull(ptr); }
bool mem_isOuterAllocated(qcomp**  ptr) { return isNonNull(ptr); }
bool mem_isOuterAllocated(qcomp*** ptr) { return isNonNull(ptr); }


// slow checks that all nested pointers in the heap structure 
// are non-NULL, implying all of them point to valid, existing
// heap memory. This is used by validation after allocation.

bool mem_isAllocated(int* heapflag)   { return isNonNull(heapflag); }
bool mem_isAllocated(qcomp* array )   { return isNonNull(array); }
bool mem_isAllocated(PauliStr* array) { return isNonNull(array); }

bool mem_isAllocated(qcomp** matrix, qindex numRows) {

    if (matrix == nullptr)
        return false;

    // avoid recursing for insignificant speedup
    for (qindex r=0; r<numRows; r++)
        if (matrix[r] == nullptr)
            return false;

    return true;
}

bool mem_isAllocated(qcomp*** matrixList, qindex numRows, int numMatrices) {

    if (matrixList == nullptr)
        return false;

    // fine to recurse because we expect few matrices
    for (qindex n=0; n<numMatrices; n++)
        if (!mem_isAllocated(matrixList[n], numRows))
            return false;

    return true;
}
