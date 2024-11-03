/** @file
 * Miscellaneous operations used by all deployment modes for 
 * fast, common algebra of indices and/or Pauli elements
 */

#ifndef MISCELLANEOUS_HPP
#define MISCELLANEOUS_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"

#include <vector>

using std::vector;



/* 
 * PERFORMANCE-CRITICAL FUNCTIONS
 *
 * which are called in hot loops loops (like by OpenMP threads and
 * CUDA kernels) so are aggressively inlined.
 */


INLINE qindex misc_getLocalIndexOfDiagonalAmp(
    qindex localIndOfBasisState, qindex localIndOfFirstDiagAmp, qindex numAmpsPerCol
) {
    // next diagonal is 1 column across and 1 row down
    qindex interDiagSpace = 1 + numAmpsPerCol;

    return localIndOfFirstDiagAmp + (localIndOfBasisState * interDiagSpace);
}



/*
 * CONVENIENCE FUNCTIONS
 */


qindex misc_getLocalIndexOfFirstDiagonalAmp(Qureg qureg) {

    return qureg.rank * powerOf2(qureg.logNumColsPerNode);
}


qindex misc_getNumLocalDiagonalsWithBits(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    // a corresponding bra-qubit in the prefix with an inconsistent outcome means the node
    // contains no diagonal basis states consistent with the given outcomes
    for (size_t i=0; i<qubits.size(); i++)
        if (!util_isBraQubitInSuffix(qubits[i], qureg))
            if (util_getRankBitOfBraQubit(qubits[i], qureg) != outcomes[i])
                return 0;

    // otherwise, every 2^#qubits local diagonal is consistent with outcomes
    qindex numColsPerNode = powerOf2(qureg.logNumColsPerNode);
    qindex numDiags = numColsPerNode / powerOf2(qubits.size());
    return numDiags;
}



#endif // MISCELLANEOUS_HPP