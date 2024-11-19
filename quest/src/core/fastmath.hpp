/** @file
 * Oerations used by all deployment modes for fast,
 * low-level maths, inlined and callable within hot
 * loops (i.e OpenMP loops and CUDA kernels)
 */

#ifndef FASTMATH_HPP
#define FASTMATH_HPP

#include "quest/include/precision.h"
#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"



/*
 * ARITHMETIC
 */


INLINE int fast_getPlusOrMinusOne(int isMinus) {

    return 1 - 2 * isMinus;
}



/*
 * INDEX ALGEBRA
 */


INLINE qindex fast_getGlobalRowFromFlatIndex(qindex globalInd, qindex numAmpsPerCol) {

    return globalInd % numAmpsPerCol;
}


INLINE qindex fast_getGlobalColFromFlatIndex(qindex globalInd, qindex numAmpsPerCol) {

    return globalInd / numAmpsPerCol; // floors
}


INLINE qindex fast_getLocalIndexOfDiagonalAmp(qindex localIndOfBasisState, qindex localIndOfFirstDiagAmp, qindex numAmpsPerCol) {

    qindex interDiagSpace = 1 + numAmpsPerCol;
    return localIndOfFirstDiagAmp + (localIndOfBasisState * interDiagSpace);
}


INLINE qindex fast_getLocalFlatIndex(qindex row, qindex localCol, qindex numAmpsPerCol) {

    return row +  localCol * numAmpsPerCol;
}



/*
 * PAULI ALGEBRA
 */


INLINE qcomp fast_getPauliStrCoeff(qindex i, util_pauliStrData data) {

    // (str)|i> = (coeff)|j>
    int par = getBitMaskParity(i & data.allMaskYZ);
    int fac = fast_getPlusOrMinusOne(par);
    return fac * data.powI;
}



#endif // FASTMATH_HPP