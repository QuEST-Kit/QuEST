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



/*
 * ARITHMETIC
 */


INLINE int fast_getPlusOrMinusOne(int isMinus) {

    return 1 - 2 * isMinus;
}


INLINE int fast_getPlusOrMinusMaskedBitParity(qindex num, qindex mask) {

    qindex bits = num & mask;
    int pari = getBitMaskParity(bits);
    int sign = fast_getPlusOrMinusOne(pari);
    return sign;
}



/*
 * INDEX ALGEBRA
 */


INLINE qindex fast_getGlobalRowFromFlatIndex(qindex localOrGlobalInd, qindex numAmpsPerCol) {

    return localOrGlobalInd % numAmpsPerCol;
}


INLINE qindex fast_getGlobalColFromFlatIndex(qindex globalInd, qindex numAmpsPerCol) {

    return globalInd / numAmpsPerCol; // floors
}


INLINE qindex fast_getLocalIndexOfDiagonalAmp(qindex localIndOfBasisState, qindex localIndOfFirstDiagAmp, qindex numAmpsPerCol) {

    qindex interDiagSpace = 1 + numAmpsPerCol; // constant and optimised away
    return localIndOfFirstDiagAmp + (localIndOfBasisState * interDiagSpace);
}


INLINE qindex fast_getLocalFlatIndex(qindex row, qindex localCol, qindex numAmpsPerCol) {

    return row +  localCol * numAmpsPerCol;
}



/*
 * PAULI ALGEBRA
 */


// 'qcomp' cannot be used inside CUDA kernels/thrust, so below functions are 
// parsed as cu_qcomp in the GPU backend, which will prior define USE_CU_QCOMP
#ifdef USE_CU_QCOMP
    #define QCOMP_ALIAS cu_qcomp
#else
    #define QCOMP_ALIAS qcomp
#endif


INLINE QCOMP_ALIAS fast_getLowerPauliStrElem(PauliStr str, qindex row, qindex col) {

    // this function is deliberately NOT named getPauliStrElem():
    // the size of PauliStr's LOWER mask (of the two) is 64 bits = 32 Paulis.
    // this function is only ever called to populate a density matrix, which
    // will assuredly contain fewer than 32 qubits (=64 qubit statevector), so
    // we actually do not ever need to consult the HIGHER mas. Ergo, the
    // precondition is that the str.highPaulis == 0

    // regrettably duplicated from paulis.cpp, inaccessible here
    constexpr int MAX_NUM_PAULIS_PER_MASK = sizeof(PAULI_MASK_TYPE) * 8 / 2;

    static const QCOMP_ALIAS matrices[][2][2] = {
        {{1,0},{0,1}},
        {{0,1},{1,0}},
        {{0,-1_i},{1_i,0}},
        {{1,0},{0,-1}}};

    QCOMP_ALIAS elem = 1;

    // could be compile-time unrolled into 32 iterations
    for (int t=0; t<MAX_NUM_PAULIS_PER_MASK; t++) {
        int p = getTwoAdjacentBits(str.lowPaulis, 2*t);
        int i = getBit(row, t);
        int j = getBit(col, t);
        elem *= matrices[p][i][j];
    }

    // to crudely safe-guard against the erroneous scenario where str.highPaulis!=0,
    // without compromising the inline performance, we will efficiently sabotage the
    // result so that the resulting bug is not insidious. We'd prefer to multiply
    // NaN (assuming 0 * NaN = 0) but it there is no platform agnostic efficient literal
    elem *= 1 + str.highPaulis * 1E200;

    return elem;
}


INLINE QCOMP_ALIAS fast_getLowerPauliStrSumElem(QCOMP_ALIAS* coeffs, PauliStr* strings, qindex numTerms, qindex row, qindex col) {

    // this function accepts unpacked PauliStrSum fields since a PauliStrSum cannot 
    // be directly processed in CUDA kernels/thrust due to its 'qcomp' field.
    // it also assumes str.highPaulis==0 for all str in strings, as per above func.

    QCOMP_ALIAS elem = 0;

    // this loop is expected exponentially smaller than caller's loop
    for (qindex n=0; n<numTerms; n++)
        elem += coeffs[n] * fast_getLowerPauliStrElem(strings[n], row, col);

    return elem;
}


// avoid exposing temporary macro
#undef QCOMP_ALIAS



#endif // FASTMATH_HPP