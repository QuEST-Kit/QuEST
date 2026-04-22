/** @file
 * Oerations used by all deployment modes for fast,
 * low-level maths, inlined and callable within hot
 * loops (i.e OpenMP loops and CUDA kernels).
 * 
 * Note this file uses the backend OpenMP/CUDA-agnostic
 * base_qcomp, in lieu of the C++-user-facing qcomp (i.e.
 * std::complex), to avoid its compiler-specific performance
 * pitfalls.
 * 
 * @author Tyson Jones
 */

#ifndef FASTMATH_HPP
#define FASTMATH_HPP

#include "quest/include/precision.h"
#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/base_qcomp.hpp"



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


INLINE qindex fast_getQuregGlobalRowFromFlatIndex(qindex localOrGlobalInd, qindex numAmpsPerCol) {

    return localOrGlobalInd % numAmpsPerCol;
}


INLINE qindex fast_getQuregGlobalColFromFlatIndex(qindex globalInd, qindex numAmpsPerCol) {

    return globalInd / numAmpsPerCol; // floors
}


INLINE qindex fast_getQuregLocalIndexOfDiagonalAmp(qindex localIndOfBasisState, qindex localIndOfFirstDiagAmp, qindex numAmpsPerCol) {

    qindex interDiagSpace = 1 + numAmpsPerCol; // constant and optimised away
    return localIndOfFirstDiagAmp + (localIndOfBasisState * interDiagSpace);
}


INLINE qindex fast_getQuregLocalFlatIndex(qindex row, qindex localCol, qindex numAmpsPerCol) {

    // Qureg density matrices are column-major
    return row + localCol * numAmpsPerCol;
}

INLINE qindex fast_getMatrixFlatIndex(qindex row, qindex col, qindex numAmpsPerCol) {

    // non-distributed matrices are row-major
    return col + row * numAmpsPerCol;
}


INLINE void fast_getSubQuregValues(qindex basisStateIndex, int* numQubitsPerSubQureg, int numSubQuregs, bool areSigned, qindex* outValues) {

    qindex remainingValue = basisStateIndex;

    // find unsigned integer value of each var by partitioning index bits between vars (never unrolled)
    for (int n=0; n<numSubQuregs; n++) {
        outValues[n] =  getBitsRightOfIndex(remainingValue, numQubitsPerSubQureg[n]);
        remainingValue = getBitsLeftOfIndex(remainingValue, numQubitsPerSubQureg[n]-1);
    }

    // two's-complement signed integers are negated if the leftmost variable sign bit is 1 
    if (areSigned)
        for (int v=0; v<numSubQuregs; v++)
            if (getBit(outValues[v], numQubitsPerSubQureg[v]-1))
                outValues[v] -= powerOf2(numQubitsPerSubQureg[v] - 1);
}



/*
 * PAULI ALGEBRA
 */


INLINE base_qcomp fast_getPauliStrElem(PauliStr str, qindex row, qindex col) {

    // this function is called by both fullstatediagmatr_setElemsToPauliStrSum()
    // and densmatr_setAmpsToPauliStrSum_sub(). The former's PauliStr can have
    // Paulis on any of the 64 sites, but the latter's PauliStr is always
    // constrainted to the lower 32 sites (because a 32-qubit density matrix
    // is already too large for the world's computers). As such, the latter
    // scenario can be optimised since str.highPaulis == 0, making the second
    // loop below redundant. Avoiding this loop can at most half the runtime,
    // though opens the risk that the former caller erroneously has its upper
    // Paulis ignore. We forego this optimisation in defensive design, and
    // because this function is only invoked during data structure initilisation
    // and ergo infrequently.

    // regrettably duplicated from paulis.cpp which is inaccessible here
    constexpr int numPaulisPerMask = sizeof(PAULI_MASK_TYPE) * 8 / 2;

    // T-agnostic complex literals
    base_qcomp p0, p1,n1, pI,nI;
    p0 = {0,  0}; //  0
    p1 = {+1, 0}; //  1
    n1 = {-1, 0}; // -1
    pI = {0, +1}; //  i
    nI = {0, -1}; // -i

    // 'matrices' below is not declared constexpr or static const, even though
    // it is fixed/known at compile-time, because this makes it incompatible
    // with CUDA kernels/thrust. It is instead left as runtime innitialisation
    // but this poses no real slowdown; this function, and its caller, are inlined
    // so these 16 amps are re-processed one for each full enumeration of the
    // PauliStrSum which is expected to have significantly more terms/coeffs
    base_qcomp matrices[][2][2] = {
        {{p1,p0},{p0,p1}},  // I
        {{p0,p1},{p1,p0}},  // X
        {{p0,nI},{pI,p0}},  // Y
        {{p1,p0},{p0,n1}}}; // Z

    base_qcomp elem = p1; // 1

    // could be compile-time unrolled into 32 iterations
    for (int t=0; t<numPaulisPerMask; t++) {
        int p = getTwoAdjacentBits(str.lowPaulis, 2*t);
        int i = getBit(row, t);
        int j = getBit(col, t);
        elem *= matrices[p][i][j];
    }

    // could be compile-time unrolled into 32 iterations
    for (int t=0; t<numPaulisPerMask; t++) {
        int p = getTwoAdjacentBits(str.highPaulis, 2*t);
        int i = getBit(row, t + numPaulisPerMask);
        int j = getBit(col, t + numPaulisPerMask);
        elem *= matrices[p][i][j];
    }

    return elem;
}


INLINE base_qcomp fast_getPauliStrSumElem(base_qcomp* coeffs, PauliStr* strings, qindex numTerms, qindex row, qindex col) {

    // this function accepts unpacked PauliStrSum fields since a PauliStrSum cannot 
    // be directly processed in CUDA kernels/thrust due to its 'qcomp' field.
    // it also assumes str.highPaulis==0 for all str in strings, as per above func.

    base_qcomp elem = {0, 0};

    // this loop is expected exponentially smaller than caller's loop
    for (qindex n=0; n<numTerms; n++)
        elem += coeffs[n] * fast_getPauliStrElem(strings[n], row, col);

    return elem;
}


#endif // FASTMATH_HPP