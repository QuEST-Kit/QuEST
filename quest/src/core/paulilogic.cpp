/** @file
 * Internal functions which process Pauli strings
 * and their weighted sums
 * 
 * @author Tyson Jones
 */

#include "quest/include/paulis.h"
#include "quest/include/qureg.h"

#include "quest/src/core/paulilogic.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"

#include <utility>
#include <vector>
#include <array>

using std::vector;



/*
 * PRIVATE UTILITIES
 */


int getPauliFromMaskAt(PAULI_MASK_TYPE mask, int ind) {

    return getTwoAdjacentBits(mask, 2*ind); // bits at (ind+1, ind)
}



/*
 * PauliStr
 */


bool paulis_isIdentity(PauliStr str) {

    return 
        (str.lowPaulis  == 0) && 
        (str.highPaulis == 0);
}


int paulis_getPauliAt(PauliStr str, int ind) {

    return (ind < MAX_NUM_PAULIS_PER_MASK)?
        getPauliFromMaskAt(str.lowPaulis,  ind) :
        getPauliFromMaskAt(str.highPaulis, ind - MAX_NUM_PAULIS_PER_MASK);
}


int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str) {

    int ind   = (str.highPaulis == 0)? 0 : MAX_NUM_PAULIS_PER_MASK;
    auto mask = (str.highPaulis == 0)? str.lowPaulis : str.highPaulis;

    while (mask) {
        mask >>= 2;
        ind++;
    }

    return ind - 1;
}


int paulis_getIndOfLefmostNonIdentityPauli(PauliStr* strings, qindex numStrings) {

    int maxInd = 0;

    for (qindex i=0; i<numStrings; i++) {
        int ind = paulis_getIndOfLefmostNonIdentityPauli(strings[i]);
        if (ind > maxInd)
            maxInd = ind;
    }

    return maxInd;
}


bool paulis_containsXOrY(PauliStr str) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    for (int i=0; i<=maxInd; i++) {
        int pauli = paulis_getPauliAt(str, i);

        if (pauli == 1 || pauli == 2)
            return true;
    }

    return false;
}


int paulis_getSignOfPauliStrConj(PauliStr str) {

    // determine parity of Y count in str
    bool odd = false;
    for (int targ=0; targ < MAX_NUM_PAULIS_PER_STR; targ++) 
        if (paulis_getPauliAt(str, targ) == 2)
            odd = !odd;

    // conj(Y) = -Y, conj(YY) = YY
    return odd? -1 : 1;
}


int paulis_getPrefixZSign(Qureg qureg, vector<int> prefixZ) {

    int sign = 1;

    // each Z contributes +- 1
    for (int qubit : prefixZ)
        sign *= util_getRankBitOfQubit(qubit, qureg)? -1 : 1;

    return sign;
}


qcomp paulis_getPrefixPaulisElem(Qureg qureg, vector<int> prefixY, vector<int> prefixZ) {

    // each Z contributes +- 1
    qcomp elem = paulis_getPrefixZSign(qureg, prefixZ);

    // each Y contributes -+ i
    for (int qubit : prefixY)
        elem *= 1_i * (util_getRankBitOfQubit(qubit, qureg)? 1 : -1);

    return elem;
}


vector<int> paulis_getTargetInds(PauliStr str) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    vector<int> inds(0);
    inds.reserve(maxInd+1);

    for (int i=0; i<=maxInd; i++)
        if (paulis_getPauliAt(str, i) != 0) // Id
            inds.push_back(i);

    return inds;
}


qindex paulis_getTargetBitMask(PauliStr str) {
    
    /// @todo 
    /// would compile-time MAX_NUM_PAULIS_PER_STR bound be faster here,
    /// since this function is invoked upon every PauliStrSum element?
    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    qindex mask = 0;

    for (int i=0; i<=maxInd; i++)
        if (paulis_getPauliAt(str, i) != 0) // Id
            mask = flipBit(mask, i);

    return mask;
}


std::array<vector<int>,3> paulis_getSeparateInds(PauliStr str) {

    vector<int> iXYZ = paulis_getTargetInds(str);
    vector<int> iX, iY, iZ;

    vector<int>* ptrs[] = {&iX, &iY, &iZ};

    for (int i : iXYZ)
        ptrs[paulis_getPauliAt(str, i) - 1]->push_back(i);

    return {iX, iY, iZ};
}


PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift) {

    if (pauliShift <= 0 || pauliShift >= MAX_NUM_PAULIS_PER_MASK)
        error_pauliStrShiftedByIllegalAmount();

    int numBitsPerPauli = 2;
    int numMaskBits = numBitsPerPauli * MAX_NUM_PAULIS_PER_MASK;
    int bitShift    = numBitsPerPauli * pauliShift;

    // record the bits we will lose from lowPaulis, to move to highPaulis
    PAULI_MASK_TYPE lostBits = getBitsLeftOfIndex(str.lowPaulis, numMaskBits - bitShift - 1);

    // ensure we actually lose these bits from lowPaulis
    PAULI_MASK_TYPE lowerBits = getBitsRightOfIndex(str.lowPaulis, numMaskBits - bitShift) << bitShift;

    // and add them to highPaulis; we don't have to force lose upper bits of high paulis
    PAULI_MASK_TYPE upperBits = concatenateBits(str.highPaulis, lostBits, bitShift);

    // return a new stack PauliStr instance (avoiding C++20 initialiser)
    PauliStr out;
    out.lowPaulis = lowerBits;
    out.highPaulis = upperBits;
    return out;
}


PauliStr paulis_getTensorProdOfPauliStr(PauliStr left, PauliStr right, int numQubits) {

    // computes left (tensor) right, assuming right is smaller than numQubits
    PauliStr shifted = paulis_getShiftedPauliStr(left, numQubits);

    // return a new stack PauliStr instance (avoiding C++20 initialiser)
    PauliStr out;
    out.lowPaulis  = right.lowPaulis  | shifted.lowPaulis;
    out.highPaulis = right.highPaulis | shifted.highPaulis;
    return out;
}


PauliStr paulis_getKetAndBraPauliStr(PauliStr str, Qureg qureg) {

    return paulis_getTensorProdOfPauliStr(str, str, qureg.numQubits);
}


PAULI_MASK_TYPE paulis_getKeyOfSameMixedAmpsGroup(PauliStr str) {

    PAULI_MASK_TYPE key = 0;

    // in theory, we can reduce the number of involved operations by bit-shifting
    // str left by 1, XOR'ing this with str, and retaining every 2nd bit, producing
    // e.g. key=0110 from str=IXYZ. However, this is an insignificant speedup which
    // risks sneaky bugs related to handling str's two masks.

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    for (int i=0; i<=maxInd; i++) {
        int pauli = paulis_getPauliAt(str, i);
        int isXY = (pauli == 1 || pauli == 2);
        key |= (isXY << i);
    }

    return key;
}


std::pair<qcomp,PauliStr> paulis_getPauliStrProd(PauliStr strA, PauliStr strB) {

    // a . b = coeff * (a ^ b)
    PauliStr strOut;
    strOut.lowPaulis  = strA.lowPaulis  ^ strB.lowPaulis;
    strOut.highPaulis = strA.highPaulis ^ strB.highPaulis;

    // coeff = product of single-site product coeffs
    qcomp coeff = 1;
    for (int i=0; i<MAX_NUM_PAULIS_PER_STR; i++) {
        int pA = paulis_getPauliAt(strA, i);
        int pB = paulis_getPauliAt(strB, i);
        
        // I.P = P.I = P and P.P = I contribute factor=1
        if (pA == 0 || pB == 0 || pA == pB)
            continue;

        // XY,YZ,ZX=i, XZ,YX,ZY=-i
        int dif = pB - pA;
        coeff *= qcomp(0, (dif == 1 || dif == -2)? 1 : -1);
    }
    
    return {coeff, strOut};
}



/*
 * PauliStrSum
 */


int paulis_getIndOfLefmostNonIdentityPauli(PauliStrSum sum) {

    return paulis_getIndOfLefmostNonIdentityPauli(sum.strings, sum.numTerms);
}


bool paulis_containsXOrY(PauliStrSum sum) {

    for (qindex i=0; i<sum.numTerms; i++)
        if (paulis_containsXOrY(sum.strings[i]))
            return true;

    return false;
}


qindex paulis_getTargetBitMask(PauliStrSum sum) {

    qindex mask = 0;

    // mask has 1 where any str has a != Id
    for (int t=0; t<sum.numTerms; t++)
        mask |= paulis_getTargetBitMask(sum.strings[t]);

    return mask;
}


void paulis_setPauliStrSumToScaledTensorProdOfConjWithSelf(PauliStrSum out, qreal factor, PauliStrSum in, int numQubits) {

    // sets out = factor * conj(in) (x) in, where in has dim of numQubits
    if (paulis_getIndOfLefmostNonIdentityPauli(in) >= numQubits)
        error_pauliStrSumHasMoreQubitsThanSpecifiedInTensorProd();
    if (out.numTerms != in.numTerms * in.numTerms)
        error_pauliStrSumTensorProdHasIncorrectNumTerms();

    // conj(in) (x) in = sum_jk conj(c_j) c_k conj(P_j) (x) P_k...
    qindex i = 0;
    for (qindex j=0; j<in.numTerms; j++) {
        for (qindex k=0; k<in.numTerms; k++) {

            // ... where conj(P_j) = sign_j P_j
            out.strings[i] = paulis_getTensorProdOfPauliStr(in.strings[j], in.strings[k], numQubits);
            out.coeffs[i] = factor * std::conj(in.coeffs[j]) * in.coeffs[k] * paulis_getSignOfPauliStrConj(in.strings[j]);
            i++;
        }
    }
}


qindex paulis_getNumTermsInPauliStrSumProdOfAdjointWithSelf(PauliStrSum in) {

    // adj(in).in has fewer terms than the numTerms^2 bound, since 
    // a.a = I (causing -n and +1 below) and a.b ~ b.a (causing /2);
    // we do not however consider any cancellations of coefficients
    int n = in.numTerms;
    return 1 + (n*n - n)/2;
}


void paulis_setPauliStrSumToScaledProdOfAdjointWithSelf(PauliStrSum out, qreal factor, PauliStrSum in) {

    // sets out = factor * adj(in) . in, permitting duplicate strings
    if (out.numTerms != paulis_getNumTermsInPauliStrSumProdOfAdjointWithSelf(in))
        error_pauliStrSumProdHasIncorrectNumTerms();

    // since out definitely contains an identity (when neglecting coeff cancellation)
    // which is contributed toward by all j=k iterations below, we keep it at i=0
    out.strings[0] = getPauliStr("I");
    out.coeffs[0] = 0;
    qindex i = 1;

    // we leverage that sum_jk a_j^* a_k P_j P_k...
    for (qindex j=0; j<in.numTerms; j++) {

        // = sum_j ( |a_j|^2 Id + sum_k<j ...)
        out.coeffs[0] += factor * std::norm(in.coeffs[j]);

        // containing sum_k<j (a_j^* a_k P_j P_k + a_k^* a_j P_k P_j)
        for (qindex k=0; k<j; k++) {

            // = (a_j^* a_k b_jk + a_k^* a_j b_jk^*) P'
            auto [coeff, str] = paulis_getPauliStrProd(in.strings[j], in.strings[k]);

            // = (x + x^*) P' = 2 Re[x] P'
            out.strings[i] = str;
            out.coeffs[i] = factor * 2 * std::real(std::conj(in.coeffs[j]) * in.coeffs[k] * coeff);
            i++;
        }
    }
}


void paulis_setPauliStrSumToShiftedConj(PauliStrSum out, PauliStrSum in, int numQubits) {

    // sets out = conj(in) (x) I
    if (paulis_getIndOfLefmostNonIdentityPauli(in) >= numQubits)
        error_pauliStrSumHasMoreQubitsThanSpecifiedInConjShift();
    if (out.numTerms != in.numTerms)
        error_pauliStrSumConjHasIncorrectNumTerms();

    // where conj(c P) = conj(c) sign P
    for (qindex i=0; i<out.numTerms; i++) {
        out.strings[i] = paulis_getShiftedPauliStr(in.strings[i], numQubits);
        out.coeffs[i] = std::conj(in.coeffs[i]) * paulis_getSignOfPauliStrConj(in.strings[i]);
    }
}
