/** @file
 * Miscellaneous utility functions needed internally,
 * which are not performance critical (i.e. not used
 * in hot loops). These include qubit and state index
 * logic, matrix algebra, and channel parameters.
 * 
 * @author Tyson Jones
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"
#include "quest/include/precision.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <functional>
#include <algorithm>
#include <complex>
#include <cmath>
#include <vector>
#include <array>
#include <list>
#include <new>

using std::vector;



/*
 * QUBIT PROCESSING
 */

int util_getPrefixInd(int qubit, Qureg qureg) {
    if (qubit < qureg.logNumAmpsPerNode)
        error_utilsGetPrefixIndGivenSuffixQubit();

    return qubit - qureg.logNumAmpsPerNode;
}

int util_getBraQubit(int ketQubit, Qureg qureg) {
    if (!qureg.isDensityMatrix)
        error_utilsGetBraIndGivenNonDensMatr();

    return ketQubit + qureg.numQubits;
}

int util_getPrefixBraInd(int ketQubit, Qureg qureg) {
    if (!qureg.isDensityMatrix)
        error_utilsGetPrefixBraIndGivenNonDensMatr();
    if (ketQubit < qureg.logNumColsPerNode)
        error_utilsGetPrefixBraIndGivenSuffixQubit();
    
    // equivalent to util_getPrefixInd of util_getBraQubit
    return ketQubit - qureg.logNumColsPerNode;
}

bool util_isQubitInSuffix(int qubit, Qureg qureg) {

    return qubit < qureg.logNumAmpsPerNode;
}

bool util_areAllQubitsInSuffix(vector<int> qubits, Qureg qureg) {

    for (int q : qubits)
        if (!util_isQubitInSuffix(q, qureg))
            return false;
        
    return true;
}

bool util_isBraQubitInSuffix(int ketQubit, Qureg qureg) {
    if (!qureg.isDensityMatrix)
        error_utilsIsBraQubitInSuffixGivenNonDensMatr();

    return ketQubit < qureg.logNumColsPerNode;
}

vector<int> getPrefixOrSuffixQubits(vector<int> qubits, Qureg qureg, bool getSuffix) {

    // note that when the qureg is local/duplicated, 
    // all qubits will be suffix, none will be prefix

    vector<int> subQubits(0);
    subQubits.reserve(qubits.size());

    for (int qubit : qubits)
        if (util_isQubitInSuffix(qubit, qureg) == getSuffix)
            subQubits.push_back(qubit);

    return subQubits;
}

std::array<vector<int>,2> util_getPrefixAndSuffixQubits(vector<int> qubits, Qureg qureg) {
    return {
        getPrefixOrSuffixQubits(qubits, qureg, false), 
        getPrefixOrSuffixQubits(qubits, qureg, true)
    };
}

int util_getRankBitOfQubit(int ketQubit, Qureg qureg) {

    int rankInd = util_getPrefixInd(ketQubit, qureg);
    int rankBit = getBit(qureg.rank, rankInd);
    return rankBit;
}

int util_getRankBitOfBraQubit(int ketQubit, Qureg qureg) {
    
    int rankInd = util_getPrefixBraInd(ketQubit, qureg);
    int rankBit = getBit(qureg.rank, rankInd);
    return rankBit;
}

int util_getRankWithQubitFlipped(int prefixKetQubit, Qureg qureg) {

    int rankInd = util_getPrefixInd(prefixKetQubit, qureg);
    int rankFlip = flipBit(qureg.rank, rankInd);
    return rankFlip;
}

int util_getRankWithQubitsFlipped(vector<int> prefixQubits,  Qureg qureg) {

    int rank = qureg.rank;
    for (int qubit : prefixQubits)
        rank = flipBit(rank, util_getPrefixInd(qubit, qureg));

    return rank;
}

int util_getRankWithBraQubitFlipped(int ketQubit, Qureg qureg) {

    int rankInd = util_getPrefixBraInd(ketQubit, qureg);
    int rankFlip = flipBit(qureg.rank, rankInd);
    return rankFlip;
}

int util_getRankWithBraQubitsFlipped(vector<int> ketQubits, Qureg qureg) {

    int rank = qureg.rank;
    for (int qubit : ketQubits)
        rank = flipBit(rank, util_getPrefixBraInd(qubit, qureg));

    return rank;
}

vector<int> util_getBraQubits(vector<int> ketQubits, Qureg qureg) {

    vector<int> braInds(0);
    braInds.reserve(ketQubits.size());

    for (int qubit : ketQubits)
        braInds.push_back(util_getBraQubit(qubit, qureg));

    return braInds;
}

vector<int> util_getNonTargetedQubits(int* targets, int numTargets, int numQubits) {
    
    qindex mask = getBitMask(targets, numTargets);

    vector<int> nonTargets;
    nonTargets.reserve(numQubits - numTargets);

    for (int i=0; i<numQubits; i++)
        if (getBit(mask, i) == 0)
            nonTargets.push_back(i);

    return nonTargets;
}

vector<int> util_getConcatenated(vector<int> list1, vector<int> list2) {

    // modify the copy of list1
    list1.insert(list1.end(), list2.begin(), list2.end());
    return list1;
}

vector<int> util_getSorted(vector<int> qubits) {

    vector<int> copy = qubits;
    std::sort(copy.begin(), copy.end());
    return copy;
}

vector<int> util_getSorted(vector<int> ctrls, vector<int> targs) {

    return util_getSorted(util_getConcatenated(ctrls, targs));
}

qindex util_getBitMask(vector<int> qubits) {

    // inserts qubits in state 1
    return getBitMask(qubits.data(), qubits.size());
}

qindex util_getBitMask(vector<int> qubits, vector<int> states) {

    return getBitMask(qubits.data(), states.data(), states.size());
}

qindex util_getBitMask(vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, vector<int> targStates) {

    auto qubits = util_getConcatenated(ctrls, targs);
    auto states = util_getConcatenated(ctrlStates, targStates);
    return util_getBitMask(qubits, states);
}

vector<int> util_getVector(int* qubits, int numQubits) {

    // permit qubits=nullptr, overriding numQubits (might be non-zero)
    if (qubits == nullptr)
        return {};

    return vector<int> (qubits, qubits + numQubits);
}



/*
 * INDEX ALGEBRA
 */

qindex util_getGlobalIndexOfFirstLocalAmp(Qureg qureg) {

    return qureg.rank * qureg.numAmpsPerNode;
}

qindex util_getGlobalColumnOfFirstLocalAmp(Qureg qureg) {
    assert_utilsGivenDensMatr(qureg);

    return qureg.rank * powerOf2(qureg.logNumColsPerNode);
}

qindex util_getLocalIndexOfGlobalIndex(Qureg qureg, qindex globalInd) {

    // equivalent to below, but clearer
    if (!qureg.isDistributed)
        return globalInd;

    // defensive-design integrity check
    qindex globalStart = util_getGlobalIndexOfFirstLocalAmp(qureg);
    qindex globalEnd = globalStart + qureg.numAmpsPerNode;
    if (globalInd < globalStart || globalInd >= globalEnd)
        error_utilsGivenGlobalIndexOutsideNode();

    return globalInd % qureg.numAmpsPerNode;
}

qindex util_getLocalIndexOfFirstDiagonalAmp(Qureg qureg) {
    assert_utilsGivenDensMatr(qureg);

    return qureg.rank * powerOf2(qureg.logNumColsPerNode);
}

qindex util_getGlobalFlatIndex(Qureg qureg, qindex globalRow, qindex globalCol) {
    assert_utilsGivenDensMatr(qureg);

    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    return (globalCol * numAmpsPerCol) + globalRow;
}

int util_getRankContainingIndex(Qureg qureg, qindex globalInd) {

    // when not distributed, root contains the index (as incidentally do all nodes)
    if (!qureg.isDistributed)
        return ROOT_RANK;

    // accepts flat density matrix index too
    return globalInd / qureg.numAmpsPerNode; // floors
}
int util_getRankContainingIndex(FullStateDiagMatr matr, qindex globalInd) {

    // when not distributed, root contains the index (as incidentally do all nodes)
    if (!matr.isDistributed)
        return ROOT_RANK;

    return globalInd / matr.numElemsPerNode; // floors
}

int util_getRankContainingColumn(Qureg qureg, qindex globalCol) {
    assert_utilsGivenDensMatr(qureg);

    /// when not distributed, root contains the column (as incidentally do all nodes)
    if (!qureg.isDistributed)
        return ROOT_RANK;

    qindex numColsPerNode = powerOf2(qureg.logNumColsPerNode);
    return globalCol / numColsPerNode; // floors
}

qindex util_getNextPowerOf2(qindex number) {

    int nextExponent = static_cast<int>(std::ceil(std::log2(number)));
    return powerOf2(nextExponent);
}

qcomp util_getElemFromNestedPtrs(void* in, qindex* inds, int numInds) {
    qindex ind = inds[0];

    if (numInds == 1)
        return ((qcomp*) in)[ind];

    qcomp* ptr = ((qcomp**) in)[ind];
    return util_getElemFromNestedPtrs(ptr, &inds[1], numInds-1); // compiler may optimise tail-recursion
}



/*
 * SCALAR ALGEBRA
 */

bool util_isStrictlyInteger(qreal num) {

    return std::trunc(num) == num;
}

bool util_isApproxReal(qcomp num, qreal eps) {

    return std::abs(std::imag(num)) <= eps;
}

qcomp util_getPowerOfI(size_t exponent) {

    // seems silly, but at least it's precision agnostic!
    qcomp values[] = {1, 1_i, -1, -1_i};
    return values[exponent % 4];
}



/*
 * VECTOR REDUCTION
 */

qreal util_getSum(vector<qreal> list) {

    qreal sum = 0;
    qreal y, t, c=0;
    
    // complex Kahan summation
    for (auto& x : list) {
        y = x - c;
        t = sum + y;
        c = ( t - sum ) - y;
        sum = t;
    }

    return sum;
}



/*
 * MATRIX CONJUGATION
 */

// type T can be qcomp** or qcomp*[]
template <typename T>
void setDenseElemsConj(T elems, qindex dim) {
    for (qindex i=0; i<dim; i++)
        for (qindex j=0; j<dim; j++)
           elems[i][j] = std::conj(elems[i][j]);
}

// diagonals don't need templating because arrays decay to pointers, yay!
void setDiagElemsConj(qcomp* elems, qindex dim) {
    for (qindex i=0; i<dim; i++)
        elems[i] = std::conj(elems[i]);
}

CompMatr1 util_getConj(CompMatr1 matrix) {
    CompMatr1 conj = matrix;
    setDenseElemsConj(conj.elems, matrix.numRows);
    return conj;
}
CompMatr2 util_getConj(CompMatr2 matrix) {
    CompMatr2 conj = matrix;
    setDenseElemsConj(conj.elems, matrix.numRows);
    return conj;
}

DiagMatr1 util_getConj(DiagMatr1 matrix) {
    DiagMatr1 conj = matrix;
    setDiagElemsConj(conj.elems, matrix.numElems);
    return conj;
}
DiagMatr2 util_getConj(DiagMatr2 matrix) {
    DiagMatr2 conj = matrix;
    setDiagElemsConj(conj.elems, matrix.numElems);
    return conj;
}

void util_setConj(CompMatr matrix) {
    setDenseElemsConj(matrix.cpuElems, matrix.numRows);
}
void util_setConj(DiagMatr matrix) {
    setDiagElemsConj(matrix.cpuElems, matrix.numElems);
}



/*
 * MATRIX UNITARITY
 */

bool isApprox(qreal a, qreal b, qreal eps) {
    return std::abs(a-b) <= eps;
}
bool isApprox(qcomp a, qcomp b, qreal eps) {
    return std::norm(a-b) <= eps;
}

// type T can be qcomp** or qcomp*[]
template <typename T>
bool getUnitarity(T elems, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    /// @todo
    /// consider multithreading or GPU-accelerating this
    /// when caller is big and e.g. has GPU memory

    // check m * dagger(m) == identity
    for (qindex r=0; r<dim; r++) {
        for (qindex c=0; c<dim; c++) {

            // compute m[r,...] * dagger(m)[...,c]
            qcomp elem = 0;
            for (qindex i=0; i<dim; i++)
                elem += elems[r][i] * std::conj(elems[c][i]);

            // check if further than epsilon from identity[r,c]
            if (!isApprox(elem, qcomp(r==c,0), eps))
                return false;
        }
    }

    return true;
}

// diagonal version doesn't need templating because array decays to pointer, yay!
bool getUnitarity(qcomp* diags, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    /// @todo
    /// consider multithreading or GPU-accelerating this
    /// when caller is big and e.g. has GPU memory

    // check every element has unit magnitude
    for (qindex i=0; i<dim; i++)
        if (!isApprox(std::abs(diags[i]), 1, eps))
            return false;

    return true;
}

// unitarity of fixed-size matrices is always computed afresh
bool util_isUnitary(CompMatr1 m, qreal eps) { return getUnitarity(m.elems, m.numRows,  eps); }
bool util_isUnitary(CompMatr2 m, qreal eps) { return getUnitarity(m.elems, m.numRows,  eps); }
bool util_isUnitary(DiagMatr1 m, qreal eps) { return getUnitarity(m.elems, m.numElems, eps); }
bool util_isUnitary(DiagMatr2 m, qreal eps) { return getUnitarity(m.elems, m.numElems, eps); }

// unitarity of heap matrices is cached
bool util_isUnitary(CompMatr m, qreal eps) {

    // compute and record unitarity if not already known
    if (*(m.isApproxUnitary) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isApproxUnitary) = getUnitarity(m.cpuElems, m.numRows, eps);

    // eps may have been ignored
    return *(m.isApproxUnitary);
}
bool util_isUnitary(DiagMatr m, qreal eps) {

    // compute and record unitarity if not already known
    if (*(m.isApproxUnitary) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isApproxUnitary) = getUnitarity(m.cpuElems, m.numElems, eps);

    // eps may have been ignored
    return *(m.isApproxUnitary);
}
bool util_isUnitary(FullStateDiagMatr m, qreal eps) {

    // compute and record unitarity if not already known
    if (*(m.isApproxUnitary) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isApproxUnitary) = getUnitarity(m.cpuElems, m.numElemsPerNode, eps);

    // communication may be necessary
    if (m.isDistributed)
        *(m.isApproxUnitary) = comm_isTrueOnAllNodes(*(m.isApproxUnitary));

    // eps may have been ignored
    return *(m.isApproxUnitary);
}



/*
 * MATRIX HERMITICITY
 */

// type T can be qcomp** or qcomp*[]
template <typename T>
bool getHermiticity(T elems, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    /// @todo
    /// consider multithreading or GPU-accelerating this
    /// when caller is big and e.g. has GPU memory

    // check adjoint(elems) == elems
    for (qindex r=0; r<dim; r++)
        for (qindex c=0; c<r; c++)
            if (!isApprox(elems[r][c], std::conj(elems[c][r]), eps))
                return false;

    return true;
}

// diagonal version doesn't need templating because array decays to pointer, yay!
bool getHermiticity(qcomp* diags, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    /// @todo
    /// consider multithreading or GPU-accelerating this
    /// when caller is big and e.g. has GPU memory

    // check every element has a zero (or <eps) imaginary component
    for (qindex i=0; i<dim; i++)
        if (!isApprox(std::imag(diags[i]), 0, eps))
            return false;

    return true;
}

// hermiticity of fixed-size matrices is always computed afresh
bool util_isHermitian(CompMatr1 m, qreal eps) { return getHermiticity(m.elems, m.numRows,  eps); }
bool util_isHermitian(CompMatr2 m, qreal eps) { return getHermiticity(m.elems, m.numRows,  eps); }
bool util_isHermitian(DiagMatr1 m, qreal eps) { return getHermiticity(m.elems, m.numElems, eps); }
bool util_isHermitian(DiagMatr2 m, qreal eps) { return getHermiticity(m.elems, m.numElems, eps); }

// hermiticity of heap matrices is cached
bool util_isHermitian(CompMatr m, qreal eps) {

    // compute and record hermiticity if not already known
    if (*(m.isApproxHermitian) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isApproxHermitian) = getHermiticity(m.cpuElems, m.numRows, eps);

    // eps may have been ignored
    return *(m.isApproxHermitian);
}
bool util_isHermitian(DiagMatr m, qreal eps) {

    // compute and record hermiticity if not already known
    if (*(m.isApproxHermitian) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isApproxHermitian) = getHermiticity(m.cpuElems, m.numElems, eps);

    // eps may have been ignored
    return *(m.isApproxHermitian);
}
bool util_isHermitian(FullStateDiagMatr m, qreal eps) {

    // compute and record hermiticity if not already known
    if (*(m.isApproxHermitian) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isApproxHermitian) = getHermiticity(m.cpuElems, m.numElemsPerNode, eps);

    // communication may be necessary
    if (m.isDistributed)
        *(m.isApproxHermitian) = comm_isTrueOnAllNodes(*(m.isApproxHermitian));

    // eps may have been ignored
    return *(m.isApproxHermitian);
}



/*
 * EXPONENTIABLE MATRIX IS NON-ZERO
 */

bool getWhetherNonZero(qcomp* diags, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    /// @todo
    /// consider multithreading or GPU-accelerating this
    /// when caller is big and e.g. has GPU memory

    for (qindex i=0; i<dim; i++) {

        // check each complex element has non-zero abs
        if (isApprox(std::abs(diags[i]), 0, eps))
            return false;

        // note that calc-expec functions which assume
        // hermiticity will only consult the real component,
        // so its magnitude alone should determine divergence.
        // But alas an elem = eps*i will pass the above validation
        // yet also be validly Hermitian (imag <= eps), and cause
        // real(elem)=0 to be accepted within the matrix. This
        // will cause a divergence or divison-by-zero error when
        // the matrix is raised to a negative exponent; as this
        // function was supposed to detect and prevent! Fixing
        // this thoroughly would necessitate creating another
        // matrix field, separating when .absIsApproxNonZero and
        // .realIsApproxNonZero. But this is revolting and we 
        // simply accept the above strange scenario as a
        // non-validated dge-case.
    }

    return true;
}

// non-zeroness of fixed-size matrices is always computed afresh
bool util_isApproxNonZero(DiagMatr1 m, qreal eps) { return getWhetherNonZero(m.elems, m.numElems, eps); }
bool util_isApproxNonZero(DiagMatr2 m, qreal eps) { return getWhetherNonZero(m.elems, m.numElems, eps); }

// non-zeroness of heap matrices is cached
bool util_isApproxNonZero(DiagMatr matrix, qreal eps) {

    // compute and record whether matrix >= 0 if not already known
    if (*(matrix.isApproxNonZero) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(matrix.isApproxNonZero) = getWhetherNonZero(matrix.cpuElems, matrix.numElems, eps);

    // eps may be ignored
    return *(matrix.isApproxNonZero);
}
bool util_isApproxNonZero(FullStateDiagMatr matrix, qreal eps) {

    // compute and record whether matrix >= 0 if not already known
    if (*(matrix.isApproxNonZero) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(matrix.isApproxNonZero) = getWhetherNonZero(matrix.cpuElems, matrix.numElemsPerNode, eps);

    // communication may be necessary
    if (matrix.isDistributed)
        *(matrix.isApproxNonZero) = comm_isTrueOnAllNodes(*(matrix.isApproxNonZero));

    return *(matrix.isApproxNonZero);
}



/*
 * EXPONENTIABLE MATRIX IS NON-NEGATIVE
 */

bool getWhetherRealsAreStrictlyNonNegative(qcomp* diags, qindex dim) {

    /// @todo
    /// consider multithreading or GPU-accelerating this
    /// when caller is big and e.g. has GPU memory

    // it may seem like this function should be combined with
    // getWhetherRealsAreApproxNonZero() above, or indeed with
    // getHermiticity(), since all three are relevant to vali-
    // dation of diagonl matrices which can be exponentiated.
    // Alas, such as design is complicated by this particular
    // property (non-negativeness) being independent of the
    // validation epsilon. For example, sometimes it needs to
    // be computed even when numerical validation is disabled,
    // and does not need recomputing with the epsilon changes.
    
    // check every real component is strictly >= 0
    for (qindex i=0; i<dim; i++)
        if (std::real(diags[i]) < 0)
            return false;

    return true;
}

// non-negativity of fixed-size matrices is always computed afresh
bool util_isStrictlyNonNegative(DiagMatr1 m) { return getWhetherRealsAreStrictlyNonNegative(m.elems, m.numElems); }
bool util_isStrictlyNonNegative(DiagMatr2 m) { return getWhetherRealsAreStrictlyNonNegative(m.elems, m.numElems); }

bool util_isStrictlyNonNegative(DiagMatr m) {

    // compute and record whether m >= 0 if not already known
    if (*(m.isStrictlyNonNegative) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isStrictlyNonNegative) = getWhetherRealsAreStrictlyNonNegative(m.cpuElems, m.numElems);

    // eps may be ignored
    return *(m.isStrictlyNonNegative);
}

bool util_isStrictlyNonNegative(FullStateDiagMatr m) {

    // compute and record whether m >= 0 if not already known
    if (*(m.isStrictlyNonNegative) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(m.isStrictlyNonNegative) = getWhetherRealsAreStrictlyNonNegative(m.cpuElems, m.numElemsPerNode);

    // communication may be necessary
    if (m.isDistributed)
        *(m.isStrictlyNonNegative) = comm_isTrueOnAllNodes(*(m.isStrictlyNonNegative));

    return *(m.isStrictlyNonNegative);
}



/*
 * PAULI STR SUM HERMITICITY
 */

bool util_isHermitian(PauliStrSum sum, qreal eps) {

    // check whether all coefficients are real (just like a diagonal matrix)
    if (*(sum.isApproxHermitian) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(sum.isApproxHermitian) = getHermiticity(sum.coeffs, sum.numTerms, eps);

    // eps may have been ignored
    return *(sum.isApproxHermitian);
}



/*
 * KRAUS MAPS
 */

bool util_isCPTP(KrausMap map, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    // use pre-computed CPTP if it exists
    if (*(map.isApproxCPTP) != validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        return *(map.isApproxCPTP);

    /// @todo
    /// if KrausMap is GPU-accelerated, we should maybe
    /// instead perform this calculation using the GPU.
    /// otherwise, if matrix is large, we should potentially
    /// use a multithreaded routine

    *(map.isApproxCPTP) = 1;

    // check whether each element satisfies Identity = sum dagger(m)*m
    for (qindex r=0; r<map.numRows; r++) {
        for (qindex c=0; c<map.numRows; c++) {

            // calculate (r,c)-th element of sum dagger(m)*m
            qcomp elem = 0;
            for (int n=0; n<map.numMatrices; n++)
                for (qindex k=0; k<map.numRows; k++)
                    elem += std::conj(map.matrices[n][k][r]) * map.matrices[n][k][c];

            // fail if too distant from Identity element...
            qreal distSquared = std::norm(elem - (r==c));
            if (distSquared > eps) {

                // by recording the result and returning immediately
                *(map.isApproxCPTP) = 0;
                return *(map.isApproxCPTP);
            }
        }
    }

    // always true by this point
    return *(map.isApproxCPTP);
}

// T can be qcomp*** or vector<vector<vector<qcomp>>>
template <typename T> 
void setSuperoperator(qcomp** superop, T matrices, int numMatrices, qindex logMatrixDim) {

    /// @todo
    /// we initialise the superoperator completely serially, under the assumption that the
    /// superoperator will be small in size and initialised infrequently. Still, it would
    /// be better to provide backend initialisation functions (OpenMP and CUDA accelerated),
    /// called when the superoperator size is above some threshold!

    qindex matrixDim = powerOf2(logMatrixDim);
    qindex superopDim = matrixDim * matrixDim;

    // clear superoperator
    for (qindex r=0; r<superopDim; r++)
        for (qindex c=0; c<superopDim; c++)
            superop[r][c] = 0;

    // add each matrix's contribution to the superoperator
    for (int n=0; n<numMatrices; n++) {
        auto matrix = matrices[n];
        
        // superop += conj(matrix) (tensor) matrix
        for (qindex i=0; i<matrixDim; i++)
            for (qindex j=0; j<matrixDim; j++)
                for (qindex k=0; k<matrixDim; k++)
                    for (qindex l=0; l<matrixDim; l++) {
                        qindex r = i*matrixDim + k;
                        qindex c = j*matrixDim + l;
                        superop[r][c] += std::conj(matrix[i][j]) * matrix[k][l];
                    }
    }
}
void util_setSuperoperator(qcomp** superop, vector<vector<vector<qcomp>>> matrices, int numMatrices, int numQubits) {
    setSuperoperator(superop, matrices, numMatrices, numQubits);
}
void util_setSuperoperator(qcomp** superop, qcomp*** matrices, int numMatrices, int numQubits) {
    setSuperoperator(superop, matrices, numMatrices, numQubits);
}



/*
 * STRUCT PROPERTY CACHING
 */

std::list<int*> globalStructFieldPtrs(0);

void util_setFlagToUnknown(int* ptr) {

    *ptr = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
}

int* util_allocEpsilonSensitiveHeapFlag() {

    int* ptr = cpu_allocHeapFlag(); // may be nullptr

    // if failed to alloc, do not add to global list;
    // caller will handle validation/error messaging
    if (!mem_isAllocated(ptr))
        return ptr;

    // caller should set ptr to the default "unknown"
    // value, but we do so here too just to be safe
    util_setFlagToUnknown(ptr);

    // store the pointer so that we can reset the
    // value to "unknown" when epsilon is changed
    globalStructFieldPtrs.push_back(ptr);
    return ptr;
}

void util_deallocEpsilonSensitiveHeapFlag(int* ptr) {

    // nothing to do if the heap flag wasn't allocated;
    // it was never added to the list nor needs freeing
    if (!mem_isAllocated(ptr))
        return;

    globalStructFieldPtrs.remove(ptr);
    cpu_deallocHeapFlag(ptr);
}

void util_setEpsilonSensitiveHeapFlagsToUnknown() {

    for (auto ptr : globalStructFieldPtrs)
        util_setFlagToUnknown(ptr);
}



/*
 * DISTRIBUTED ELEMENTS INDEXING
 */

bool util_areAnyVectorElemsWithinNode(int rank, qindex numElemsPerNode, qindex elemStartInd, qindex numInds) {

    qindex elemEndIndExcl = elemStartInd + numInds;
    qindex nodeStartInd = numElemsPerNode * rank;
    qindex nodeEndIndExcl = nodeStartInd + numElemsPerNode;

    // 'no' if all targeted elems occur after this node
    if (elemStartInd >= nodeEndIndExcl)
        return false;

    // 'no' if all targeted elems occur before this node
    if (elemEndIndExcl <= nodeStartInd)
        return false;

    // otherwise yes; this node MUST contain one or more targeted elems
    return true;
}

util_VectorIndexRange util_getLocalIndRangeOfVectorElemsWithinNode(int rank, qindex numElemsPerNode, qindex elemStartInd, qindex numInds) {

    if (!util_areAnyVectorElemsWithinNode(rank, numElemsPerNode, elemStartInd, numInds))
        error_nodeUnexpectedlyContainedNoElems();

    // global indices of the user's targeted elements
    qindex elemEndInd = elemStartInd + numInds;

    // global indices of all elements contained in the node
    qindex nodeStartInd = numElemsPerNode * rank;
    qindex nodeEndInd   = nodeStartInd + numElemsPerNode;

    // global indices of user's targeted elements which are contained within node
    qindex globalRangeStartInd = std::max(elemStartInd, nodeStartInd);
    qindex globalRangeEndInd   = std::min(elemEndInd,   nodeEndInd);

    util_VectorIndexRange out;
    out.numElems = globalRangeEndInd - globalRangeStartInd;           // number of local elems in range
    out.localDistribStartInd = globalRangeStartInd % numElemsPerNode; // local inds of user's targeted elems to overwrite
    out.localDuplicStartInd  = globalRangeStartInd - elemStartInd;    // local inds of user's passed elems that correspond to above
    return out;
}



/*
 * GATE PARAMETERS
 */

qreal util_getPhaseFromGateAngle(qreal angle) {

    return - angle / 2;
}



/*
 * DECOHERENCE FACTORS
 */

qreal util_getOneQubitDephasingFactor(qreal prob) {

    return 1 - (2 * prob);
}

qreal util_getTwoQubitDephasingTerm(qreal prob) {

    return - 4 * prob / 3;
}

util_Scalars util_getOneQubitDepolarisingFactors(qreal prob) {

    util_Scalars out;

    // effected where braQubit == ketQubit
    out.c1 = 1 - (2 * prob / 3); // AA
    out.c2 = 2 * prob / 3;       // BB

    // effected where braQubit != ketQubit
    out.c3 = 1 - (4 * prob / 3); // AB

    // not used
    out.c4 = 0;

    return out;
}

util_Scalars util_getTwoQubitDepolarisingFactors(qreal prob) {

    util_Scalars out;

    out.c1 = 1 - (4 * prob / 5);
    out.c2 = 4 * prob / 15;
    out.c3 = - (16 * prob / 15);

    out.c4 = 0; // not used

    return out;
}

util_Scalars util_getOneQubitPauliChannelFactors(qreal pI, qreal pX, qreal pY, qreal pZ) {

    util_Scalars out;

    // effected where braQubit == ketQubit
    out.c1 = pI + pZ; // AA
    out.c2 = pX + pY; // BB

    // effected where braQubit != ketQubit
    out.c3 = pI - pZ; // AB
    out.c4 = pX - pY; // BA

    return out;
}

util_Scalars util_getOneQubitDampingFactors(qreal prob) {

    util_Scalars out;

    // we assume 0 < prob < 1 (true even of the inverse channel), so c1 is always real
    out.c1 = std::sqrt(1 - prob);
    out.c2 = 1 - prob;

    // not used
    out.c3 = 0;
    out.c4 = 0;

    return out;
}

qreal util_getMaxProbOfOneQubitDephasing() {
    return 1/2.;
}

qreal util_getMaxProbOfTwoQubitDephasing() {
    return 3/4.;
}

qreal util_getMaxProbOfOneQubitDepolarising() {
    return 3/4.;
}

qreal util_getMaxProbOfTwoQubitDepolarising() {
    return 15/16.;
}

// no equivalent function for oneQubitDamping, which
// can accept any valid probability, because we permit
// it to exceed maximal-mixing and induce purity



/*
 * TEMPORARY MEMORY ALLOCATION
 */

template <typename T>
void tryAllocVector(vector<T> &vec, qindex size, std::function<void()> errFunc) {

    // this function resizes the vector, not only reserving it,
    // such that vec.size() will subsequently return 'size'

    if (size == 0)
        return;

    try {
        vec.resize(size);

    } catch (std::bad_alloc &e) { 
        errFunc();
    } catch (std::length_error &e) {
        errFunc();
    }
}

void util_tryAllocVector(vector<qreal>    &vec, qindex size, std::function<void()> errFunc) { tryAllocVector(vec, size, errFunc); }
void util_tryAllocVector(vector<qcomp>    &vec, qindex size, std::function<void()> errFunc) { tryAllocVector(vec, size, errFunc); }
void util_tryAllocVector(vector<qcomp*>   &vec, qindex size, std::function<void()> errFunc) { tryAllocVector(vec, size, errFunc); }
void util_tryAllocVector(vector<unsigned> &vec, qindex size, std::function<void()> errFunc) { tryAllocVector(vec, size, errFunc); }

// cuQuantum needs a vector<double> overload, which we additionally define when qreal!=double. Gross!
#if FLOAT_PRECISION != 2
    void util_tryAllocVector(vector<double> &vec, qindex size, std::function<void()> errFunc) { tryAllocVector(vec, size, errFunc); }
#endif


void util_tryAllocMatrix(vector<vector<qcomp>> &matr, qindex numRows, qindex numCols, std::function<void()> errFunc) {

    // this function resizes the matrix, not only reserving it,
    // such that matr.size() will subsequently return 'numRows'
    // (unless numCols=0), and matr[0].size() will return 'numCols'

    if (numRows == 0 || numCols == 0)
        return;

    try {
        // alloc span of rows
        matr.resize(numRows);

        // alloc each row (serially enumerate since we expected matrix is small/tractable)
        for (qindex r=0; r<numRows; r++)
            matr[r].resize(numCols);

    } catch (std::bad_alloc &e) { 
        errFunc();
    } catch (std::length_error &e) {
        errFunc();
    }
}
