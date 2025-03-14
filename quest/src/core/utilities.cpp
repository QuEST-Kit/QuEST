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
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <algorithm>
#include <complex>
#include <vector>
#include <array>
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

    // when not distributed, all nodes (each believing themselves root) contain the index
    if (!qureg.isDistributed)
        return 0;

    // accepts flat density matrix index too
    return globalInd / qureg.numAmpsPerNode; // floors
}
int util_getRankContainingIndex(FullStateDiagMatr matr, qindex globalInd) {

    // when not distributed, all nodes (each believing themselves root) contain the index
    if (!matr.isDistributed)
        return 0;

    return globalInd / matr.numElemsPerNode; // floors
}

int util_getRankContainingColumn(Qureg qureg, qindex globalCol) {
    assert_utilsGivenDensMatr(qureg);

    // when not distributed, all nodes (each believing themselves root) contain the index
    if (!qureg.isDistributed)
        return 0;

    qindex numColsPerNode = powerOf2(qureg.logNumColsPerNode);
    return globalCol / numColsPerNode; // floors
}

qindex util_getNextPowerOf2(qindex number) {

    int nextExponent = static_cast<int>(std::ceil(std::log2(number)));
    return powerOf2(nextExponent);
}



/*
 * COMPLEX ALGEBRA
 */

qcomp util_getPowerOfI(size_t exponent) {

    // seems silly, but at least it's precision agnostic!
    qcomp values[] = {1, 1_i, -1, -1_i};
    return values[exponent % 4];
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

// type T can be qcomp** or qcomp*[]
template <typename T>
bool isUnitary(T elems, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    // check m * dagger(m) == identity
    for (qindex r=0; r<dim; r++) {
        for (qindex c=0; c<dim; c++) {

            // compute m[r,...] * dagger(m)[...,c]
            qcomp elem = 0;
            for (qindex i=0; i<dim; i++)
                elem += elems[r][i] * std::conj(elems[c][i]);

            // check if further than epsilon from identity[r,c]
            qcomp dif = elem - qcomp(r == c, 0);
            qreal distSq = std::norm(dif);
            if (distSq > eps)
                return false;
        }
    }

    return true;
}

// diagonal version doesn't need templating because array decays to pointer, yay!
bool isUnitary(qcomp* diags, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    // check every element has unit magnitude
    for (qindex i=0; i<dim; i++) {
        qreal mag = std::abs(diags[i]);
        qreal dif = std::abs(1 - mag);

        if (dif > eps)
            return false;
    }

    return true;
}

bool util_isUnitary(CompMatr1 matrix, qreal eps) {
    return isUnitary(matrix.elems, matrix.numRows, eps);
}
bool util_isUnitary(CompMatr2 matrix, qreal eps) {
    return isUnitary(matrix.elems, matrix.numRows, eps);
}
bool util_isUnitary(CompMatr matrix, qreal eps) {

    /// @todo
    /// if matrix is GPU-accelerated, we should maybe
    /// instead perform this calculation using the GPU.
    /// otherwise, if matrix is large, we should potentially
    /// use a multithreaded routine

    return isUnitary(matrix.cpuElems, matrix.numRows, eps);
}

bool util_isUnitary(DiagMatr1 matrix, qreal eps) {
    return isUnitary(matrix.elems, matrix.numElems, eps);
}
bool util_isUnitary(DiagMatr2 matrix, qreal eps) {
    return isUnitary(matrix.elems, matrix.numElems, eps);
}
bool util_isUnitary(DiagMatr matrix, qreal eps) {

    /// @todo
    /// if matrix is GPU-accelerated, we should maybe
    /// instead perform this calculation using the GPU.
    /// otherwise, if matrix is large, we should potentially
    /// use a multithreaded routine

    return isUnitary(matrix.cpuElems, matrix.numElems, eps);
}

bool util_isUnitary(FullStateDiagMatr matrix, qreal eps) {

    /// @todo
    /// we should definitely be using an accelerated routine
    /// here, e.g. GPU-acceleration or multithreading

    // we must check all node's sub-diagonals satisfy unitarity
    bool res = isUnitary(matrix.cpuElems, matrix.numElems, eps);
    if (matrix.isDistributed)
        res = comm_isTrueOnAllNodes(res);

    return res;
}



/*
 * MATRIX HERMITICITY
 */

// type T can be qcomp** or qcomp*[]
template <typename T>
bool isHermitian(T elems, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    // check adj(elems) == elems
    for (qindex r=0; r<dim; r++) {
        for (qindex c=0; c<r; c++) {

            qcomp dif = elems[r][c] - std::conj(elems[c][r]);
            qreal distSq = std::norm(dif);
            if (distSq > eps)
                return false;
        }
    }

    return true;
}

// diagonal version doesn't need templating because array decays to pointer, yay!
bool isHermitian(qcomp* diags, qindex dim, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    // check every element has a zero (or <eps) imaginary component
    for (qindex i=0; i<dim; i++)
        if (std::abs(std::imag(diags[i])) > eps)
            return false;

    return true;
}

bool util_isHermitian(CompMatr1 matrix, qreal eps) {
    return isHermitian(matrix.elems, matrix.numRows, eps);
}
bool util_isHermitian(CompMatr2 matrix, qreal eps) {
    return isHermitian(matrix.elems, matrix.numRows, eps);
}
bool util_isHermitian(CompMatr matrix, qreal eps) {

    /// @todo
    /// if matrix is GPU-accelerated, we should maybe
    /// instead perform this calculation using the GPU.
    /// otherwise, if matrix is large, we should potentially
    /// use a multithreaded routine

    return isHermitian(matrix.cpuElems, matrix.numRows, eps);
}

bool util_isHermitian(DiagMatr1 matrix, qreal eps) {
    return isHermitian(matrix.elems, matrix.numElems, eps);
}
bool util_isHermitian(DiagMatr2 matrix, qreal eps) {
    return isHermitian(matrix.elems, matrix.numElems, eps);
}
bool util_isHermitian(DiagMatr matrix, qreal eps) {

    /// @todo
    /// if matrix is GPU-accelerated, we should maybe
    /// instead perform this calculation using the GPU.
    /// otherwise, if matrix is large, we should potentially
    /// use a multithreaded routine

    return isHermitian(matrix.cpuElems, matrix.numElems, eps);
}

bool util_isHermitian(FullStateDiagMatr matrix, qreal eps) {

    /// @todo
    /// we should definitely be using an accelerated routine
    /// here, e.g. GPU-acceleration or multithreading

    // we must check all node's sub-diagonals satisfy unitarity
    bool res = isHermitian(matrix.cpuElems, matrix.numElems, eps);
    if (matrix.isDistributed)
        res = comm_isTrueOnAllNodes(res);

    return res;
}



/*
 * PAULI STR SUM HERMITICITY
 */

bool util_isHermitian(PauliStrSum sum, qreal eps) {

    // check whether all coefficients are real (just like a diagonal matrix)
    return isHermitian(sum.coeffs, sum.numTerms, eps);
}



/*
 * KRAUS MAPS
 */

bool util_isCPTP(KrausMap map, qreal eps) {
    assert_utilsGivenNonZeroEpsilon(eps);

    /// @todo
    /// if KrausMap is GPU-accelerated, we should maybe
    /// instead perform this calculation using the GPU.
    /// otherwise, if matrix is large, we should potentially
    /// use a multithreaded routine

    // each whether each element satisfies Identity = sum dagger(m)*m
    for (qindex r=0; r<map.numRows; r++) {
        for (qindex c=0; c<map.numRows; c++) {

            // calculate (r,c)-th element of sum dagger(m)*m
            qcomp elem = 0;
            for (int n=0; n<map.numMatrices; n++)
                for (qindex k=0; k<map.numRows; k++)
                    elem += std::conj(map.matrices[n][k][r]) * map.matrices[n][k][c];

            // fail if too distant from Identity element
            qreal distSquared = std::norm(elem - (r==c));
            if (distSquared > eps)   
                return false;
        }
    }
    
    return true;
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
    qindex numLocalElems       = globalRangeEndInd - globalRangeStartInd;

    // local indices of user's targeted elements to overwrite
    qindex localRangeStartInd = globalRangeStartInd % numElemsPerNode;
    
    // local indices of user's passed elements that correspond to above
    qindex localOffsetInd = globalRangeStartInd - elemStartInd;
    
    return {
        .localDistribStartInd = localRangeStartInd,
        .localDuplicStartInd = localOffsetInd,
        .numElems = numLocalElems
    };
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

    // effected where braQubit == ketQubit
    qreal facAA = 1 - (2 * prob / 3);
    qreal facBB = 2 * prob / 3;

    // effected where braQubit != ketQubit
    qreal facAB  = 1 - (4 * prob / 3);

    return {.c1=facAA, .c2=facBB, .c3=facAB, .c4=0}; // c4 ignored
}

util_Scalars util_getTwoQubitDepolarisingFactors(qreal prob) {

    return {
        .c1 = 1 - (4 * prob / 5), 
        .c2 = 4 * prob / 15, 
        .c3 = - (16 * prob / 15), 
        .c4 = 0 // ignored
    };
}

util_Scalars util_getOneQubitPauliChannelFactors(qreal pI, qreal pX, qreal pY, qreal pZ) {

    // effected where braQubit == ketQubit
    qreal facAA = pI + pZ;
    qreal facBB = pX + pY;

    // effected where braQubit != ketQubit
    qreal facAB = pI - pZ;
    qreal facBA = pX - pY;

    return {.c1=facAA, .c2=facBB, .c3=facAB, .c4=facBA};
}

util_Scalars util_getOneQubitDampingFactors(qreal prob) {

    // we assume 0 < prob < 1 (true even of the inverse channel), so c1 is always real
    qreal c1 = std::sqrt(1 - prob);
    qreal c2 = 1 - prob;

    return {.c1=c1, .c2=c2, .c3=0, .c4=0}; //c3 and c4 ignored
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
void tryAllocVector(vector<T> &vec, qindex size, void (*errFunc)()) {

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

void util_tryAllocVector(vector<qreal > &vec, qindex size, void (*errFunc)()) { tryAllocVector(vec, size, errFunc); }
void util_tryAllocVector(vector<qcomp > &vec, qindex size, void (*errFunc)()) { tryAllocVector(vec, size, errFunc); }
void util_tryAllocVector(vector<qcomp*> &vec, qindex size, void (*errFunc)()) { tryAllocVector(vec, size, errFunc); }

// cuQuantum needs a vector<double> overload, which we additionally define when qreal!=double. Gross!
#if FLOAT_PRECISION != 2
    void util_tryAllocVector(vector<double> &vec, qindex size, void (*errFunc)()) { tryAllocVector(vec, size, errFunc); }
#endif


void util_tryAllocMatrix(vector<vector<qcomp>> &matr, qindex numRows, qindex numCols, void (*errFunc)()) {

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
