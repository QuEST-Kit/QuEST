/** @file
 * Miscellaneous utility functions needed internally.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <algorithm>
#include <complex>
#include <vector>

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

bool util_isBraQubitInSuffix(int ketQubit, Qureg qureg) {
    if (!qureg.isDensityMatrix)
        error_utilsIsBraQubitInSuffixGivenNonDensMatr();

    return ketQubit < qureg.logNumColsPerNode;
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

int util_getRankWithQubitFlipped(int ketQubit, Qureg qureg) {

    int rankInd = util_getPrefixInd(ketQubit, qureg);
    int rankFlip = flipBit(qureg.rank, rankInd);
    return rankFlip;
}

int util_getRankWithBraQubitFlipped(int ketQubit, Qureg qureg) {

    int rankInd = util_getPrefixBraInd(ketQubit, qureg);
    int rankFlip = flipBit(qureg.rank, rankInd);
    return rankFlip;
}

vector<int> util_getBraQubits(vector<int> ketQubits, Qureg qureg) {

    vector<int> braInds(0);
    braInds.reserve(ketQubits.size());

    for (int qubit : ketQubits)
        braInds.push_back(util_getBraQubit(qubit, qureg));

    return braInds;
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

    return vector<int> (qubits, qubits + numQubits);
}



/*
 * MATRIX CONJUGATION
 */

// type T can be qcomp** or qcomp*[]
template <typename T>
void setDenseElemsConj(T elems, qindex dim) {
    for (qindex i=0; i<dim; i++)
        for (qindex j=0; j<dim; j++)
           elems[i][j] = conj(elems[i][j]);
}

// diagonals don't need templating because arrays decay to pointers, yay!
void setDiagElemsConj(qcomp* elems, qindex dim) {
    for (qindex i=0; i<dim; i++)
        elems[i] = conj(elems[i]);
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
bool isUnitary(T elems, qindex dim) {

    qreal epsSq = VALIDATION_EPSILON * VALIDATION_EPSILON;

    // check m * dagger(m) == identity
    for (qindex r=0; r<dim; r++) {
        for (qindex c=0; c<dim; c++) {

            // compute m[r,...] * dagger(m)[...,c]
            qcomp elem = 0;
            for (qindex i=0; i<dim; i++)
                elem += elems[r][i] * conj(elems[c][i]);

            // check if further than epsilon from identity[r,c]
            qcomp dif = elem - qcomp(r == c, 0);
            qreal dist = real(dif)*real(dif) + imag(dif)*imag(dif);
            if (dist > epsSq)
                return false;
        }
    }

    return true;
}

// diagonal version doesn't need templating because array decays to pointer, yay!
bool isUnitary(qcomp* diags, qindex dim) {

    // check every element has unit magnitude
    for (qindex i=0; i<dim; i++) {
        qreal mag = std::abs(diags[i]);
        qreal dif = std::abs(1 - mag);

        if (dif > VALIDATION_EPSILON)
            return false;
    }

    return true;
}

bool util_isUnitary(CompMatr1 matrix) {
    return isUnitary(matrix.elems, matrix.numRows);
}
bool util_isUnitary(CompMatr2 matrix) {
    return isUnitary(matrix.elems, matrix.numRows);
}
bool util_isUnitary(CompMatr matrix) {
    return isUnitary(matrix.cpuElems, matrix.numRows);
}

bool util_isUnitary(DiagMatr1 matrix) {
    return isUnitary(matrix.elems, matrix.numElems);
}
bool util_isUnitary(DiagMatr2 matrix) {
    return isUnitary(matrix.elems, matrix.numElems);
}
bool util_isUnitary(DiagMatr matrix) {
    return isUnitary(matrix.cpuElems, matrix.numElems);
}

bool util_isUnitary(FullStateDiagMatr matrix) {

    // we must check all node's sub-diagonals satisfy unitarity
    bool res = isUnitary(matrix.cpuElems, matrix.numElems);
    if (comm_isInit())
        res = comm_isTrueOnAllNodes(res);

    return res;
}



/*
 * KRAUS MAPS
 */

bool util_isCPTP(KrausMap map) {

    qreal epsSquared = VALIDATION_EPSILON * VALIDATION_EPSILON;

    // each whether each element satisfies Identity = sum dagger(m)*m
    for (qindex r=0; r<map.numRows; r++) {
        for (qindex c=0; c<map.numRows; c++) {

            // calculate (r,c)-th element of sum dagger(m)*m
            qcomp elem = 0;
            for (int n=0; n<map.numMatrices; n++)
                for (qindex k=0; k<map.numRows; k++)
                    elem += conj(map.matrices[n][k][r]) * map.matrices[n][k][c];

            // fail if too distant from Identity element
            qreal distSquared = norm(elem - (r==c));
            if (distSquared > epsSquared)   
                return false;
        }
    }
    
    return true;
}

// T can be qcomp*** or vector<vector<vector<qcomp>>>
template <typename T> 
void setSuperoperator(qcomp** superop, T matrices, int numMatrices, qindex logMatrixDim) {

    // TODO:
    // we initialise the superoperator completely serially, under the assumption that the
    // superoperator will be small in size and initialised infrequently. Still, it would
    // be better to provide backend initialisation functions (OpenMP and CUDA accelerated),
    // called when the superoperator size is above some threshold!

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
                        superop[r][c] += conj(matrix[i][j]) * matrix[k][l];
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

bool util_areAnyElemsWithinThisNode(int numElemsPerNode, qindex elemStartInd, qindex numInds) {

    qindex elemEndIndExcl = elemStartInd + numInds;
    qindex nodeStartInd = comm_getRank() * numElemsPerNode;
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

util_IndexRange util_getLocalIndRangeOfElemsWithinThisNode(int numElemsPerNode, qindex elemStartInd, qindex numInds) {

    if (!util_areAnyElemsWithinThisNode(numElemsPerNode, elemStartInd, numInds))
        error_nodeUnexpectedlyContainedNoElems();

    // global indices of the user's targeted elements
    qindex elemEndInd = elemStartInd + numInds;

    // global indices of all elements contained in the node
    qindex nodeStartInd = numElemsPerNode * comm_getRank();
    qindex nodeEndInd   = nodeStartInd + numElemsPerNode;

    // global indices of user's targeted elements which are contained within node
    qindex globalRangeStartInd = std::max(elemStartInd, nodeStartInd);
    qindex globalRangeEndInd   = std::min(elemEndInd,   nodeEndInd);
    qindex numLocalElems       = globalRangeEndInd - globalRangeStartInd;

    // local indices of user's targeted elements to overwrite
    qindex localRangeStartInd = globalRangeStartInd % numElemsPerNode;
    
    // local indices of user's passed elements that correspond to above
    qindex localOffsetInd = globalRangeStartInd - elemStartInd;
    
    return (util_IndexRange) {
        .localDistribStartInd = localRangeStartInd,
        .localDuplicStartInd = localOffsetInd,
        .numElems = numLocalElems
    };
}



/*
 * OPERATOR PARAMETERS
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
    qreal c1 = sqrt(1 - prob);
    qreal c2 = 1 - prob;

    return {.c1=c1, .c2=c2, .c3=0, .c4=0}; //c3 and c4 ignored
}
