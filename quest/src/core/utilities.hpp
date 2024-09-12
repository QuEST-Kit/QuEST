/** @file
 * Miscellaneous utility functions needed internally.
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"

#include <type_traits>
#include <string>
#include <vector>

using std::is_same_v;
using std::vector;



/*
 * QUBIT PROCESSING
 */

bool util_isQubitInSuffix(int qubit, Qureg qureg);
bool util_isBraQubitInSuffix(int ketQubit, Qureg qureg);

int util_getBraQubit(int ketQubit, Qureg qureg);

int util_getPrefixInd(int qubit, Qureg qureg);
int util_getPrefixBraInd(int ketQubit, Qureg qureg);

int util_getRankBitOfQubit(int ketQubit, Qureg qureg);
int util_getRankBitOfBraQubit(int ketQubit, Qureg qureg);

int util_getRankWithQubitFlipped(int ketQubit, Qureg qureg);
int util_getRankWithBraQubitFlipped(int ketQubit, Qureg qureg);

vector<int> util_getBraQubits(vector<int> ketQubits, Qureg qureg);

vector<int> util_getVector(int* qubits, int numQubits);

vector<int> util_getConcatenated(vector<int> list1, vector<int> list2);

vector<int> util_getSorted(vector<int> list);
vector<int> util_getSorted(vector<int> ctrls, vector<int> targs);

qindex util_getBitMask(vector<int> qubits);
qindex util_getBitMask(vector<int> qubits, vector<int> states);
qindex util_getBitMask(vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, vector<int> targStates);



/*
 * MATRIX TYPING
 *
 * defined here in the header since templated, and which use compile-time inspection.
 */

template<class T>
constexpr bool util_isDenseMatrixType() {

    // CompMatr, SuperOp and (in are sense) KrausMaps are "dense", storing all 2D elements
    if constexpr (
        is_same_v<T, CompMatr1> ||
        is_same_v<T, CompMatr2> ||
        is_same_v<T, CompMatr>  ||
        is_same_v<T, KrausMap>  ||
        is_same_v<T, SuperOp>
    )
        return true;

    // DiagMatr are "sparse", storing only the diagonals
    if constexpr (
        is_same_v<T, DiagMatr1> ||
        is_same_v<T, DiagMatr2> ||
        is_same_v<T, DiagMatr>  ||
        is_same_v<T, FullStateDiagMatr>
    )
        return false;

    // this line is unreachable but throwing errors in a template expansion is ludicrous;
    // above type checks are explicit in case we add more matrix types later
    return false;
}

template<class T>
constexpr bool util_isDiagonalMatrixType() {

    return !util_isDenseMatrixType<T>();
}

template<class T>
constexpr bool util_isFixedSizeMatrixType() {

    return (
        is_same_v<T, CompMatr1> ||
        is_same_v<T, CompMatr2> ||
        is_same_v<T, DiagMatr1> ||
        is_same_v<T, DiagMatr2>
    );
}

template<class T>
constexpr bool util_isHeapMatrixType() {

    // all non-fixed size matrices are stored in the heap (never the stack)
    return ! util_isFixedSizeMatrixType<T>();
}

template<class T>
constexpr bool util_isDistributableMatrixType() {

    return (is_same_v<T, FullStateDiagMatr>);
}

template<class T>
bool util_isDistributedMatrix(T matr) {

    if constexpr (util_isDistributableMatrixType<T>())
        return matr.isDistributed;

    return false;
}

template<class T>
std::string util_getMatrixTypeName() {
    
    if constexpr (is_same_v<T, CompMatr1>) return "CompMatr1";
    if constexpr (is_same_v<T, CompMatr2>) return "CompMatr2";
    if constexpr (is_same_v<T, CompMatr >) return "CompMatr" ;
    if constexpr (is_same_v<T, DiagMatr1>) return "DiagMatr1";
    if constexpr (is_same_v<T, DiagMatr2>) return "DiagMatr2";
    if constexpr (is_same_v<T, DiagMatr >) return "DiagMatr" ;
    if constexpr (is_same_v<T, FullStateDiagMatr>)
        return "FullStateDiagMatr";

    // these types do not need to have this function called, but
    // we include them for completeness
    if constexpr (is_same_v<T, KrausMap >) return "KrausMap" ;
    if constexpr (is_same_v<T, SuperOp >) return "SuperOp" ;

    // no need to create a new error for this situation
    return "UnrecognisedMatrix";
}

template<class T>
qindex util_getMatrixDim(T matr) {
    
    if constexpr (util_isDenseMatrixType<T>())
        return matr.numRows;
    else
        return matr.numElems;
}

// T can be CompMatr, DiagMatr, FullStateDiagMatr, KrausMap, SuperOp (i.e. heap-based non-Qureg structures)
template<class T>
qcomp util_getFirstLocalElem(T obj) {

    // Kraus map elements never reach the GPU; their superoperator elements do, so query those
    if constexpr (is_same_v<T, KrausMap>)
        return util_getFirstLocalElem(obj.superop);

    // otherwise the elems field dimension depends on the matrix type
    else if constexpr (util_isDenseMatrixType<T>())
        return obj.cpuElems[0][0];
    else
        return obj.cpuElems[0];
}

// T can be CompMatr, DiagMatr, FullStateDiagMatr, SuperOp (but NOT KrausMap)
template<class T>
qcomp* util_getGpuMemPtr(T matr) {

    // 2D CUDA structures are always stored as 1D
    if constexpr (util_isDenseMatrixType<T>())
        return matr.gpuElemsFlat;
    else
        return matr.gpuElems;
}



/*
 * MATRIX CONJUGATION
 */

CompMatr1 util_getConj(CompMatr1 matrix);
CompMatr2 util_getConj(CompMatr2 matrix);
DiagMatr1 util_getConj(DiagMatr1 matrix);
DiagMatr2 util_getConj(DiagMatr2 matrix);

void util_setConj(CompMatr matrix);
void util_setConj(DiagMatr matrix);



/*
 * MATRIX UNITARITY
 */

bool util_isUnitary(CompMatr1 matrix);
bool util_isUnitary(CompMatr2 matrix);
bool util_isUnitary(CompMatr matrix);
bool util_isUnitary(DiagMatr1 matrix);
bool util_isUnitary(DiagMatr2 matrix);
bool util_isUnitary(DiagMatr matrix);
bool util_isUnitary(FullStateDiagMatr matrix);



/*
 * KRAUS MAPS AND SUPEROPERATORS
 */

bool util_isCPTP(KrausMap map);

void util_setSuperoperator(qcomp** superop, vector<vector<vector<qcomp>>> matrices, int numMatrices, int numQubits);
void util_setSuperoperator(qcomp** superop, qcomp*** matrices, int numMatrices, int numQubits);



/*
 * DISTRIBUTED ELEMENTS INDEXING
 */

typedef struct {

    // the starting local index among the node's distributed elements
    qindex localDistribStartInd;

    // the corresponding local index of in the duplicated (non-distributed) elements 
    qindex localDuplicStartInd;

    // the number of elements in the range, all of which are in this node's distributed elements
    qindex numElems;

} util_IndexRange;

bool util_areAnyElemsWithinThisNode(int numElemsPerNode, qindex startInd, qindex numInds);

util_IndexRange util_getLocalIndRangeOfElemsWithinThisNode(int numElemsPerNode, qindex elemStartInd, qindex numInds);



/*
 * OPERATOR PARAMETERS
 */

typedef struct { qreal c1; qreal c2; qreal c3; qreal c4; } util_Scalars;

qreal util_getOneQubitDephasingFactor(qreal prob);

qreal util_getTwoQubitDephasingTerm(qreal prob);

util_Scalars util_getOneQubitDepolarisingFactors(qreal prob);

util_Scalars util_getTwoQubitDepolarisingFactors(qreal prob);

util_Scalars util_getOneQubitDampingFactors(qreal prob);

util_Scalars util_getOneQubitPauliChannelFactors(qreal pI, qreal pX, qreal pY, qreal pZ);



#endif // UTILITIES_HPP