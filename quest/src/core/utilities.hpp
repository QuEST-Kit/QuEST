/** @file
 * Miscellaneous utility functions needed internally.
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"
#include "quest/include/environment.h"

#include <type_traits>
#include <string>
#include <vector>
#include <array>

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

std::array<vector<int>,2> util_getPrefixAndSuffixQubits(vector<int> qubits, Qureg qureg);

int util_getRankBitOfQubit(int ketQubit, Qureg qureg);
int util_getRankBitOfBraQubit(int ketQubit, Qureg qureg);

int util_getRankWithQubitFlipped(int ketQubit, Qureg qureg);
int util_getRankWithBraQubitFlipped(int ketQubit, Qureg qureg);
int util_getRankWithQubitsFlipped(vector<int> prefixQubits, Qureg qureg);

vector<int> util_getBraQubits(vector<int> ketQubits, Qureg qureg);

vector<int> util_getVector(int* qubits, int numQubits);

vector<int> util_getConcatenated(vector<int> list1, vector<int> list2);

vector<int> util_getSorted(vector<int> list);
vector<int> util_getSorted(vector<int> ctrls, vector<int> targs);

qindex util_getBitMask(vector<int> qubits);
qindex util_getBitMask(vector<int> qubits, vector<int> states);
qindex util_getBitMask(vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, vector<int> targStates);



/*
 * INDEX ALGEBRA
 */

qindex util_getGlobalIndexOfFirstLocalAmp(Qureg qureg);

qindex util_getLocalIndexOfGlobalIndex(Qureg qureg, qindex globalInd);

qindex util_getLocalIndexOfFirstDiagonalAmp(Qureg qureg);

qindex util_getNumLocalDiagonalAmpsWithBits(Qureg qureg, vector<int> qubits, vector<int> outcomes);

qindex util_getGlobalFlatIndex(Qureg qureg, qindex globalRow, qindex globalCol);

int util_getRankContainingIndex(Qureg qureg, qindex globalInd);
int util_getRankContainingColumn(Qureg qureg, qindex globalCol);
int util_getRankContainingIndex(FullStateDiagMatr matr, qindex globalInd);



/*
 * COMPLEX ALGEBRA
 */

qcomp util_getPowerOfI(size_t exponent);



/*
 * STRUCT TYPING
 *
 * defined here in the header since templated, and which use compile-time inspection.
 */

template <class T> constexpr bool util_isQuregType() { return is_same_v<T, Qureg>; }
template <class T> constexpr bool util_isCompMatr1() { return is_same_v<T, CompMatr1>; }
template <class T> constexpr bool util_isCompMatr2() { return is_same_v<T, CompMatr2>; }
template <class T> constexpr bool util_isCompMatr () { return is_same_v<T, CompMatr >; }
template <class T> constexpr bool util_isDiagMatr1() { return is_same_v<T, DiagMatr1>; }
template <class T> constexpr bool util_isDiagMatr2() { return is_same_v<T, DiagMatr2>; }
template <class T> constexpr bool util_isDiagMatr () { return is_same_v<T, DiagMatr >; }
template <class T> constexpr bool util_isFullStateDiagMatr () { return is_same_v<T, FullStateDiagMatr >; }

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

    // this line is reached if the type is not a matrix
    return false;
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
constexpr bool util_isDistributableType() {

    return (is_same_v<T, FullStateDiagMatr> || is_same_v<T, Qureg>);
}

template<class T>
bool util_isDistributedMatrix(T matr) {

    if constexpr (util_isDistributableType<T>())
        return matr.isDistributed;

    return false;
}

template<class T>
bool util_isGpuAcceleratedMatrix(T matr) {

    if constexpr (util_isFullStateDiagMatr<T>())
        return matr.isGpuAccelerated;

    if constexpr (util_isHeapMatrixType<T>())
        return getQuESTEnv().isGpuAccelerated;

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

bool util_isUnitary(CompMatr1 matrix, qreal epsilon);
bool util_isUnitary(CompMatr2 matrix, qreal epsilon);
bool util_isUnitary(CompMatr  matrix, qreal epsilon);
bool util_isUnitary(DiagMatr1 matrix, qreal epsilon);
bool util_isUnitary(DiagMatr2 matrix, qreal epsilon);
bool util_isUnitary(DiagMatr  matrix, qreal epsilon);
bool util_isUnitary(FullStateDiagMatr matrix, qreal epsilon);



/*
 * MATRIX HERMITICITY
 */

bool util_isHermitian(CompMatr1 matrix, qreal epsilon);
bool util_isHermitian(CompMatr2 matrix, qreal epsilon);
bool util_isHermitian(CompMatr  matrix, qreal epsilon);
bool util_isHermitian(DiagMatr1 matrix, qreal epsilon);
bool util_isHermitian(DiagMatr2 matrix, qreal epsilon);
bool util_isHermitian(DiagMatr  matrix, qreal epsilon);
bool util_isHermitian(FullStateDiagMatr matrix, qreal epsilon);



/*
 * PAULI STR SUM HERMITICITY
 */

bool util_isHermitian(PauliStrSum sum, qreal epsilon);



/*
 * KRAUS MAPS AND SUPEROPERATORS
 */

bool util_isCPTP(KrausMap map, qreal epsilon);

void util_setSuperoperator(qcomp** superop, vector<vector<vector<qcomp>>> matrices, int numMatrices, int numQubits);
void util_setSuperoperator(qcomp** superop, qcomp*** matrices, int numMatrices, int numQubits);



/*
 * DISTRIBUTED ELEMENTS INDEXING
 */

struct util_VectorIndexRange {

    // the first local index of this node's amps which are in the queried distributed range
    qindex localDistribStartInd;

    // the corresponding local index of the non-distributed (i.e. duplicated on every node) data structure
    qindex localDuplicStartInd;

    // the number of this node's amps which are within the queried distributed range
    qindex numElems;
};

bool util_areAnyVectorElemsWithinNode(int rank, qindex numElemsPerNode, qindex startInd, qindex numInds);

util_VectorIndexRange util_getLocalIndRangeOfVectorElemsWithinNode(int rank, qindex numElemsPerNode, qindex elemStartInd, qindex numInds);



/*
 * GATE PARAMETERS
 */

qreal util_getPhaseFromGateAngle(qreal angle);



/*
 * DECOHERENCE FACTORS
 */

struct util_Scalars { qreal c1; qreal c2; qreal c3; qreal c4; };

qreal util_getOneQubitDephasingFactor(qreal prob);

qreal util_getTwoQubitDephasingTerm(qreal prob);

util_Scalars util_getOneQubitDepolarisingFactors(qreal prob);

util_Scalars util_getTwoQubitDepolarisingFactors(qreal prob);

util_Scalars util_getOneQubitDampingFactors(qreal prob);

util_Scalars util_getOneQubitPauliChannelFactors(qreal pI, qreal pX, qreal pY, qreal pZ);

qreal util_getMaxProbOfOneQubitDephasing();

qreal util_getMaxProbOfTwoQubitDephasing();

qreal util_getMaxProbOfOneQubitDepolarising();

qreal util_getMaxProbOfTwoQubitDepolarising();



/*
 * TEMPORARY MEMORY ALLOCATION
 */

void util_tryAllocVector(vector<qreal > &vec, qindex size, void (*errFunc)());
void util_tryAllocVector(vector<qcomp > &vec, qindex size, void (*errFunc)());
void util_tryAllocVector(vector<qcomp*> &vec, qindex size, void (*errFunc)());

// cuQuantum needs a vector<double> overload, which we additionally define when qreal!=double. Gross!
#if FLOAT_PRECISION != 2
void util_tryAllocVector(vector<double> &vec, qindex size, void (*errFunc)());
#endif

void util_tryAllocMatrix(vector<vector<qcomp>> &vec, qindex numRows, qindex numCols, void (*errFunc)());



#endif // UTILITIES_HPP