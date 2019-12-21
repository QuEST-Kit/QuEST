/** @file
 * Simple implementations of matrix operations used by QuEST_unit_tests
 *
 * @author Tyson Jones
 */
 
 
 // @TODO: ADD THESE DOCS to doxygen with @group testing

#ifndef QUEST_MATRICES_H
#define QUEST_MATRICES_H

#include "QuEST.h"
#include "QuEST_complex.h"
#include "catch.hpp"
#include <vector>

typedef std::vector<std::vector<qcomp>> QMatrix;
typedef std::vector<qcomp> QVector;

/* modifying data-types */
QVector operator + (const QVector& v1, const QVector& v2);
QVector operator - (const QVector& v1, const QVector& v2);
QVector operator * (const qcomp& a, const QVector& v);
QVector operator * (const QVector& v, const qcomp& a);
QVector operator / (const QVector& v, const qcomp& a);
qcomp operator * (const QVector &v1, const QVector& v2);
void operator += (QVector& v1, const QVector& v2);
void operator -= (QVector& v1, const QVector& v2);
void operator *= (QVector& v1, const qcomp& a);
void operator /= (QVector& v1, const qcomp& a);
QMatrix operator + (const QMatrix& m1, const QMatrix& m2);
QMatrix operator - (const QMatrix& m1, const QMatrix& m2);
QMatrix operator * (const qcomp& a, const QMatrix& m);
QMatrix operator * (const QMatrix& m, const qcomp& a);
QMatrix operator / (const QMatrix& m, const qcomp& a);
QMatrix operator * (const QMatrix& m1, const QMatrix& m2);
void operator += (QMatrix& m1, const QMatrix& m2);
void operator -= (QMatrix& m1, const QMatrix& m2);
void operator *= (QMatrix& m1, const qreal& a);
void operator /= (QMatrix& m1, const qreal& a);
void operator *= (QMatrix& m1, const QMatrix& m2);
QVector operator * (const QMatrix& m, const QVector& v);

/* converting data-types */
QVector toQVector(Qureg qureg);
QMatrix toQMatrix(Complex alpha, Complex beta);
QMatrix toQMatrix(Qureg qureg);
ComplexMatrix2 toComplexMatrix2(QMatrix qm);
ComplexMatrix4 toComplexMatrix4(QMatrix qm);
void toComplexMatrixN(QMatrix qm, ComplexMatrixN cm);

/* building unitary matrices */
QMatrix getZeroMatrix(size_t dim);
QMatrix getExponentialDiagonalMatrix(QMatrix a);
QMatrix getExponentialPauliMatrix(qreal angle, QMatrix a);
QMatrix getKroneckerProduct(QMatrix a, QMatrix b);
QMatrix getFullOperatorMatrix(int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits); 

/* generating random inputs */
int getRandomInt(int min, int max); //  exclusive max
qreal getRandomReal(qreal min, qreal max);
QMatrix getRandomUnitary(int numQb);
std::vector<QMatrix> getRandomKrausMap(int numQb, int numOps);

/* applying operations to reference state-vectors and density matrices. 
 * These can be non-unitary, but will still be applied as U rho U^dagger to density 
 * matrices (useful for Kraus operators).
 */
void applyReferenceOp(QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);
void applyReferenceOp(QMatrix &state, int* ctrls, int numCtrls, int targ1, int targ2, QMatrix op);
void applyReferenceOp(QMatrix &state, int* ctrls, int numCtrls, int target, QMatrix op);
void applyReferenceOp(QMatrix &state, int *targs, int numTargs, QMatrix op);
void applyReferenceOp(QMatrix &state, int ctrl, int targ, QMatrix op);
void applyReferenceOp(QMatrix &state, int ctrl, int* targs, int numTargs, QMatrix op);
void applyReferenceOp(QMatrix &state, int ctrl, int targ1, int targ2, QMatrix op);
void applyReferenceOp(QMatrix &state, int targ, QMatrix op);
void applyReferenceOp(QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);
void applyReferenceOp(QVector &state, int* ctrls, int numCtrls, int targ1, int targ2, QMatrix op);
void applyReferenceOp(QVector &state, int* ctrls, int numCtrls, int target, QMatrix op);
void applyReferenceOp(QVector &state, int *targs, int numTargs, QMatrix op);
void applyReferenceOp(QVector &state, int ctrl, int targ, QMatrix op);
void applyReferenceOp(QVector &state, int ctrl, int* targs, int numTargs, QMatrix op);
void applyReferenceOp(QVector &state, int ctrl, int targ1, int targ2, QMatrix op);
void applyReferenceOp(QVector &state, int targ, QMatrix op);

/* comparing quregs to reference data-types */
bool areEqual(Qureg qureg1, Qureg qureg2);
bool areEqual(Qureg qureg, QVector vec);
bool areEqual(Qureg qureg, QMatrix matr);
bool areEqual(Qureg qureg1, Qureg qureg2, qreal precision);
bool areEqual(Qureg qureg, QVector vec, qreal precision);
bool areEqual(Qureg qureg, QMatrix matr, qreal precision);

/** returns log2 of numbers which must be gauranteed to be 2^n */
unsigned int calcLog2(unsigned int res);

/* generating qubit lists */
template<class T> using CatchGen = Catch::Generators::GeneratorWrapper<T>;
CatchGen<int*> sublists(int* list, int len, int sublen);
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps, const int* exclude, int numExclude);
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps, int excluded);
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps);
CatchGen<int*> bitsets(int numBits);
CatchGen<int*> sequences(int base, int numDigits);
CatchGen<pauliOpType*> pauliseqs(int numPaulis);

#endif // QUEST_MATRICES_H
