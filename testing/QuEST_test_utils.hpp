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

/* converting data-types */
QVector toQVector(Qureg qureg);
QMatrix toQMatrix(Complex alpha, Complex beta);
QMatrix toQMatrix(Qureg qureg);
ComplexMatrix2 toComplexMatrix2(QMatrix qm);
ComplexMatrix4 toComplexMatrix4(QMatrix qm);
void toComplexMatrixN(QMatrix qm, ComplexMatrixN cm);

/* building unitary matrices */
QMatrix getMatrixSum(QMatrix a, QMatrix b);
QMatrix getScalarMatrixProduct(qcomp scalar, QMatrix matr);
QMatrix getExponentialDiagonalMatrix(QMatrix a);
QMatrix getKroneckerProduct(QMatrix a, QMatrix b);
QMatrix getFullOperatorMatrix(int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits); 

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
bool areEqual(Qureg qureg, QVector vec);
bool areEqual(Qureg qureg, QMatrix matr);

/* generating random inputs */
QMatrix getRandomUnitary(int numQb);
qreal getRandomReal(qreal min, qreal max);
int getRandomInt(int min, int max);

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
