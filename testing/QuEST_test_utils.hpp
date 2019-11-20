/** @file
 * Simple implementations of matrix operations used by QuEST_unit_tests
 *
 * @author Tyson Jones
 */

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

/* applying unitary operations to reference state-vectors and density matrices */
void applyUnitaryOp(QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);
void applyUnitaryOp(QMatrix &state, int* ctrls, int numCtrls, int target, QMatrix op);
void applyUnitaryOp(QMatrix &state, int *targs, int numTargs, QMatrix op);
void applyUnitaryOp(QMatrix &state, int ctrl, int targ, QMatrix op);
void applyUnitaryOp(QMatrix &state, int targ, QMatrix op);
void applyUnitaryOp(QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);
void applyUnitaryOp(QVector &state, int* ctrls, int numCtrls, int target, QMatrix op);
void applyUnitaryOp(QVector &state, int *targs, int numTargs, QMatrix op);
void applyUnitaryOp(QVector &state, int ctrl, int targ, QMatrix op);
void applyUnitaryOp(QVector &state, int targ, QMatrix op);

/* comparing quregs to reference data-types */
bool areEqual(Qureg qureg, QVector vec);
bool areEqual(Qureg qureg, QMatrix matr);

/* generating random inputs */
QMatrix getRandomUnitary(int numQb);

/* generating qubit lists */
template<class T> using CatchGen = Catch::Generators::GeneratorWrapper<T>;
CatchGen<int*> sublists(int* list, int len, int sublen);
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps, int* exclude, int numExclude);
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps, int excluded);

#endif // QUEST_MATRICES_H
