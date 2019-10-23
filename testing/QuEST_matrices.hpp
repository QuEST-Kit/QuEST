/** @file
 * Simple implementations of matrix operations used by QuEST_unit_tests
 *
 * @author Tyson Jones
 */

#ifndef QUEST_MATRICES_H
#define QUEST_MATRICES_H

#include "QuEST.h"
#include "QuEST_complex.h"
#include <vector>

typedef std::vector<std::vector<qcomp>> QMatrix;
typedef std::vector<qcomp> QVector;

QMatrix toQMatrix(Complex alpha, Complex beta);

void applyQUnitary(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
);
void applyQUnitary(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
);

bool areEqual(QVector vec, Qureg qureg);

bool areEqual(QMatrix matr, Qureg qureg);

#endif // QUEST_MATRICES_H
