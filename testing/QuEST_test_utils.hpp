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

QVector toQVector(Qureg qureg);

QMatrix toQMatrix(Complex alpha, Complex beta);

QMatrix toQMatrix(Qureg qureg);

void applyQUnitary(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
);
void applyQUnitary(
    QMatrix &state, int *targs, int numTargs, QMatrix op, int numQubits
);
void applyQUnitary(
    QMatrix &state, int ctrl, int targ, QMatrix op, int numQubits
);
void applyQUnitary(
    QMatrix &state, int targ, QMatrix op, int numQubits
);
void applyQUnitary(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
);
void applyQUnitary(
    QVector &state, int *targs, int numTargs, QMatrix op, int numQubits
);
void applyQUnitary(
    QVector &state, int ctrl, int targ, QMatrix op, int numQubits
);
void applyQUnitary(
    QVector &state, int targ, QMatrix op, int numQubits
);

bool areEqual(QVector vec, Qureg qureg);

bool areEqual(QMatrix matr, Qureg qureg);

#endif // QUEST_MATRICES_H
