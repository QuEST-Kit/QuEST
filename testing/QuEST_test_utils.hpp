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

void applyUnitaryOp(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
);
void applyUnitaryOp(
    QMatrix &state, int *targs, int numTargs, QMatrix op
);
void applyUnitaryOp(
    QMatrix &state, int ctrl, int targ, QMatrix op
);
void applyUnitaryOp(
    QMatrix &state, int targ, QMatrix op
);
void applyUnitaryOp(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
);
void applyUnitaryOp(
    QVector &state, int *targs, int numTargs, QMatrix op
);
void applyUnitaryOp(
    QVector &state, int ctrl, int targ, QMatrix op
);
void applyUnitaryOp(
    QVector &state, int targ, QMatrix op
);

bool areEqual(Qureg qureg, QVector vec);

bool areEqual(Qureg qureg, QMatrix matr);





// DEBUG: DO NOT COMMIT 
void printMatrix(QMatrix a);
void printVector(QVector a);
void printQureg(Qureg q);
void printDif(QMatrix a, Qureg q);




#endif // QUEST_MATRICES_H
