#ifndef LINALG_HPP
#define LINALG_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


int getLog2(qindex);
qindex getPow2(int);
qcomp getExpI(qreal);

qvector getNormalised(qvector);
qvector getDisceteFourierTransform(qvector);

qcomp   getInnerProduct(const qvector& bra, const qvector& ket);
qmatrix getOuterProduct(const qvector& ket, const qvector& bra);

qvector operator * (const qmatrix&, const qvector&);

bool isDiagonal(qmatrix);
bool isApproxUnitary(qmatrix);

qmatrix getTranspose(qmatrix);
qmatrix getConjugateTranspose(qmatrix);
qmatrix getPowerOfDiagonalMatrix(qmatrix diag, qcomp power);
qmatrix getExponentialOfDiagonalMatrix(qmatrix diag);
qmatrix getExponentialOfPauliMatrix(qreal arg, qmatrix pauli);
qmatrix getExponentialOfNormalisedPauliVector(qreal arg, qreal x, qreal y, qreal z);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getKroneckerProduct(qmatrix, qmatrix);
qmatrix getKroneckerProduct(qmatrix, int);

bool isCompletelyPositiveTracePreserving(vector<qmatrix>);


#endif // LINALG_HPP