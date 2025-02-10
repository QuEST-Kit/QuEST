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

qcomp   getInnerProduct(const qvector &bra, const qvector& ket);
qmatrix getOuterProduct(const qvector &ket, const qvector &bra);

qvector operator * (const qmatrix&, const qvector&);

bool isDiagonal(qmatrix);
bool isUnitary(qmatrix);

qmatrix getTranspose(qmatrix);
qmatrix getConjugateTranspose(qmatrix);
qmatrix getPowerOfDiagonalMatrix(qmatrix, qcomp);
qmatrix getExponentialOfDiagonalMatrix(qmatrix);
qmatrix getExponentialOfPauliMatrix(qreal, qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getKroneckerProduct(qmatrix, qmatrix);

bool isCompletelyPositiveTracePreserving(vector<qmatrix>);


#endif // LINALG_HPP