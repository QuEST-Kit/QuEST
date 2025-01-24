#ifndef LINALG_HPP
#define LINALG_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


int log2(qindex);
qindex pow2(int);
qcomp expI(qreal);

qvector getNormalised(qvector);
qvector getDisceteFourierTransform(qvector);

qcomp   getInnerProduct(const qvector &bra, const qvector& ket);
qmatrix getOuterProduct(const qvector &ket, const qvector &bra);

qvector operator * (const qmatrix&, const qvector&);

bool isDiagonal(qmatrix);
bool isUnitary(qmatrix);

qmatrix getTranspose(qmatrix);
qmatrix getConjugateTranspose(qmatrix);
qmatrix getExponentialOfDiagonalMatrix(qmatrix);
qmatrix getExponentialOfPauliMatrix(qreal, qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getKroneckerProduct(qmatrix, qmatrix);

bool isCompletelyPositiveTracePreserving(vector<qmatrix>);


#endif // LINALG_HPP