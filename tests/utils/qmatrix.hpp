#ifndef QMATRIX_HPP
#define QMATRIX_HPP

#include "quest.h"
#include "qvector.hpp"

#include <vector>
using std::vector;


typedef vector<vector<qcomp>> qmatrix;


qmatrix getZeroMatrix(size_t dim);
qmatrix getIdentityMatrix(size_t dim);
qmatrix getDiagonalMatrix(qvector v);

qmatrix operator * (const qcomp&,   const qmatrix& );
qmatrix operator * (const qmatrix&, const qcomp&);
qmatrix operator * (const qreal&,   const qmatrix&);
qmatrix operator * (const qmatrix&, const qreal&);
qmatrix operator *= (qmatrix&, const qcomp&);
qmatrix operator *= (qmatrix&, const qreal&);

qmatrix operator / (const qmatrix&, const qcomp&);
qmatrix operator / (const qmatrix&, const qreal&) ;
qmatrix operator /= (qmatrix&, const qcomp&);
qmatrix operator /= (qmatrix&, const qreal&);

qmatrix operator + (const qmatrix&, const qmatrix&);
qmatrix operator += (qmatrix&, const qmatrix&);

qmatrix operator - (const qmatrix&, const qmatrix&);
qmatrix operator -= (qmatrix&, const qmatrix&);

qmatrix operator * (const qmatrix&, const qmatrix&);
qmatrix operator *= (qmatrix&, const qmatrix&);

void setSubMatrix(qmatrix &dest, qmatrix sub, size_t r, size_t c);
void setToDebugState(qmatrix &m);

qvector getDiagonals(qmatrix m);


#endif // QMATRIX_HPP