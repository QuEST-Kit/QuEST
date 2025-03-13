/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsqmatrix qmatrix
 * @ingroup testutils
 * @brief
 * Testing utilities which define 'qmatrix', used
 * to perform reference complex matrix algebra, and
 * as a reference proxy to a quantum density matrix.
 * @{
 */

#ifndef QMATRIX_HPP
#define QMATRIX_HPP

#include "quest/include/quest.h"
#include "qvector.hpp"

#include <vector>
using std::vector;


typedef vector<vector<qcomp>> qmatrix;


qmatrix getZeroMatrix(size_t dim);
qmatrix getConstantMatrix(size_t dim, qcomp elem);
qmatrix getIdentityMatrix(size_t dim);
qmatrix getDiagonalMatrix(qvector v);
qmatrix getPauliMatrix(int id);

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
void setSubMatrix(qmatrix &dest, qvector sub, size_t flatInd);
void setToDebugState(qmatrix &m);

qvector getDiagonals(qmatrix m);


#endif // QMATRIX_HPP

/** @} (end defgroup) */
