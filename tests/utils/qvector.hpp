/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsqvector qvector
 * @ingroup testutils
 * @brief
 * Testing utilities which define 'qvector', used
 * as a reference proxy to a quantum statevector.
 * @{
 */

#ifndef QVECTOR_HPP
#define QVECTOR_HPP

#include "quest/include/quest.h"
#include "macros.hpp"
#include <vector>


typedef std::vector<qcomp> qvector;


qvector getZeroVector(size_t dim);
qvector getConstantVector(size_t dim, qcomp elem);

qvector operator * (const qcomp&, const qvector&);
qvector operator * (const qvector&, const qcomp&);
qvector operator * (const qreal&, const qvector&);
qvector operator * (const qvector&, const qreal&);
qvector operator *= (qvector&, const qcomp&);
qvector operator *= (qvector&, const qreal&);

qvector operator / (const qvector&, const qcomp&);
qvector operator / (const qvector&, const qreal&);
qvector operator /= (qvector&, const qcomp&);
qvector operator /= (qvector&, const qreal&);

qvector operator + (const qvector&, const qvector&);
qvector operator += (qvector&, const qvector&);

qvector operator - (const qvector&, const qvector&);
qvector operator -= (qvector&, const qvector&);

void setSubVector(qvector &dest, qvector sub, size_t i);
void setToDebugState(qvector &v);


#endif // QVECTOR_HPP

/** @} (end defgroup) */
