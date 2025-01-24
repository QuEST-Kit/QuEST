#ifndef QVECTOR_HPP
#define QVECTOR_HPP

#include "quest.h"
#include "macros.hpp"
#include <vector>


typedef std::vector<qcomp> qvector;


qvector getZeroVector(size_t dim);

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


#endif // QVECTOR_HPP