/** @file
 * Testing utilities which define 'qvector', used
 * as a reference proxy to a quantum statevector.
 *
 * @author Tyson Jones
 */

#include "qvector.hpp"
#include "macros.hpp"


/*
 * VECTOR CREATION
 */

qvector getZeroVector(size_t dim) {
    // permit dim = 0
    
    return qvector(dim, 0);
}

qvector getConstantVector(size_t dim, qcomp elem) {
    DEMAND( dim >= 1 );
    
    return qvector(dim, elem);
}


/*
 * SCALAR MULTIPLICATION
 */

qvector operator * (const qcomp& a, const qvector& v) {
    qvector out = v;

    for (auto& x : out)
        x *= a;

    return out;
}

qvector operator * (const qvector& v, const qcomp& a) {
    return a * v;
}

qvector operator *= (qvector& v, const qcomp& a) {
    v = a * v;
    return v;
}

qvector operator * (const qreal& a, const qvector& v) {
    return qcomp(a,0) * v;
}

qvector operator * (const qvector& v, const qreal& a) {
    return a * v;
}

qvector operator *= (qvector& v, const qreal& a) {
    v = a * v;
    return v;
}


/*
 * SCALAR DIVISION
 */

qvector operator / (const qvector& v, const qcomp& a) {
    DEMAND( std::abs(a) != 0 );

    return (1/a) * v;
}

qvector operator /= (qvector& v, const qcomp& a) {
    v = v / a;
    return v;
}

qvector operator / (const qvector& v, const qreal& a) {
    return v / qcomp(a,0);
}

qvector operator /= (qvector& v, const qreal& a) {
    v = v / a;
    return v;
}


/*
 * VECTOR ADDITION
 */

qvector operator + (const qvector& v1, const qvector& v2) {
    DEMAND( v1.size() == v2.size() );

    qvector out = v1;
    
    for (size_t i=0; i<v2.size(); i++)
        out[i] += v2[i];

    return out;
}

qvector operator += (qvector& v1, const qvector& v2) {
    v1 = v1 + v2;
    return v1;
}


/*
 * VECTOR SUBTRACTION
 */

qvector operator - (const qvector& v1, const qvector& v2) {
    return v1 + (-1 * v2);
}

qvector operator -= (qvector& v1, const qvector& v2) {
    v1 = v1 - v2;
    return v1;
}


/*
 * SETTERS
 */

void setSubVector(qvector &dest, qvector sub, size_t i) {
    DEMAND( sub.size() + i <= dest.size() );

    for (size_t j=0; j<sub.size(); j++)
        dest[j+i] = sub[j];
}

void setToDebugState(qvector &v) {
    DEMAND( !v.empty() );

    for (size_t i=0; i<v.size(); i++)
        v[i] = qcomp(2*i/10., (2*i+1)/10.);
}