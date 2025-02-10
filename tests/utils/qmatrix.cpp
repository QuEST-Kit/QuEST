#include "qmatrix.hpp"
#include "qvector.hpp"
#include "macros.hpp"


/*
 * MATRIX CREATION
 */

qmatrix getZeroMatrix(size_t dim) {
    DEMAND( dim > 1 );

    qmatrix out = qmatrix(dim);

    for (auto& row : out)
        row.resize(dim);

    return out;
}

qmatrix getIdentityMatrix(size_t dim) {
    DEMAND( dim > 1 );

    qmatrix out = getZeroMatrix(dim);

    for (size_t i=0; i<dim; i++)
        out[i][i] = 1;
    
    return out;
}

qmatrix getDiagonalMatrix(qvector v) {
    DEMAND( v.size() > 1 );

    qmatrix out = getZeroMatrix(v.size());

    for (size_t i=0; i<v.size(); i++)
        out[i][i] = v[i];
    
    return out;
}


/*
 * SCALAR MULTIPLICATION
 */

qmatrix operator * (const qcomp& a, const qmatrix& m) {
    qmatrix out = m;

    for (auto& row : out)
        for (auto& elem : row)
            elem *= a;
    
    return out;
}

qmatrix operator * (const qmatrix& m, const qcomp& a) {
    return a * m;
}

qmatrix operator *= (qmatrix& m, const qcomp& a) {
    m = a * m;
    return m;
}

qmatrix operator * (const qreal& a, const qmatrix& m) {
    return qcomp(a,0) * m;
}

qmatrix operator * (const qmatrix& m, const qreal& a) {
    return a * m;
}

qmatrix operator *= (qmatrix& m, const qreal& a) {
    m = a * m;
    return m;
}


/*
 * SCALAR DIVISION
 */

qmatrix operator / (const qmatrix& m, const qcomp& a) {
    DEMAND( abs(a) != 0 );

    return (1/a) * m;
}

qmatrix operator /= (qmatrix& m, const qcomp& a) {
    m = m / a;
    return m;
}

qmatrix operator / (const qmatrix& m, const qreal& a) {
    return m / qcomp(a,0);
}

qmatrix operator /= (qmatrix& m, const qreal& a) {
    m = m / a;
    return m;
}


/*
 * MATRIX ADDITION
 */

qmatrix operator + (const qmatrix& m1, const qmatrix& m2) {
    DEMAND( m1.size() == m2.size() );

    qmatrix out = m1;

    for (size_t r=0; r<m1.size(); r++)
        for (size_t c=0; c<m1.size(); c++)
            out[r][c] += m2[r][c];

    return out;
}

qmatrix operator += (qmatrix& m1, const qmatrix& m2) {
    m1 = m1 + m2;
    return m1;
}


/*
 * MATRIX SUBTRACTION
 */

qmatrix operator - (const qmatrix& m1, const qmatrix& m2) {
    return m1 + (-1 * m2);
}

qmatrix operator -= (qmatrix& m1, const qmatrix& m2) {
    m1 = m1 - m2;
    return m1;
}


/*
 * MATRIX MULTIPLICATION
 */

qmatrix operator * (const qmatrix& m1, const qmatrix& m2) {
    DEMAND( m1.size() == m2.size() );

    qmatrix out = getZeroMatrix(m1.size());

    for (size_t r=0; r<m1.size(); r++)
        for (size_t c=0; c<m1.size(); c++)
            for (size_t k=0; k<m1.size(); k++)
                out[r][c] += m1[r][k] * m2[k][c];
    
    return out;
}

qmatrix operator *= (qmatrix& m1, const qmatrix& m2) {
    m1 = m1 * m2;
    return m1;
}


/*
 * SETTERS
 */

void setSubMatrix(qmatrix &dest, qmatrix sub, size_t r, size_t c) {
    DEMAND( sub.size() + r <= dest.size() );
    DEMAND( sub.size() + c <= dest.size() );

    for (size_t i=0; i<sub.size(); i++)
        for (size_t j=0; j<sub.size(); j++)
            dest[r+i][c+j] = sub[i][j];
}

void setToDebugState(qmatrix &m) {
    DEMAND( !m.empty() );

    size_t i = 0;

    // iterate column-wise
    for (size_t c=0; c<m.size(); c++) {
        for (size_t r=0; r<m.size(); r++) {
            m[r][c] = qcomp(2*i/10., (2*i+1)/10.);
            i++;
        }
    }
}


/*
 * GETTERS
 */

qvector getDiagonals(qmatrix m) {
    DEMAND( m.size() > 1 );

    qvector out = getZeroVector(m.size());

    for (size_t i=0; i<out.size(); i++)
        out[i] = m[i][i];
    
    return out;
}
