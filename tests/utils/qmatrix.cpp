/** @file
 * Testing utilities which define 'qmatrix', used
 * to perform reference complex matrix algebra, and
 * as a reference proxy to a quantum density matrix.
 *
 * @author Tyson Jones
 */

#include "qmatrix.hpp"
#include "qvector.hpp"
#include "macros.hpp"


/*
 * MATRIX CREATION
 */

qmatrix getZeroMatrix(size_t dim) {
    DEMAND( dim >= 1 );

    return qmatrix(dim, qvector(dim, 0));
}

qmatrix getConstantMatrix(size_t dim, qcomp elem) {
    DEMAND( dim >= 1 );

    return qmatrix(dim, qvector(dim, elem));
}

qmatrix getIdentityMatrix(size_t dim) {
    DEMAND( dim >= 1 );

    qmatrix out = getZeroMatrix(dim);

    for (size_t i=0; i<dim; i++)
        out[i][i] = 1;
    
    return out;
}

qmatrix getDiagonalMatrix(qvector v) {
    DEMAND( v.size() >= 1 );

    qmatrix out = getZeroMatrix(v.size());

    for (size_t i=0; i<v.size(); i++)
        out[i][i] = v[i];
    
    return out;
}

qmatrix getPauliMatrix(int id) {
    DEMAND( id >= 0 );
    DEMAND( id <= 3 );

    if (id == 0)
        return {{1 ,0}, {0, 1}};

    if (id == 1)
        return {{0, 1}, {1, 0}};

    if (id == 2)
        return {{0, -1_i}, {1_i, 0}};

    if (id == 3)
        return {{1, 0}, {0, -1}};

    // unreachable
    return {{-1}};
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
    DEMAND( std::abs(a) != 0 );

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
    DEMAND( m1.size() > 0 );
    DEMAND( m2.size() > 0 );
    DEMAND( m1[0].size() == m2.size() );

    // unlike most functions which assume qmatrix is square,
    // we cheekily permit m1 and m2 to be non-square (and
    // ergo differing), necessary to calculate partial traces
    qmatrix out = qmatrix(m1.size(), qvector(m2[0].size(),0));

    for (size_t r=0; r<out.size(); r++)
        for (size_t c=0; c<out[0].size(); c++)
            for (size_t k=0; k<m2.size(); k++)
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
    DEMAND( sub.size() > 0 );
    DEMAND( sub   .size() + r <= dest.size() );
    DEMAND( sub[0].size() + c <= dest.size() );

    // this function cheekily supports when 'sub' is non-square,
    // which is inconsistent with the preconditions of most of
    // the qmatrix functions, but is needed by setDensityQuregAmps()

    for (size_t i=0; i<sub.size(); i++)
        for (size_t j=0; j<sub[i].size(); j++)
            dest[r+i][c+j] = sub[i][j];
}

void setSubMatrix(qmatrix &dest, qvector sub, size_t flatInd) {
    DEMAND( sub.size() + flatInd <= dest.size()*dest.size() );

    for (size_t i=0; i<sub.size(); i++) {
        size_t r = (i + flatInd) / dest.size(); // floors
        size_t c = (i + flatInd) % dest.size();
        dest[r][c] = sub[i];
    }
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
