#include "qvector.hpp"
#include "qmatrix.hpp"
#include "linalg.hpp"
#include "macros.hpp"

#include <vector>
using std::vector;



/*
 * SCALAR OPERATIONS
 */


int getLog2(qindex a) {
    DEMAND( a >= 0 );
    DEMAND( (a & (a - 1)) == 0 );  // is pow2

    int n = 0;
    while (a >>= 1)
        n++;

    return n;
}


qindex getPow2(int a) {
    DEMAND( a >= 0 );

    return ((qindex) 1) << a;
}


qcomp getExpI(qreal x) {
    return qcomp(cos(x), sin(x));
}



/*
 * VECTOR OPERATIONS
 */


qvector getNormalised(qvector vec) {

    qreal norm = 0;
    qreal y, t, c=0;
    
    // compute norm via Kahan summation
    for (auto& x : vec) {
        y = real(x)*real(x) - c;
        t = norm + y;
        c = ( t - norm ) - y;
        norm = t;
        
        y = imag(x)*imag(x) - c;
        t = norm + y;
        c = ( t - norm ) - y;
        norm = t;
    }
    
    // normalise vector
    for (auto& x : vec)
        x /= sqrt(norm);

    return vec;
}


qvector getDisceteFourierTransform(qvector in) {
    REQUIRE( in.size() > 0 );
    
    size_t dim = in.size();
    qvector out = getZeroVector(dim);

    qreal a = 1 / sqrt(dim);
    qreal b = 2 * M_PI / dim;
    
    for (size_t x=0; x<dim; x++)
        for (size_t y=0; y<dim; y++)
            out[x] += a * getExpI(b * x * y) * in[y];

    return out;
}



/*
 * VECTOR & VECTOR OPERATIONS
 */


qcomp getInnerProduct(const qvector &bra, const qvector& ket) {
    DEMAND( bra.size() == ket.size() );

    qcomp out = 0;

    for (size_t i=0; i<bra.size(); i++)
        out += conj(bra[i]) * ket[i];

    return out;
}


qmatrix getOuterProduct(const qvector &ket, const qvector &bra) {
    DEMAND( bra.size() == ket.size() );

    qmatrix out = getZeroMatrix(bra.size());

    for (size_t i=0; i<ket.size(); i++)
        for (size_t j=0; j<ket.size(); j++)
            out[i][j] = ket[i] * conj(bra[j]);

    return out;
}



/*
 * MATRIX OPERATIONS
 */


bool isDiagonal(qmatrix m) {

    for (size_t r=0; r<m.size(); r++)
        for (size_t c=0; c<m.size(); c++)
            if (r!=c && m[r][c] != 0_i)
                return false;

    return true;
}


bool isApproxUnitary(qmatrix m) {

    // should be identity
    qmatrix md = m * getConjugateTranspose(m);

    for (size_t r=0; r<md.size(); r++)
        for (size_t c=0; c<md.size(); c++)
            if (abs(md[r][c] - (r==c)) > TEST_EPSILON)
                return false;

    return true;
}


qmatrix getTranspose(qmatrix m) {

    qmatrix out = getZeroMatrix(m.size());

    for (size_t r=0; r<m.size(); r++)
        for (size_t c=0; c<m.size(); c++)
            out[r][c] = m[c][r];

    return out;
}


qmatrix getConjugateTranspose(qmatrix m) {

    qmatrix out = getZeroMatrix(m.size());

    for (size_t r=0; r<m.size(); r++)
        for (size_t c=0; c<m.size(); c++)
            out[r][c] = conj(m[c][r]);

    return out;
}


qmatrix getPowerOfDiagonalMatrix(qmatrix m, qcomp p) {
    DEMAND( isDiagonal(m) );

    qmatrix out = getZeroMatrix(m.size());

    for (size_t i=0; i<m.size(); i++)
        out[i][i] = pow(m[i][i], p);

    return out;
}


qmatrix getExponentialOfDiagonalMatrix(qmatrix m) {
    DEMAND( isDiagonal(m) );

    qmatrix out = getZeroMatrix(m.size());

    for (size_t i=0; i<m.size(); i++)
        out[i][i] = exp(m[i][i]);

    return out;
}


qmatrix getExponentialOfPauliMatrix(qreal arg, qmatrix m) {
    
    // exp(-i arg/2 m) where m = prod(paulis)
    qmatrix id = getIdentityMatrix(m.size());
    qmatrix out = cos(arg/2)*id - 1_i*sin(arg/2)*m;
    return out;
}


qmatrix getExponentialOfNormalisedPauliVector(qreal arg, qreal x, qreal y, qreal z) {

    // exp(-arg/2 i [x^ X + y^ Y + z^ Z])
    qreal n = sqrt(x*x + y*y + z*z);
    x /= n;
    y /= n;
    z /= n;

    qmatrix id = getIdentityMatrix(2);
    qmatrix out = cos(arg/2)*id - 1_i*sin(arg/2)*(
        x * getPauliMatrix(1) +
        y * getPauliMatrix(2) +
        z * getPauliMatrix(3));

    return out;
}


qmatrix getOrthonormalisedRows(qmatrix matr) {

    // perform the Gram-Schmidt process, processing each row of matr in-turn
    for (size_t i=0; i<matr.size(); i++) {
        qvector row = matr[i];
        
        // compute new orthogonal row by subtracting proj row onto prevs
        for (int k=i-1; k>=0; k--) {

            // compute inner_product(row, prev) = row . conj(prev)
            qcomp prod = getInnerProduct(matr[k], row);
                            
            // subtract (proj row onto prev) = (prod * prev) from final row
            matr[i] -= prod * matr[k];
        }
    
        // normalise the row 
        matr[i] = getNormalised(matr[i]);
    }
    
    // return the new orthonormal matrix
    return matr;
}



/*
 * MATRIX & VECTOR OPERATIONS
 */


qvector operator * (const qmatrix& m, const qvector& v) {
    DEMAND( m.size() == v.size() );

    qvector out = getZeroVector(v.size());

    for (size_t r=0; r<v.size(); r++)
        for (size_t c=0; c<v.size(); c++)
            out[r] += m[r][c] * v[c];

    return out;
}



/*
 * MATRIX & MATRIX OPERATIONS
 */


qmatrix getKroneckerProduct(qmatrix a, qmatrix b) {

    qmatrix out = getZeroMatrix(a.size() * b.size());

    for (size_t r=0; r<b.size(); r++)
        for (size_t c=0; c<b.size(); c++)
            for (size_t i=0; i<a.size(); i++)
                for (size_t j=0; j<a.size(); j++)
                    out[r+b.size()*i][c+b.size()*j] = a[i][j] * b[r][c];

    return out;
}



/*
 * MATRIX & SCALAR OPERATIONS
 */


qmatrix getKroneckerProduct(qmatrix m, int count) {
    DEMAND( count >= 1 );

    qmatrix out = getIdentityMatrix(1);

    for (int n=0; n<count; n++)
        out = getKroneckerProduct(out, m);

    return out;
}



/*
 * MATRIX COLLECTIONS
 */


bool isCompletelyPositiveTracePreserving(vector<qmatrix> matrices) {
    DEMAND( matrices.size() >= 1 );

    size_t dim = matrices[0].size();
    qmatrix id = getIdentityMatrix(dim);
    qmatrix sum = getZeroMatrix(dim);

    for (auto& m : matrices)
        sum += getConjugateTranspose(m) * m;

    return sum == id;
}
