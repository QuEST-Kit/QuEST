/** @file
 * Testing utilities which perform linear algebra
 * routines upon reference qvector and qmatrix. 
 * These are slow, serial, un-optimised, defensively-
 * designed routines.
 *
 * @author Tyson Jones
 */

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "linalg.hpp"
#include "macros.hpp"
#include "compare.hpp"

#include <algorithm>
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


int getBitAt(qindex num, int ind) {
    DEMAND( num >= 0 );

    return (num >> ind) & 1;
}


vector<int> getBits(qindex num, int numBits) {
    DEMAND( numBits > 0 );

    // out ordered least to most significant
    vector<int> out(numBits);

    for (int i=0; i<numBits; i++)
        out[i] = getBitAt(num, i);

    return out;
}


qindex getBitsAt(qindex num, vector<int> inds) {
    DEMAND( num >= 0 );

    qindex out = 0;

    for (size_t i=0; i<inds.size(); i++)
        out |= getBitAt(num, inds[i]) << i;
    
    return out;
}


qindex setBitAt(qindex num, int ind, int bit) {
    DEMAND( num >= 0 );

    qindex one = 1;
    return (num & ~(one << ind)) | (bit << ind);
}


qindex setBitsAt(qindex num, vector<int> inds, qindex bits) {

    for (size_t i=0; i<inds.size(); i++)
        num = setBitAt(num, inds[i], getBitAt(bits, i));

    return num;
}


int getNumPermutations(int n, int k) {
    DEMAND( n >= k );
    DEMAND( n <= 11 ); // else int overflow

    // P(n, k) = n! / (n-k)!
    qindex p = 1;
    for (int t=n-k+1; t<=n; t++)
        p *= t;

    return p;
}



/*
 * VECTOR OPERATIONS
 */


qcomp getSum(qvector vec) {

    qcomp sum = qcomp(0,0);
    qcomp y, t, c=sum;
    
    // complex Kahan summation
    for (auto& x : vec) {
        y = x - c;
        t = sum + y;
        c = ( t - sum ) - y;
        sum = t;
    }

    return sum;
}


qreal getSum(vector<qreal> vec) {

    // in = real(vec)
    qvector in = getZeroVector(vec.size());
    for (size_t i=0; i<in.size(); i++)
        in[i] = qcomp(vec[i],0);

    return std::real(getSum(in));
}


qvector getNormalised(qvector vec) {

    // prob[i] = abs(vec[i])^2
    vector<qreal> probs(vec.size());
    for (size_t i=0; i<vec.size(); i++)
        probs[i] = std::norm(vec[i]);

    // normalise vector
    qreal norm = getSum(probs);
    qreal fac = 1 / std::sqrt(norm);
    for (auto& x : vec)
        x *= fac;

    return vec;
}


qvector getDisceteFourierTransform(qvector in) {
    DEMAND( in.size() > 0 );
    
    size_t dim = in.size();
    qvector out = getZeroVector(dim);

    // PI must be accurate here
    qreal pi = 3.14159265358979323846;
    qreal a = 1 / std::sqrt(dim);
    qreal b = 2 * pi / dim;
    
    for (size_t x=0; x<dim; x++)
        for (size_t y=0; y<dim; y++)
            out[x] += a * std::exp(b * x * y * 1_i) * in[y];

    return out;
}


qvector getDisceteFourierTransform(qvector in, vector<int> targs) {
    DEMAND( in.size() > 0 );

    size_t dim = in.size();
    qvector out = getZeroVector(dim);
    
    qindex len = getPow2(targs.size());
    qreal pi = 3.14159265358979323846;
    qreal a = 1 / std::sqrt(len);
    qreal b = 2 * pi / len;

    for (size_t i=0; i<dim; i++) {
        size_t x = getBitsAt(i, targs);
        for (size_t y=0; y<len; y++) {
            qindex j = setBitsAt(i, targs, y);
            out[j] += a * std::exp(b * x * y * 1_i) * in[i];
        }
    }

    return out;
}



/*
 * VECTOR & VECTOR OPERATIONS
 */


qcomp getInnerProduct(qvector bra, qvector ket) {
    DEMAND( bra.size() == ket.size() );

    qcomp out = 0;

    for (size_t i=0; i<bra.size(); i++)
        out += std::conj(bra[i]) * ket[i];

    return out;
}


qmatrix getOuterProduct(qvector ket, qvector bra) {
    DEMAND( bra.size() == ket.size() );

    qmatrix out = getZeroMatrix(bra.size());

    for (size_t i=0; i<ket.size(); i++)
        for (size_t j=0; j<ket.size(); j++)
            out[i][j] = ket[i] * std::conj(bra[j]);

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
    qmatrix id = getIdentityMatrix(m.size());
    return doMatricesAgree(md, id);
}


qcomp getTrace(qmatrix m) {

    qcomp out = 0;
    for (size_t r=0; r<m.size(); r++)
        out += m[r][r];

    return out;
}


qmatrix getTranspose(qmatrix m) {

    qmatrix out = getZeroMatrix(m.size());

    for (size_t r=0; r<m.size(); r++)
        for (size_t c=0; c<m.size(); c++)
            out[r][c] = m[c][r];

    return out;
}


qmatrix getConjugate(qmatrix m) {

    for (auto& row : m)
        for (auto& elem : row)
            elem = std::conj(elem);

    return m;
}


qmatrix getConjugateTranspose(qmatrix m) {
    DEMAND( m.size() > 0 );

    // unlike most functions which assume qmatrix
    // is square, this one cheekily handles when
    // 'm' is non-square, since necessary for
    // computing partial traces

    qmatrix out(m[0].size(), qvector(m.size()));

    for (size_t r=0; r<out.size(); r++)
        for (size_t c=0; c<out[0].size(); c++)
            out[r][c] = std::conj(m[c][r]);

    return out;
}


qmatrix getPowerOfDiagonalMatrix(qmatrix m, qcomp p) {
    DEMAND( isDiagonal(m) );

    qmatrix out = getZeroMatrix(m.size());

    // pow(qcomp,qcomp) introduces wildly erroneous
    // imaginary components when both base is real
    // and negative, and exponent is real and integer
    // (so ergo does not produce complex numbers).
    // We divert to real-pow in that scenario!

    for (size_t i=0; i<m.size(); i++) {
        bool mIsRe  = std::imag(m[i][i]) == 0;
        bool mIsNeg = std::real(m[i][i]) < 0;
        bool pIsRe  = std::imag(p) == 0;
        bool pIsInt = std::trunc(std::real(p)) == std::real(p);

        // use pow(qreal,qreal) or pow(qcomp,qcomp)
        out[i][i] = (mIsRe && mIsNeg && pIsRe && pIsInt)?
            qcomp(std::pow(std::real(m[i][i]), std::real(p)),0):
            std::pow(m[i][i], p);
    }

    return out;
}


qmatrix getExponentialOfDiagonalMatrix(qmatrix m) {
    DEMAND( isDiagonal(m) );

    qmatrix out = getZeroMatrix(m.size());

    for (size_t i=0; i<m.size(); i++)
        out[i][i] = std::exp(m[i][i]);

    return out;
}


qmatrix getExponentialOfPauliMatrix(qreal arg, qmatrix m) {
    
    // exp(-i arg/2 m) where m = prod(paulis)
    qmatrix id = getIdentityMatrix(m.size());
    qmatrix out = std::cos(arg/2)*id - 1_i*std::sin(arg/2)*m;
    return out;
}


qmatrix getExponentialOfNormalisedPauliVector(qreal arg, qreal x, qreal y, qreal z) {

    // exp(-arg/2 i [x^ X + y^ Y + z^ Z])
    qreal n = std::sqrt(x*x + y*y + z*z);
    x /= n;
    y /= n;
    z /= n;

    qmatrix id = getIdentityMatrix(2);
    qmatrix out = std::cos(arg/2)*id - 1_i*std::sin(arg/2)*(
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


qmatrix getProjector(int outcome) {
    DEMAND( outcome == 0 || outcome == 1 );

    qmatrix out = getZeroMatrix(2);
    out[outcome][outcome] = 1.;
    return out;
}


qmatrix getProjector(vector<int> targets, vector<int> outcomes, int numQubits) {
    DEMAND( targets.size() == outcomes.size() );
    DEMAND( numQubits > *std::max_element(targets.begin(), targets.end()) );

    // prepare { |0><0|, I, I, |1><1|, ... }
    vector<qmatrix> matrices(numQubits, getIdentityMatrix(2));
    for (size_t i=0; i<targets.size(); i++)
        matrices[targets[i]] = getProjector(outcomes[i]);

    return getKroneckerProduct(matrices);
}


qmatrix getPartialTrace(qmatrix in, vector<int> targets) {
    DEMAND( in.size() > getPow2(targets.size()) );

    auto numTargs = targets.size();
    auto numQubits = getLog2(in.size());
    auto numTargVals = getPow2(numTargs);

    qmatrix out = getZeroMatrix(getPow2(numQubits - numTargs));

    for (qindex v=0; v<numTargVals; v++) {

        // prepare { |0>, I, I, |1>, ... }
        vector<qmatrix> matrices(numQubits, getIdentityMatrix(2));
        for (size_t t=0; t<numTargs; t++) {
            int bit = getBitAt(v, t);
            matrices[targets[t]] = {
                {bit? 0.:1.}, 
                {bit? 1.:0.}};
        }

        qmatrix ket = getKroneckerProduct(matrices);
        qmatrix bra = getConjugateTranspose(ket);
        out += bra * in * ket;
    }

    return out;
}


qmatrix getControlledMatrix(qmatrix matrix, int numCtrls) {

    size_t dim = getPow2(numCtrls) * matrix.size();
    size_t off = dim - matrix.size();
    
    qmatrix out = getIdentityMatrix(dim);
    setSubMatrix(out, matrix, off, off);

    return out;
}


qmatrix getMixture(vector<qmatrix> densmatrs, vector<qreal> probs) {
    DEMAND( densmatrs.size() > 0 );

    qmatrix out = getZeroMatrix(densmatrs[0].size());
    for (size_t i=0; i<densmatrs.size(); i++)
        out += probs[i] * densmatrs[i];

    return out;
}


qmatrix getMixture(vector<qvector> statevecs, vector<qreal> probs) {

    vector<qmatrix> densmatrs(statevecs.size());
    for (size_t i=0; i<statevecs.size(); i++)
        densmatrs[i] = getOuterProduct(statevecs[i], statevecs[i]);

    return getMixture(densmatrs, probs);
}


qmatrix getSuperOperator(vector<qmatrix> matrices) {
    DEMAND( matrices.size() > 0 );

    size_t dim = matrices[0].size();

    // out = sum_m conj(m) (x) m
    qmatrix out = getZeroMatrix(dim * dim);
    for (auto& matr : matrices)
        out += getKroneckerProduct(getConjugate(matr), matr);

    return out;
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

    // we permit the matrices to be non-square which is
    // pretty cheeky (since qmatrix is assumed square with
    // a 2^N dimension by most other functions), but is
    // necessary for us to compute partial traces

    size_t aRows = a.size();
    size_t bRows = b.size();
    size_t aCols = a[0].size();
    size_t bCols = b[0].size();

    qmatrix out(aRows * bRows, qvector(aCols * bCols));
    
    for (size_t r=0; r<bRows; r++)
        for (size_t c=0; c<bCols; c++)
            for (size_t i=0; i<aRows; i++)
                for (size_t j=0; j<aCols; j++)
                    out[r+bRows*i][c+bCols*j] = a[i][j] * b[r][c];

    return out;
}


qmatrix getKroneckerProduct(vector<qmatrix> matrices) {

    qmatrix out = getIdentityMatrix(1);

    // matrices[n-1] (x) ... (x) matrices[0]
    for (auto& m : matrices)
        out = getKroneckerProduct(m, out);

    return out;
}


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


bool isApproxCPTP(vector<qmatrix> matrices) {
    DEMAND( matrices.size() >= 1 );

    size_t dim = matrices[0].size();
    qmatrix id = getIdentityMatrix(dim);
    qmatrix sum = getZeroMatrix(dim);

    for (auto& m : matrices)
        sum += getConjugateTranspose(m) * m;

    return doMatricesAgree(sum, id);
}
