/** @file
 * Unoptimised, analytic implementations of matrix operations used by QuEST_unit_tests
 *
 * @author Tyson Jones
 */

#include "QuEST.h"
#include "utilities.hpp"
#include "catch.hpp"
#include <random>
#include <algorithm>
#include <bitset>

#ifdef DISTRIBUTED_MODE 
#include <mpi.h>
#endif

/* (don't generate doxygen doc) 
 *
 * preconditions to the internal unit testing functions are checked using 
 * DEMAND rather than Catch2's REQUIRE, so that they are not counted in the 
 * total unit testing statistics (e.g. number of checks passed).
 */
#define DEMAND( cond ) if (!(cond)) FAIL( ); 

QVector operator + (const QVector& v1, const QVector& v2) {
    DEMAND( v1.size() == v2.size() );
    QVector out = v1;
    for (size_t i=0; i<v2.size(); i++)
        out[i] += v2[i];
    return out;
}
QVector operator - (const QVector& v1, const QVector& v2) {
    DEMAND( v1.size() == v2.size() );
    QVector out = v1;
    for (size_t i=0; i<v2.size(); i++)
        out[i] -= v2[i];
    return out;
}
QVector operator * (const qcomp& a, const QVector& v) {
    QVector out = v;
    for (size_t i=0; i<v.size(); i++)
        out[i] *= a;
    return out;
}
QVector operator * (const QVector& v, const qcomp& a) {
    return a * v;
}
QVector operator / (const QVector& v, const qcomp& a) {
    DEMAND( abs(a) != 0 );
    QVector out = v;
    for (size_t i=0; i<v.size(); i++)
        out[i] /= a;
    return out;
}
qcomp operator * (const QVector &v1, const QVector& v2) {
    // this is sum_i v1_i conj(v2_i)
    DEMAND( v1.size() == v2.size() );
    qcomp out = 0;
    for (size_t i=0; i<v1.size(); i++)
        out += v1[i] * conj(v2[i]);
    return out;
}
void operator += (QVector& v1, const QVector& v2) { // these violate += returns (fine)
    v1 = v1 + v2;
}
void operator -= (QVector& v1, const QVector& v2) {
    v1 = v1 - v2;
}
void operator *= (QVector& v1, const qcomp& a) {
    v1 = v1 * a;
}
void operator /= (QVector& v1, const qcomp& a) {
    v1 = v1 / a;
}
QMatrix operator + (const QMatrix& m1, const QMatrix& m2) {
    DEMAND( m1.size() == m2.size() );
    QMatrix out = m1;
    for (size_t r=0; r<m1.size(); r++)
        for (size_t c=0; c<m1.size(); c++)
            out[r][c] += m2[r][c];
    return out;
}
QMatrix operator - (const QMatrix& m1, const QMatrix& m2) {
    DEMAND( m1.size() == m2.size() );
    QMatrix out = m1;
    for (size_t r=0; r<m1.size(); r++)
        for (size_t c=0; c<m1.size(); c++)
            out[r][c] -= m2[r][c];
    return out;
}
QMatrix operator * (const qcomp& a, const QMatrix& m) {
    QMatrix out = m;
    for (size_t r=0; r<m.size(); r++)
        for (size_t c=0; c<m.size(); c++)
            out[r][c] *= a;
    return out;
}
QMatrix operator * (const QMatrix& m, const qcomp& a) {
    return a * m;
}
QMatrix operator / (const QMatrix& m, const qcomp& a) {
    QMatrix out = m;
    for (size_t r=0; r<m.size(); r++)
        for (size_t c=0; c<m.size(); c++)
            out[r][c] /= a;
    return out;
}
QMatrix operator * (const QMatrix& m1, const QMatrix& m2) {
    QMatrix prod = m1; // will be completely overwritten
    for (size_t r=0; r<m1.size(); r++)
        for (size_t c=0; c<m1.size(); c++) {
            prod[r][c] = 0;
            for (size_t k=0; k<m1.size(); k++)
                prod[r][c] += m1[r][k] * m2[k][c];
        }
    return prod;
}
void operator += (QMatrix& m1, const QMatrix& m2) {
    m1 = m1 + m2;
}
void operator -= (QMatrix& m1, const QMatrix& m2) {
    m1 = m1 - m2;
}
void operator *= (QMatrix& m1, const qreal& a) {
    m1 = m1 * a;
}
void operator /= (QMatrix& m1, const qreal& a) {
    m1 = m1 / a;
}
void operator *= (QMatrix& m1, const QMatrix& m2) {
    m1 = m1 * m2;
}
QVector operator * (const QMatrix& m, const QVector& v) {
    DEMAND( m.size() == v.size() );
    QVector prod = QVector(v.size());
    for (size_t r=0; r<v.size(); r++)
        for (size_t c=0; c<v.size(); c++)
            prod[r] += m[r][c] * v[c];
    return prod;
}

QVector getKroneckerProduct(QVector b, QVector a) {

    QVector prod = QVector(a.size() * b.size());
    
    for (size_t i=0; i<prod.size(); i++)
        prod[i] = b[i / a.size()] * a[i % a.size()];
    
    return prod;
}

QMatrix getZeroMatrix(size_t dim) {
    DEMAND( dim > 1 );
    QMatrix matr = QMatrix(dim);
    for (size_t i=0; i<dim; i++)
        matr[i].resize(dim);
    return matr;
}

QMatrix getIdentityMatrix(size_t dim) {
    DEMAND( dim > 1 );
    QMatrix matr = getZeroMatrix(dim);
    for (size_t i=0; i<dim; i++)
        matr[i][i] = 1;
    return matr;
}

QMatrix getKetBra(QVector ket, QVector bra) {
    DEMAND( ket.size() == bra.size() );
    QMatrix mat = getZeroMatrix(ket.size());
    
    for (size_t r=0; r<ket.size(); r++)
        for (size_t c=0; c<ket.size(); c++)
            mat[r][c] = ket[r] * conj(bra[c]);
    return mat;
}

QMatrix getKroneckerProduct(QMatrix a, QMatrix b) {
    QMatrix prod = getZeroMatrix(a.size() * b.size());
    for (size_t r=0; r<b.size(); r++)
        for (size_t c=0; c<b.size(); c++)
            for (size_t i=0; i<a.size(); i++)
                for (size_t j=0; j<a.size(); j++)
                    prod[r+b.size()*i][c+b.size()*j] = a[i][j] * b[r][c];
    return prod;
}

QMatrix getConjugateTranspose(QMatrix a) {
    QMatrix b = a;
    for (size_t r=0; r<a.size(); r++)
        for (size_t c=0; c<a.size(); c++)
            b[r][c] = conj(a[c][r]);
    return b;
}

QMatrix getExponentialOfDiagonalMatrix(QMatrix a) {
    
    // ensure diagonal
    for (size_t r=0; r<a.size(); r++)
        for (size_t c=0; c<a.size(); c++) {
            if (r == c)
                continue;
            DEMAND( real(a[r][c]) == 0. );
            DEMAND( imag(a[r][c]) == 0. );
        }

    // exp(diagonal) = diagonal(exp)
    QMatrix diag = a;
    for (size_t i=0; i<a.size(); i++)
        diag[i][i] = exp(diag[i][i]);
        
    return diag;
}

QMatrix getExponentialOfPauliMatrix(qreal angle, QMatrix a) {
    QMatrix iden = getIdentityMatrix(a.size());
    QMatrix expo = (cos(angle/2) * iden) + ((qcomp) -1i * sin(angle/2) * a);
    return expo;
}

void setSubMatrix(QMatrix &dest, QMatrix sub, size_t r, size_t c) {
    DEMAND( sub.size() + r <= dest.size() );
    DEMAND( sub.size() + c <= dest.size() );
    for (size_t i=0; i<sub.size(); i++)
        for (size_t j=0; j<sub.size(); j++)
            dest[r+i][c+j] = sub[i][j];
}

QMatrix getSwapMatrix(int qb1, int qb2, int numQb) {
    DEMAND( numQb > 1 );
    DEMAND( (qb1 >= 0 && qb1 < numQb) );
    DEMAND( (qb2 >= 0 && qb2 < numQb) );
    
    if (qb1 > qb2)
        std::swap(qb1, qb2);
        
    if (qb1 == qb2)
        return getIdentityMatrix(1 << numQb);

    QMatrix swap;
    
    if (qb2 == qb1 + 1) {
        // qubits are adjacent
        swap = QMatrix{{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};
        
    } else {
        // qubits are distant
        int block = 1 << (qb2 - qb1);
        swap = getZeroMatrix(block*2);
        QMatrix iden = getIdentityMatrix(block/2);
        
        // Lemma 3.1 of arxiv.org/pdf/1711.09765.pdf
        QMatrix p0{{1,0},{0,0}};
        QMatrix l0{{0,1},{0,0}};
        QMatrix l1{{0,0},{1,0}};
        QMatrix p1{{0,0},{0,1}};
        
        /* notating a^(n+1) = identity(1<<n) (otimes) a, we construct the matrix
         * [ p0^(N) l1^N ]
         * [ l0^(N) p1^N ]
         * where N = qb2 - qb1 */
        setSubMatrix(swap, getKroneckerProduct(iden, p0), 0, 0);
        setSubMatrix(swap, getKroneckerProduct(iden, l0), block, 0);
        setSubMatrix(swap, getKroneckerProduct(iden, l1), 0, block);
        setSubMatrix(swap, getKroneckerProduct(iden, p1), block, block);
    }
    
    // pad swap with outer identities
    if (qb1 > 0)
        swap = getKroneckerProduct(swap, getIdentityMatrix(1<<qb1));
    if (qb2 < numQb-1)
        swap = getKroneckerProduct(getIdentityMatrix(1<<(numQb-qb2-1)), swap);
        
    return swap;
}

/* (don't generate doxygen doc)
 *
 * iterates list1 (of length len1) and replaces element oldEl with newEl, which is 
 * gauranteed to be present at most once (between list1 AND list2), though may 
 * not be present at all. If oldEl isn't present in list1, does the same for list2. 
 * list1 is skipped if == NULL. This is used by getFullOperatorMatrix() to ensure
 * that, when qubits are swapped, their appearences in the remaining qubit lists 
 * are updated.
 */
void updateIndices(int oldEl, int newEl, int* list1, int len1, int* list2, int len2) {
    if (list1 != NULL) {
        for (int i=0; i<len1; i++) {
            if (list1[i] == oldEl) {
                list1[i] = newEl;
                return;
            }
        }
    }
    for (int i=0; i<len2; i++) {
        if (list2[i] == oldEl) {
            list2[i] = newEl;
            return;
        }
    }
}

QMatrix getFullOperatorMatrix(
    int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
) {        
    DEMAND( numCtrls >= 0 );
    DEMAND( numTargs >= 0 );
    DEMAND( numQubits >= (numCtrls+numTargs) );
    DEMAND( op.size() == (1u << numTargs) );
    
    // copy {ctrls} and {targs}to restore at end
    std::vector<int> ctrlsCopy(ctrls, ctrls+numCtrls);
    std::vector<int> targsCopy(targs, targs+numTargs);
    
    // full-state matrix of qubit swaps
    QMatrix swaps = getIdentityMatrix(1 << numQubits);
    QMatrix unswaps = getIdentityMatrix(1 << numQubits);
    QMatrix matr;
    
    // swap targs to {0, ..., numTargs-1}
    for (int i=0; i<numTargs; i++) {
        if (i != targs[i]) {
            matr = getSwapMatrix(i, targs[i], numQubits);
            swaps = matr * swaps;
            unswaps = unswaps * matr;
            
            // even if this is the last targ, ctrls might still need updating
            updateIndices(
                i, targs[i], (i < numTargs-1)? &targs[i+1] : NULL, 
                numTargs-i-1, ctrls, numCtrls);
        }
    }

    // swap ctrls to {numTargs, ..., numTargs+numCtrls-1}
    for (int i=0; i<numCtrls; i++) {
        int newInd = numTargs+i;
        if (newInd != ctrls[i]) {
            matr = getSwapMatrix(newInd, ctrls[i], numQubits);
            swaps = matr * swaps;
            unswaps = unswaps * matr;
            
            // update remaining ctrls (if any exist)
            if (i < numCtrls-1)
                updateIndices(newInd, ctrls[i], NULL, 0, &ctrls[i+1], numCtrls-i-1);
        }
    }
    
    // construct controlled-op matrix for qubits {0, ..., numCtrls+numTargs-1}
    size_t dim = 1 << (numCtrls+numTargs);
    QMatrix fullOp = getIdentityMatrix(dim);
    setSubMatrix(fullOp, op, dim-op.size(), dim-op.size());
    
    // create full-state controlled-op matrix (left-pad identities)
    if (numQubits > numCtrls+numTargs) {
        size_t pad = 1 << (numQubits - numCtrls - numTargs);
        fullOp = getKroneckerProduct(getIdentityMatrix(pad), fullOp);
    }
    
    // apply swap to either side (to swap qubits back and forth)
    fullOp = unswaps * fullOp * swaps;
    
    // restore {ctrls and targs}
    for (int i=0; i<numCtrls; i++)
        ctrls[i] = ctrlsCopy[i];
    for (int i=0; i<numTargs; i++)
        targs[i] = targsCopy[i];

    return fullOp;
}

unsigned int calcLog2(long unsigned int res) {
    unsigned int n = 0;
    while (res >>= 1)
        n++;
    return n;
}

QMatrix getRandomQMatrix(int dim) {
    DEMAND( dim > 1 );
    
    QMatrix matr = getZeroMatrix(dim);
    for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++) {
            
            // generate 2 normally-distributed random numbers via Box-Muller
            qreal a = rand()/(qreal) RAND_MAX;
            qreal b = rand()/(qreal) RAND_MAX;
            qreal r1 = sqrt(-2 * log(a)) * cos(2 * 3.14159265 * b);
            qreal r2 = sqrt(-2 * log(a)) * sin(2 * 3.14159265 * b);
            
            matr[i][j] = r1 + r2 * (qcomp) 1i;
        }
    }
    return matr;
}

bool areEqual(QVector a, QVector b) {
    DEMAND( a.size() == b.size() );
    
    for (size_t i=0; i<a.size(); i++)
        if (abs(a[i] - b[i]) > REAL_EPS)
            return false;
    return true;
}

bool areEqual(QMatrix a, QMatrix b) {
    DEMAND( a.size() == b.size() );
    
    for (size_t i=0; i<a.size(); i++)
        for (size_t j=0; j<b.size(); j++)
            if (abs(a[i][j] - b[i][j]) > REAL_EPS)
                return false;
    return true;
}

qcomp expI(qreal phase) {
    return qcomp(cos(phase), sin(phase));
}

qreal getRandomReal(qreal min, qreal max) {
    DEMAND( min <= max );
    qreal r = min + (max - min) * (rand() / (qreal) RAND_MAX);
    
    // check bounds satisfied 
    DEMAND( r >= min );
    DEMAND( r <= max );
    return r;
}

qcomp getRandomComplex() {
    return getRandomReal(-1,1) + getRandomReal(-1,1) * (qcomp) 1i;
}

QVector getRandomQVector(int dim) { 
    QVector vec = QVector(dim);
    for (int i=0; i<dim; i++)
        vec[i] = getRandomComplex();
        
    // check we didn't get the impossibly-unlikely zero-amplitude outcome 
    DEMAND( real(vec[0]) != 0 );
        
    return vec;
}

QVector getNormalised(QVector vec) {
    qreal norm = 0;
    qreal y, t, c;
    c = 0;
    
    for (size_t i=0; i<vec.size(); i++) {
        y = real(vec[i])*real(vec[i]) - c;
        t = norm + y;
        c = ( t - norm ) - y;
        norm = t;
        
        y = imag(vec[i])*imag(vec[i]) - c;
        t = norm + y;
        c = ( t - norm ) - y;
        norm = t;
    }
    
    for (size_t i=0; i<vec.size(); i++)
        vec[i] /= sqrt(norm);
    return vec;
}

QVector getRandomStateVector(int numQb) {
    return getNormalised(getRandomQVector(1<<numQb));
}

std::vector<qreal> getRandomProbabilities(int numProbs) {
    
    // generate random unnormalised scalars
    std::vector<qreal> probs;
    qreal total = 0;
    for (int i=0; i<numProbs; i++) {
        qreal prob = getRandomReal(0, 1);
        probs.push_back(prob);
        total += prob;
    }
        
    // normalise
    for (int i=0; i<numProbs; i++)
        probs[i] /= total;
        
    return probs;
}

QMatrix getRandomDensityMatrix(int numQb) {
    DEMAND( numQb > 0 );
    
    // generate random probabilities to weight random pure states
    int dim = 1<<numQb;
    std::vector<qreal> probs = getRandomProbabilities(dim);
    
    // add random pure states
    QMatrix dens = getZeroMatrix(dim);
    for (int i=0; i<dim; i++) {
        QVector pure = getRandomStateVector(numQb);
        dens += probs[i] * getKetBra(pure, pure);
    }
    
    return dens;
}

QMatrix getPureDensityMatrix(QVector state) {
    return getKetBra(state, state);
}

QMatrix getRandomPureDensityMatrix(int numQb) {
    QVector vec = getRandomStateVector(numQb);
    QMatrix mat = getPureDensityMatrix(vec);
    return mat;
}

QVector getMatrixDiagonal(QMatrix matr) {
    
    QVector vec = QVector(matr.size());
    for (size_t i=0; i<vec.size(); i++)
        vec[i] = matr[i][i];
    
    return vec;
}

int getRandomInt(int min, int max) {
    return round(getRandomReal(min, max-1));
}

QMatrix getRandomUnitary(int numQb) {
    DEMAND( numQb >= 1 );

    QMatrix matr = getRandomQMatrix(1 << numQb);

    for (size_t i=0; i<matr.size(); i++) {
        QVector row = matr[i];
        
        // compute new orthogonal row by subtracting proj row onto prevs
        for (int k=i-1; k>=0; k--) {

            // compute row . prev = sum_n row_n conj(prev_n)
            qcomp prod = 0;
            for (size_t n=0; n<row.size(); n++)
                prod += row[n] * conj(matr[k][n]);
                            
            // subtract (proj row onto prev) = (prod * prev) from final row
            for (size_t n=0; n<row.size(); n++)
                matr[i][n] -= prod * matr[k][n];
        }
    
        // compute row magnitude 
        qreal mag = 0;
        for (size_t j=0; j<row.size(); j++)
            mag += pow(abs(matr[i][j]), 2);
        mag = sqrt(mag);
        
        // normalise row
        for (size_t j=0; j<row.size(); j++)
            matr[i][j] /= mag;
    }
    
    // ensure matrix is indeed unitary 
    QMatrix conjprod = matr * getConjugateTranspose(matr);
    QMatrix iden = getIdentityMatrix(1 << numQb);
    
    // generating big unitary matrices is hard; if we fail, default to identity
    if ( numQb >= 3 && !areEqual(conjprod, iden) ) {
        
        matr = getIdentityMatrix(1 << numQb);
        conjprod = matr;
    }
    DEMAND( areEqual(conjprod, iden) );
    
    // return the new orthonormal matrix
    return matr;
}

std::vector<QMatrix> getRandomKrausMap(int numQb, int numOps) {
    DEMAND( numOps >= 1 );
    DEMAND( numOps <= 4*numQb*numQb );

    // generate random unitaries
    std::vector<QMatrix> ops;
    for (int i=0; i<numOps; i++)
        ops.push_back(getRandomUnitary(numQb));
        
    // generate random weights
    qreal weights[numOps];
    for (int i=0; i<numOps; i++)
        weights[i] = getRandomReal(0, 1);
        
    // normalise random weights
    qreal weightSum = 0;
    for (int i=0; i<numOps; i++)
        weightSum += weights[i];
    for (int i=0; i<numOps; i++)
        weights[i] = sqrt(weights[i]/weightSum);
        
    // normalise ops
    for (int i=0; i<numOps; i++)
        ops[i] *= weights[i];
        
    // check what we produced was a valid Kraus map
    QMatrix iden = getIdentityMatrix(1 << numQb);
    QMatrix prodSum = getZeroMatrix(1 << numQb);
    for (int i=0; i<numOps; i++)
        prodSum += getConjugateTranspose(ops[i]) * ops[i];
    DEMAND( areEqual(prodSum, iden) );
        
    return ops;
}

std::vector<QVector> getRandomOrthonormalVectors(int numQb, int numStates) {
    DEMAND( numQb >= 1 );
    DEMAND( numStates >= 1);
    
    // set of orthonormal vectors
    std::vector<QVector> vecs;
    
    for (size_t n=0; n<numStates; n++) {
        
        QVector vec = getRandomStateVector(numQb);
        
        // orthogonalise by substracting projections of existing vectors
        for (int m=0; m<n; m++) {
            qcomp prod = vec * vecs[m];
            vec -= (prod * vecs[m]);
        }
        
        // renormalise
        vec = getNormalised(vec);
        
        // add to orthonormal set
        vecs.push_back(vec);
    }

    return vecs;
}

QMatrix getMixedDensityMatrix(std::vector<qreal> probs, std::vector<QVector> states) {
    DEMAND( probs.size() == states.size() );
    DEMAND( probs.size() >= 1 );
    
    QMatrix matr = getZeroMatrix(states[0].size());
    
    for (size_t i=0; i<probs.size(); i++)
        matr += probs[i] * getPureDensityMatrix(states[i]);
        
    return matr;
}

QVector getDFT(QVector in) {
    REQUIRE( in.size() > 0 );
    
    size_t dim = in.size();
    qreal ampFac = 1 / sqrt( dim );
    qreal phaseFac = 2 * M_PI / dim;
    
    QVector dftVec = QVector(dim);
    
    for (size_t x=0; x<dim; x++) {
        dftVec[x] = 0;
        for (long long int y=0; y<dim; y++)
            dftVec[x] += expI(phaseFac * x * y) * in[y];
        dftVec[x] *= ampFac;
    }
    return dftVec;
}

long long int getValueOfTargets(long long int ind, int* targs, int numTargs) {
    DEMAND( ind >= 0 );
    
    long long int val = 0;
    
    for (int t=0; t<numTargs; t++)
        val += ((ind >> targs[t]) & 1) * (1LL << t);
        
    return val;
}

long long int setBit(long long int num, int bitInd, int bitVal) {
    DEMAND( (bitVal == 0 || bitVal == 1) );
    DEMAND( num >= 0 );
    DEMAND( bitInd >= 0 );
    
    return (num & ~(1UL << bitInd)) | (bitVal << bitInd);
}

long long int getIndexOfTargetValues(long long int ref, int* targs, int numTargs, int targVal) {
    // ref state is the starting index, where the targets can be in any bit state;
    // on the bits of the non-target qubits matter 
    
    for (int t=0; t<numTargs; t++) {
        int bit = (targVal >> t) & 1;
        ref = setBit(ref, targs[t], bit);    
    }
    return ref;
}

QVector getDFT(QVector in, int* targs, int numTargs) {
    
    QVector out = QVector(in.size());
    long long int targDim = (1LL << numTargs);
    
    for (size_t j=0; j<in.size(); j++) {
        
        // |j> = |x> (x) |...>, but mixed (not separated)
        long long int x = getValueOfTargets(j, targs, numTargs);
        
        for (long long int y=0; y<targDim; y++) {
            
            // modifies sum_y |y> (x) ...
            long long int outInd = getIndexOfTargetValues(j, targs, numTargs, y);
            
            qcomp elem = (in[j] / sqrt(pow(2,numTargs))) * expI(2*M_PI * x * y / pow(2,numTargs));
            out[outInd] += elem;
        }
    }
    
    return out;
}

/* (do not generate doxygen doc)
 *
 * Overloads for applyReferenceOp, to conveniently specify all families of 
 * unitary operations on state-vectors.
 */
void applyReferenceOp(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    int numQubits = calcLog2(state.size());
    QMatrix fullOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    state = fullOp * state;
}
void applyReferenceOp(
    QVector &state, int* ctrls, int numCtrls, int targ1, int targ2, QMatrix op
) {
    int targs[2] = {targ1, targ2};
    applyReferenceOp(state, ctrls, numCtrls, targs, 2, op);
}
void applyReferenceOp(
    QVector &state, int* ctrls, int numCtrls, int target, QMatrix op
) {
    int targs[1] = {target};
    applyReferenceOp(state, ctrls, numCtrls, targs, 1, op);
}
void applyReferenceOp(
    QVector &state, int *targs, int numTargs, QMatrix op
) {
    applyReferenceOp(state, NULL, 0, targs, numTargs, op);
}
void applyReferenceOp(
    QVector &state, int ctrl, int targ, QMatrix op
) {
    int ctrls[1] = {ctrl};
    int targs[1] = {targ};
    applyReferenceOp(state, ctrls, 1, targs, 1, op);
}
void applyReferenceOp(
    QVector &state, int ctrl, int* targs, int numTargs, QMatrix op
) {
    int ctrls[1] = {ctrl};
    applyReferenceOp(state, ctrls, 1, targs, numTargs, op);
}
void applyReferenceOp(
    QVector &state, int ctrl, int targ1, int targ2, QMatrix op
) {
    int ctrls[1] = {ctrl};
    int targs[2] = {targ1, targ2};
    applyReferenceOp(state, ctrls, 1, targs, 2, op);
}
void applyReferenceOp(
    QVector &state, int targ, QMatrix op
) {
    int targs[1] = {targ};
    applyReferenceOp(state, NULL, 0, targs, 1, op);
}

/* (do not generate doxygen doc)
 *
 * Overloads for applyReferenceOp, to conveniently specify all families of 
 * unitary operations on state-vectors.
 */
void applyReferenceOp(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    int numQubits = calcLog2(state.size());
    QMatrix leftOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    QMatrix rightOp = getConjugateTranspose(leftOp);
    state = leftOp * state * rightOp;
}
void applyReferenceOp(
    QMatrix &state, int* ctrls, int numCtrls, int targ1, int targ2, QMatrix op
) {
    int targs[2] = {targ1, targ2};
    applyReferenceOp(state, ctrls, numCtrls, targs, 2, op);
}
void applyReferenceOp(
    QMatrix &state, int* ctrls, int numCtrls, int target, QMatrix op
) {
    int targs[1] = {target};
    applyReferenceOp(state, ctrls, numCtrls, targs, 1, op);
}
void applyReferenceOp(
    QMatrix &state, int *targs, int numTargs, QMatrix op
) {
    applyReferenceOp(state, NULL, 0, targs, numTargs, op);
}
void applyReferenceOp(
    QMatrix &state, int ctrl, int targ, QMatrix op
) {
    int ctrls[1] = {ctrl};
    int targs[1] = {targ};
    applyReferenceOp(state, ctrls, 1, targs, 1, op);
}
void applyReferenceOp(
    QMatrix &state, int ctrl, int* targs, int numTargs, QMatrix op
) {
    int ctrls[1] = {ctrl};
    applyReferenceOp(state, ctrls, 1, targs, numTargs, op);
}
void applyReferenceOp(
    QMatrix &state, int ctrl, int targ1, int targ2, QMatrix op
) {
    int ctrls[1] = {ctrl};
    int targs[2] = {targ1, targ2};
    applyReferenceOp(state, ctrls, 1, targs, 2, op);
}
void applyReferenceOp(
    QMatrix &state, int targ, QMatrix op
) {
    int targs[1] = {targ};
    applyReferenceOp(state, NULL, 0, targs, 1, op);
}

/* (do not generate doxygen doc)
 *
 * Overloads for applyReferenceMatrix, to simply left-multiply a matrix (possibly
 * with additional control qubits) onto a state.
 */
void applyReferenceMatrix(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    // for state-vectors, the op is always just left-multiplied
    applyReferenceOp(state, ctrls, numCtrls, targs, numTargs, op);
}
void applyReferenceMatrix(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    // for density matrices, op is left-multiplied only
    int numQubits = calcLog2(state.size());
    QMatrix leftOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    state = leftOp * state;
}

bool areEqual(Qureg qureg1, Qureg qureg2, qreal precision) {
    DEMAND( qureg1.isDensityMatrix == qureg2.isDensityMatrix );
    DEMAND( qureg1.numAmpsTotal == qureg2.numAmpsTotal );
        
    copyStateFromGPU(qureg1);
    copyStateFromGPU(qureg2);
    syncQuESTEnv(QUEST_ENV);
    
    // loop terminates when areEqual = 0
    int ampsAgree = 1;
    for (long long int i=0; ampsAgree && i<qureg1.numAmpsPerChunk; i++)
        ampsAgree = (
               absReal(qureg1.stateVec.real[i] - qureg2.stateVec.real[i]) < precision
            && absReal(qureg1.stateVec.imag[i] - qureg2.stateVec.imag[i]) < precision);
            
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAmpsAgree = ampsAgree;
#ifdef DISTRIBUTED_MODE
    MPI_Allreduce(&ampsAgree, &allAmpsAgree, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif

    return allAmpsAgree;
}
bool areEqual(Qureg qureg1, Qureg qureg2) {
    return areEqual(qureg1, qureg2, REAL_EPS);
}

bool areEqual(Qureg qureg, QVector vec, qreal precision) {
    DEMAND( !qureg.isDensityMatrix );
    DEMAND( (int) vec.size() == qureg.numAmpsTotal );
    
    copyStateFromGPU(qureg);
    syncQuESTEnv(QUEST_ENV);
    
    // the starting index in vec of this node's qureg partition.
    long long int startInd = qureg.chunkId * qureg.numAmpsPerChunk;
            
    int ampsAgree = 1;
    for (long long int i=0; i<qureg.numAmpsPerChunk; i++) {
        qreal realDif = absReal(qureg.stateVec.real[i] - real(vec[startInd+i]));
        qreal imagDif = absReal(qureg.stateVec.imag[i] - imag(vec[startInd+i]));

        if (realDif > precision || imagDif > precision) {
            ampsAgree = 0;
            
            // debug
            char buff[200];
            sprintf(buff, "Disagreement at %lld of (%s) + i(%s):\n\t%s + i(%s) VS %s + i(%s)\n",
                startInd+i,
                REAL_STRING_FORMAT, REAL_STRING_FORMAT, REAL_STRING_FORMAT, 
                REAL_STRING_FORMAT, REAL_STRING_FORMAT, REAL_STRING_FORMAT);
            printf(buff,
                realDif, imagDif,
                qureg.stateVec.real[i], qureg.stateVec.imag[i],
                real(vec[startInd+i]), imag(vec[startInd+i]));
            
            break;
        }
    }
            
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAmpsAgree = ampsAgree;
#ifdef DISTRIBUTED_MODE
    MPI_Allreduce(&ampsAgree, &allAmpsAgree, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif
    
    return allAmpsAgree;
}
bool areEqual(Qureg qureg, QVector vec) {
    return areEqual(qureg, vec, REAL_EPS);
}

bool areEqual(Qureg qureg, QMatrix matr, qreal precision) {
    DEMAND( qureg.isDensityMatrix );
    DEMAND( (long long int) (matr.size()*matr.size()) == qureg.numAmpsTotal );
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    syncQuESTEnv(QUEST_ENV);
    
    // the starting index in vec of this node's qureg partition.
    long long int startInd = qureg.chunkId * qureg.numAmpsPerChunk;
    long long int globalInd, row, col, i;
    int ampsAgree;
    
    // compare each of this node's amplitude to the corresponding matr sub-matrix
    for (i=0; i<qureg.numAmpsPerChunk; i++) {
        globalInd = startInd + i;
        row = globalInd % matr.size();
        col = globalInd / matr.size();
        qreal realDif = absReal(qureg.stateVec.real[i] - real(matr[row][col]));
        qreal imagDif = absReal(qureg.stateVec.imag[i] - imag(matr[row][col]));
        ampsAgree = (realDif < precision && imagDif < precision);
        
        // DEBUG
        if (!ampsAgree) {
            
            // debug
            char buff[200];
            sprintf(buff, "[msg from utilities.cpp] node %d has a disagreement at (global) index %lld of (%s) + i(%s)\n",
                qureg.chunkId, globalInd, REAL_STRING_FORMAT, REAL_STRING_FORMAT);
            printf(buff, realDif, imagDif);
        }

        // break loop as soon as amplitudes disagree
        if (!ampsAgree)
            break;
            
        /* TODO:
         * of the nodes which disagree, the lowest-rank should send its 
         * disagreeing (i, row, col, stateVec[i]) to rank 0 which should 
         * report it immediately (before the impending DEMAND failure)
         * using FAIL_CHECK, so users can determine nature of disagreement 
         * (e.g. numerical precision).
         * Note FAIL_CHECK accepts << like cout, e.g.
         * FAIL_CHECK( "Amp at (" << row << ", " << col ") disagreed" );
         */
    }
    
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAmpsAgree = ampsAgree;
#ifdef DISTRIBUTED_MODE
    MPI_Allreduce(&ampsAgree, &allAmpsAgree, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif
        
    return allAmpsAgree;
}
bool areEqual(Qureg qureg, QMatrix matr) {
    return areEqual(qureg, matr, REAL_EPS);
}

bool areEqual(QVector vec, qreal* reals, qreal* imags) {
    
    qreal dif;
    for (size_t i=0; i<vec.size(); i++) {
        dif = absReal(real(vec[i]) - reals[i]);
        if (dif > REAL_EPS)
            return false;
        dif = absReal(imag(vec[i]) - imags[i]);
        if (dif > REAL_EPS)
            return false;
    }
    return true;
}

bool areEqual(QVector vec, qreal* reals) {
    for (size_t i=0; i<vec.size(); i++) {
        DEMAND( imag(vec[i]) == 0. );
        
        qreal dif = abs(real(vec[i]) - reals[i]);
        if (dif > REAL_EPS)
            return false;
    }
    return true;
}

/* Copies QMatrix into a CompelxMAtrix struct */
#define macro_copyQMatrix(dest, src) { \
    for (size_t i=0; i<src.size(); i++) { \
        for (size_t j=0; j<src.size(); j++) { \
            dest.real[i][j] = real(src[i][j]); \
            dest.imag[i][j] = imag(src[i][j]); \
        } \
    } \
}
ComplexMatrix2 toComplexMatrix2(QMatrix qm) {
    DEMAND( qm.size() == 2 );
    ComplexMatrix2 cm;
    macro_copyQMatrix(cm, qm);
    return cm;
}
ComplexMatrix4 toComplexMatrix4(QMatrix qm) {
    DEMAND( qm.size() == 4 );
    ComplexMatrix4 cm;
    macro_copyQMatrix(cm, qm);
    return cm;
}
void toComplexMatrixN(QMatrix qm, ComplexMatrixN cm) {
    DEMAND( qm.size() == (1u<<cm.numQubits) );
    macro_copyQMatrix(cm, qm);
}

/** Copies ComplexMatrix structures into a QMatrix */
#define macro_copyComplexMatrix(dest, src) { \
    for (size_t i=0; i<dest.size(); i++) \
        for (size_t j=0; j<dest.size(); j++) \
            dest[i][j] = qcomp(src.real[i][j], src.imag[i][j]); \
}
QMatrix toQMatrix(ComplexMatrix2 src) {
    QMatrix dest = getZeroMatrix(2);
    macro_copyComplexMatrix(dest, src);
    return dest;
}
QMatrix toQMatrix(ComplexMatrix4 src) {
    QMatrix dest = getZeroMatrix(4);
    macro_copyComplexMatrix(dest, src);
    return dest;
}
QMatrix toQMatrix(ComplexMatrixN src) {
    DEMAND( src.real != NULL );
    DEMAND( src.imag != NULL );
    QMatrix dest = getZeroMatrix(1 << src.numQubits);
    macro_copyComplexMatrix(dest, src);
    return dest;
}

QMatrix toQMatrix(Complex alpha, Complex beta) {
    qcomp a = qcomp(alpha.real, alpha.imag);
    qcomp b = qcomp(beta.real, beta.imag);
    QMatrix matr{
        {a, -conj(b)},
        {b,  conj(a)}};
    return matr;
}

QMatrix toQMatrix(Qureg qureg) {
    DEMAND( qureg.isDensityMatrix );
#ifdef DISTRIBUTED_MODE
    DEMAND( qureg.numAmpsTotal < MPI_MAX_AMPS_IN_MSG );
#endif
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    syncQuESTEnv(QUEST_ENV);
    
    qreal* fullRe;
    qreal* fullIm;
    
    // in distributed mode, give every node the full state vector
#ifdef DISTRIBUTED_MODE
    fullRe = (qreal*) malloc(qureg.numAmpsTotal * sizeof *fullRe);
    fullIm = (qreal*) malloc(qureg.numAmpsTotal * sizeof *fullIm);
    MPI_Allgather(
        qureg.stateVec.real, qureg.numAmpsPerChunk, MPI_QuEST_REAL,
        fullRe, qureg.numAmpsPerChunk, MPI_QuEST_REAL, MPI_COMM_WORLD);
    MPI_Allgather(
        qureg.stateVec.imag, qureg.numAmpsPerChunk, MPI_QuEST_REAL,
        fullIm, qureg.numAmpsPerChunk, MPI_QuEST_REAL, MPI_COMM_WORLD);
#else
    fullRe = qureg.stateVec.real;
    fullIm = qureg.stateVec.imag;
#endif
        
    // copy full state vector into a QVector
    long long int dim = (1 << qureg.numQubitsRepresented);
    QMatrix matr = getZeroMatrix(dim);
    for (long long int n=0; n<qureg.numAmpsTotal; n++)
        matr[n%dim][n/dim] = qcomp(fullRe[n], fullIm[n]);
    
    // clean up if we malloc'd the distributed array
#ifdef DISTRIBUTED_MODE
    free(fullRe);
    free(fullIm);
#endif
    return matr;
}

QVector toQVector(Qureg qureg) {
    DEMAND( !qureg.isDensityMatrix );
#ifdef DISTRIBUTED_MODE
    DEMAND( qureg.numAmpsTotal < MPI_MAX_AMPS_IN_MSG );
#endif
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    syncQuESTEnv(QUEST_ENV);
    
    qreal* fullRe;
    qreal* fullIm;
    
    // in distributed mode, give every node the full state vector
#ifdef DISTRIBUTED_MODE
    fullRe = (qreal*) malloc(qureg.numAmpsTotal * sizeof *fullRe);
    fullIm = (qreal*) malloc(qureg.numAmpsTotal * sizeof *fullIm);
            
    MPI_Allgather(
        qureg.stateVec.real, qureg.numAmpsPerChunk, MPI_QuEST_REAL,
        fullRe, qureg.numAmpsPerChunk, MPI_QuEST_REAL, MPI_COMM_WORLD);
    MPI_Allgather(
        qureg.stateVec.imag, qureg.numAmpsPerChunk, MPI_QuEST_REAL,
        fullIm, qureg.numAmpsPerChunk, MPI_QuEST_REAL, MPI_COMM_WORLD);
#else
    fullRe = qureg.stateVec.real;
    fullIm = qureg.stateVec.imag;
#endif
    
    // copy full state vector into a QVector
    QVector vec = QVector(qureg.numAmpsTotal);
    for (long long int i=0; i<qureg.numAmpsTotal; i++)
        vec[i] = qcomp(fullRe[i], fullIm[i]);
            
    // clean up if we malloc'd distrib array
#ifdef DISTRIBUTED_MODE
    free(fullRe);
    free(fullIm);
#endif
    return vec;
}

QVector toQVector(DiagonalOp op) {
    long long int totalElems = (1LL << op.numQubits);
#ifdef DISTRIBUTED_MODE
    DEMAND( totalElems < MPI_MAX_AMPS_IN_MSG );
#endif
    
    qreal* fullRe;
    qreal* fullIm;
    
    // in distributed mode, give every node the full diagonal operator
#ifdef DISTRIBUTED_MODE
    fullRe = (qreal*) malloc(totalElems * sizeof *fullRe);
    fullIm = (qreal*) malloc(totalElems * sizeof *fullIm);
            
    MPI_Allgather(
        op.real, op.numElemsPerChunk, MPI_QuEST_REAL,
        fullRe, op.numElemsPerChunk, MPI_QuEST_REAL, MPI_COMM_WORLD);
    MPI_Allgather(
        op.imag, op.numElemsPerChunk, MPI_QuEST_REAL,
        fullIm, op.numElemsPerChunk, MPI_QuEST_REAL, MPI_COMM_WORLD);
#else
    fullRe = op.real;
    fullIm = op.imag;
#endif
    
    // copy full state vector into a QVector
    QVector vec = QVector(totalElems);
    for (long long int i=0; i<totalElems; i++)
        vec[i] = qcomp(fullRe[i], fullIm[i]);
            
    // clean up if we malloc'd distrib array
#ifdef DISTRIBUTED_MODE
    free(fullRe);
    free(fullIm);
#endif
    return vec;
}

QMatrix toQMatrix(DiagonalOp op) {
    QVector vec = toQVector(op);
    QMatrix mat = getZeroMatrix(1LL << op.numQubits);
    for (size_t i=0; i<mat.size(); i++)
        mat[i][i] = vec[i];
    return mat;
}

void toQureg(Qureg qureg, QVector vec) {
    DEMAND( !qureg.isDensityMatrix );
    DEMAND( qureg.numAmpsTotal == (long long int) vec.size() );
    
    syncQuESTEnv(QUEST_ENV);
    
    for (int i=0; i<qureg.numAmpsPerChunk; i++) {
        int ind = qureg.chunkId*qureg.numAmpsPerChunk + i;
        qureg.stateVec.real[i] = real(vec[ind]);
        qureg.stateVec.imag[i] = imag(vec[ind]);
    }
    copyStateToGPU(qureg);
}
void toQureg(Qureg qureg, QMatrix mat) {
    DEMAND( qureg.isDensityMatrix );
    DEMAND( (1 << qureg.numQubitsRepresented) == (long long int) mat.size() );
    
    syncQuESTEnv(QUEST_ENV);
    
    int len = (1 << qureg.numQubitsRepresented);
    for (int i=0; i<qureg.numAmpsPerChunk; i++) {
        int ind = qureg.chunkId*qureg.numAmpsPerChunk + i;
        qureg.stateVec.real[i] = real(mat[ind%len][ind/len]);
        qureg.stateVec.imag[i] = imag(mat[ind%len][ind/len]);
    }
    copyStateToGPU(qureg);
}

void setRandomPauliSum(qreal* coeffs, pauliOpType* codes, int numQubits, int numTerms) {
    int i=0;
    for (int n=0; n<numTerms; n++) {
        coeffs[n] = getRandomReal(-5, 5);
        for (int q=0; q<numQubits; q++)
            codes[i++] = (pauliOpType) getRandomInt(0,4);
    }
}
void setRandomPauliSum(PauliHamil hamil) {
    setRandomPauliSum(hamil.termCoeffs, hamil.pauliCodes, hamil.numQubits, hamil.numSumTerms);
}

void setRandomDiagPauliHamil(PauliHamil hamil) {
    int i=0;
    for (int n=0; n<hamil.numSumTerms; n++) {
        hamil.termCoeffs[n] = getRandomReal(-5, 5);
        for (int q=0; q<hamil.numQubits; q++)
            if (getRandomReal(-1,1) > 0)
                hamil.pauliCodes[i++] = PAULI_Z;
            else
                hamil.pauliCodes[i++] = PAULI_I;
    }
}

QMatrix toQMatrix(qreal* coeffs, pauliOpType* paulis, int numQubits, int numTerms) {
    
    // produce a numTargs-big matrix 'pauliSum' by pauli-matrix tensoring and summing
    QMatrix iMatr{{1,0},{0,1}};
    QMatrix xMatr{{0,1},{1,0}};
    QMatrix yMatr{{0,-qcomp(0,1)},{qcomp(0,1),0}};
    QMatrix zMatr{{1,0},{0,-1}};
    QMatrix pauliSum = getZeroMatrix(1<<numQubits);
    
    for (int t=0; t<numTerms; t++) {
        QMatrix pauliProd = QMatrix{{1}};
        
        for (int q=0; q<numQubits; q++) {
            int i = q + t*numQubits;
            
            QMatrix fac;
            pauliOpType code = paulis[i];
            if (code == PAULI_I) fac = iMatr;
            if (code == PAULI_X) fac = xMatr;
            if (code == PAULI_Y) fac = yMatr;
            if (code == PAULI_Z) fac = zMatr;
            pauliProd = getKroneckerProduct(fac, pauliProd);
        }
        pauliSum += coeffs[t] * pauliProd;
    }
    
    // a now 2^numQubits by 2^numQubits Hermitian matrix
    return pauliSum;
}
QMatrix toQMatrix(PauliHamil hamil) {
    return toQMatrix(hamil.termCoeffs, hamil.pauliCodes, hamil.numQubits, hamil.numSumTerms);
}

long long int getTwosComplement(long long int decimal, int numBits) {
    DEMAND( decimal >= 0 );
    DEMAND( numBits >= 2 );
    DEMAND( decimal < (1LL << numBits) );
    
    long long int maxMag = 1LL << (numBits-1);
    if (decimal >= maxMag)
        return -maxMag + (decimal - maxMag);
    else
        return decimal;
}

long long int getUnsigned(long long int twosComp, int numBits) {
    DEMAND( numBits >= 2 );
    DEMAND( twosComp < (1LL << (numBits-1)) );
    DEMAND( twosComp >= - (1LL << (numBits-1)) );
    
    if (twosComp >= 0)
        return twosComp;
    else
        return (1<<numBits) + twosComp;
}

QMatrix toDiagonalQMatrix(QVector vec) {
    QMatrix mat = getZeroMatrix(vec.size());
    for (size_t i=0; i<vec.size(); i++)
        mat[i][i] = vec[i];
    return mat;
}

void setDiagMatrixOverrides(QMatrix &matr, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
    DEMAND( (encoding == UNSIGNED || encoding == TWOS_COMPLEMENT) );
    DEMAND( numRegs > 0 );
    DEMAND( numOverrides >= 0 );
    
    int totalQb = 0;
    for (int r=0; r<numRegs; r++) {
        DEMAND( numQubitsPerReg[r] > 0 );
        totalQb += numQubitsPerReg[r];
    }
    DEMAND( matr.size() == (1 << totalQb) );
    
    // record whether a diagonal index has been already overriden
    int hasBeenOverriden[1 << totalQb];
    for (int i=0; i<(1 << totalQb); i++)
        hasBeenOverriden[i] = 0;
    
    int flatInd = 0;
    for (int v=0; v<numOverrides; v++) {
        int matrInd = 0;
        int numQubitsLeft = 0;
        
        for (int r=0; r<numRegs; r++) {
            
            if (encoding == UNSIGNED)
                matrInd += overrideInds[flatInd] * (1 << numQubitsLeft);
            else if (encoding == TWOS_COMPLEMENT)
                matrInd += getUnsigned(overrideInds[flatInd], numQubitsPerReg[r]) * (1 << numQubitsLeft);
                
            numQubitsLeft += numQubitsPerReg[r];
            flatInd += 1;
        }
        
        if (!hasBeenOverriden[matrInd]) {
            matr[matrInd][matrInd] = expI(overridePhases[v]);
            hasBeenOverriden[matrInd] = 1;
        }
    }
}

static int fn_unique_suffix_id = 0;

void setUniqueFilename(char* outFn, char* prefix) {
    sprintf(outFn, "%s_%d.txt", prefix, fn_unique_suffix_id++);
}

void writeToFileSynch(char* fn, const string& contents) {
    
    // master node writes
    if (QUEST_ENV.rank == 0) {
        FILE* file = fopen(fn, "w");
        fputs(contents.c_str(), file);
        fclose(file);
    }
    
    // other nodes wait
    syncQuESTEnv(QUEST_ENV);
}

void deleteFilesWithPrefixSynch(char* prefix) {
    
    // master node deletes all files
    if (QUEST_ENV.rank == 0) {
        char cmd[200];
        sprintf(cmd, "exec rm %s*", prefix);
        system(cmd);
    }
    
    // other nodes wait 
    syncQuESTEnv(QUEST_ENV);
}

class SubListGenerator : public Catch::Generators::IGenerator<int*> {
    int* list;
    int* sublist;
    int len;
    int sublen;
    vector<bool> featured;
private:
    void createSublist() {
        
        // sublist to send to the user
        sublist = (int*) malloc(sublen * sizeof *sublist);
        
        // indicates which list members are currently in sublist
        featured = vector<bool>(len);
        fill(featured.end() - sublen, featured.end(), true);   
    }
    
    void prepareSublist() {
        
        // choose the next combination
        int j=0;
        for (int i=0; i<len; i++)
            if (featured[i])
                sublist[j++] = list[i];
                
        // prepare for permuting
        std::sort(sublist, sublist+sublen);
    }
public:
    SubListGenerator(int* elems, int numElems, int numSamps) {
        
        DEMAND( numSamps <= numElems );
                
        // make a record of all elements
        len = numElems;
        list = (int*) malloc(len * sizeof *list);
        for (int i=0; i<len; i++)
            list[i] = elems[i];
        
        // prepare sublist
        sublen = numSamps;
        createSublist();
        prepareSublist();  
    }
    
    SubListGenerator(
        Catch::Generators::GeneratorWrapper<int>&& gen, 
        int numSamps, const int* exclude, int numExclude
    ) {    
        // extract all generator elems
        vector<int> elems = vector<int>();
        do { elems.push_back(gen.get()); } while (gen.next());
        
        // make (int*) of non-excluded elems
        len = 0;
        list = (int*) malloc(elems.size() * sizeof *list);
        for (size_t i=0; i<elems.size(); i++) {
            int elem = elems[i];
            bool present = false;
            for (int j=0; j<numExclude; j++)
                if (elem == exclude[j]) {
                    present = true;
                    break;
                }
            if (!present)
                list[len++] = elem;
        }
        
        DEMAND( numSamps <= len );

        // prepare sublist
        sublen = numSamps;
        createSublist();
        prepareSublist();
    }
    
    int* const& get() const override {
        return sublist;
    }
    
    bool next() override {
        
        // offer next permutation of the current combination
        if (next_permutation(sublist, sublist+sublen))
            return true;

        // else generate the next combination
        if (next_permutation(featured.begin(), featured.end())) {
            prepareSublist();
            return true;
        }
        
        return false;
    }
    
    ~SubListGenerator() {
        free(list);
        free(sublist);
    }
};
Catch::Generators::GeneratorWrapper<int*> sublists(
    int* list, int len, int sublen
) {    
    return Catch::Generators::GeneratorWrapper<int*>(
        std::unique_ptr<Catch::Generators::IGenerator<int*>>(
            new SubListGenerator(list, len, sublen)));
}
Catch::Generators::GeneratorWrapper<int*> sublists(
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps, const int* exclude, int numExclude
) {    
    return Catch::Generators::GeneratorWrapper<int*>(
        std::unique_ptr<Catch::Generators::IGenerator<int*>>(
            new SubListGenerator(std::move(gen), numSamps, exclude, numExclude)));
}
Catch::Generators::GeneratorWrapper<int*> sublists(
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps, int excluded
) {
    int exclude[] = {excluded};  
    return Catch::Generators::GeneratorWrapper<int*>(
        std::unique_ptr<Catch::Generators::IGenerator<int*>>(
            new SubListGenerator(std::move(gen), numSamps, exclude, 1)));
}
Catch::Generators::GeneratorWrapper<int*> sublists(
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps
) {
    int exclude[] = {};  
    return Catch::Generators::GeneratorWrapper<int*>(
        std::unique_ptr<Catch::Generators::IGenerator<int*>>(
            new SubListGenerator(std::move(gen), numSamps, exclude, 0)));
}

template <typename T>
class SequenceGenerator : public Catch::Generators::IGenerator<T*> {
    T* digits;
    int len;
    T maxDigit;
    int ind;
    int seqLen;
public:
    SequenceGenerator(T maxDigit_, int numDigits) {
        ind = 0;
        len = numDigits;
        maxDigit = maxDigit_;
        seqLen = (int) pow(1 + (int) maxDigit, len);
        digits = (T*) malloc(numDigits * sizeof *digits);
        for (int i=0; i<numDigits; i++)
            digits[i] = (T) 0;
        ind++;
    }

    T* const& get() const override {
        return digits;
    }
    
    bool next() override {
        bool isNext = (ind++) < seqLen;
        if (isNext) {
            int i=0;
            while (digits[i] == maxDigit)
                digits[i++] = (T) 0;
            digits[i] = (T) ((int) digits[i] + 1);
        }
        return isNext;
    }
    
    ~SequenceGenerator() {
        free(digits);
    }
};
Catch::Generators::GeneratorWrapper<int*> bitsets(int numBits) {    
    return Catch::Generators::GeneratorWrapper<int*>(
        std::unique_ptr<Catch::Generators::IGenerator<int*>>(
            new SequenceGenerator<int>(1, numBits)));
}
Catch::Generators::GeneratorWrapper<int*> sequences(int base, int numDigits) {    
    return Catch::Generators::GeneratorWrapper<int*>(
        std::unique_ptr<Catch::Generators::IGenerator<int*>>(
            new SequenceGenerator<int>(base-1, numDigits)));
}
Catch::Generators::GeneratorWrapper<pauliOpType*> pauliseqs(int numPaulis) {    
    return Catch::Generators::GeneratorWrapper<pauliOpType*>(
        std::unique_ptr<Catch::Generators::IGenerator<pauliOpType*>>(
            new SequenceGenerator<pauliOpType>(PAULI_Z, numPaulis)));
}
