/** @file
 * Utility functions used by the ported v3 tests of QuEST's deprecated v3 API.
 *
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#define INCLUDE_DEPRECATED_FUNCTIONS 1
#define DISABLE_DEPRECATION_WARNINGS 1
#include "quest/include/quest.h"

#include "test_utilities.hpp"

#include <random>
#include <vector>
#include <algorithm>
#include <bitset>

#if COMPILE_MPI 

    #include <mpi.h>

    #if (FLOAT_PRECISION == 1)
        #define MPI_QCOMP MPI_CXX_FLOAT_COMPLEX
    #elif (FLOAT_PRECISION == 2)
        #define MPI_QCOMP MPI_CXX_DOUBLE_COMPLEX
    #elif (FLOAT_PRECISION == 4) && defined(MPI_CXX_LONG_DOUBLE_COMPLEX)
        #define MPI_QCOMP MPI_CXX_LONG_DOUBLE_COMPLEX
    #else
        #define MPI_QCOMP MPI_C_LONG_DOUBLE_COMPLEX
    #endif

    #ifdef MPI_MAX_AMPS_IN_MSG
    #undef MPI_MAX_AMPS_IN_MSG
    #endif
    #define MPI_MAX_AMPS_IN_MSG (1 << 30)

#endif

using std::vector;



/*
 * resolve reprecated absReal()
 */

#ifdef absReal
#undef absReal
#endif

// not sure where these will go! Maybe into QuEST v itself
qreal absReal(qreal x) { 
    return abs(x); // TODO: make precision agnostic
}
qreal absComp(qcomp x) {
  return abs(x); // TODO: make precision agnostic
}



/* RNG used for generating random test states,
 * independently used from C's rand() which is
 * used to generate random test data (e.g. operators)
 */

static std::mt19937 randomGenerator;



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

void setRandomTestStateSeeds() {

    // must seed both rand() and this C++ generator (used for shuffling)
    
    // obtain a seed from hardware
    std::random_device cspnrg;
    unsigned seed = cspnrg();
    
    // broadcast to ensure node consensus
#if COMPILE_MPI
    int sendRank = 0;
    MPI_Bcast(&seed, 1, MPI_UNSIGNED, sendRank, MPI_COMM_WORLD);
#endif

    // initilise both C (rand()) and C++ (randomGenerator) RNGs
    srand(seed);
    randomGenerator.seed(seed);
}

void assertQuregAndRefInDebugState(Qureg qureg, QVector ref) {
    DEMAND( qureg.isDensityMatrix == 0 );
    DEMAND( qureg.numAmps == (long long int) ref.size() );

    // assert ref is in the debug state (else initDebugState failed)
    for (size_t i=0; i<ref.size(); i++) {
        qcomp val = qcomp(.2*i, .2*i+.1);
        DEMAND( abs(ref[i] - val) < REAL_EPS );
    }

    // check qureg and ref agree
    DEMAND( areEqual(qureg, ref) );
}

void assertQuregAndRefInDebugState(Qureg qureg, QMatrix ref) {
    DEMAND( qureg.isDensityMatrix == 1 );
    DEMAND( (1LL << qureg.numQubits) == (long long int) ref.size() );

    // assert ref is in the (column-wise) debug state (else initDebugState failed)
    size_t i = 0;
    for (size_t c=0; c<ref.size(); c++) {
        for (size_t r=0; r<ref.size(); r++) {
            qcomp val = qcomp(.2*i, .2*i+.1);
            DEMAND( abs(ref[r][c] - val) < REAL_EPS );
            i++;
        }
    }

    // check qureg and ref agree
    DEMAND( areEqual(qureg, ref) );
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

QMatrix getTranspose(QMatrix a) {
    QMatrix b = a;
    for (size_t r=0; r<a.size(); r++)
        for (size_t c=0; c<a.size(); c++)
            b[r][c] = a[c][r];
    return b;
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
    vector<int> ctrlsCopy(ctrls, ctrls+numCtrls);
    vector<int> targsCopy(targs, targs+numTargs);
    
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
            if (a == 0) a = REAL_EPS; // prevent rand()=0 creation of NaN
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
    
    return vec;
}

QVector getNormalised(QVector vec) {

    // compute the vec norm via Kahan summation to suppress numerical error
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

vector<qreal> getRandomProbabilities(int numProbs) {
    
    // generate random unnormalised scalars
    vector<qreal> probs;
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
    vector<qreal> probs = getRandomProbabilities(dim);
    
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
    return (int) round(getRandomReal(min, max-1));
}

QMatrix getOrthonormalisedRows(QMatrix matr) {

    // perform the Gram-Schmidt process, processing each row of matr in-turn
    for (size_t i=0; i<matr.size(); i++) {
        QVector row = matr[i];
        
        // compute new orthogonal row by subtracting proj row onto prevs
        for (int k=i-1; k>=0; k--) {

            // compute inner_product(row, prev) = row . conj(prev)
            qcomp prod = row * matr[k];
                            
            // subtract (proj row onto prev) = (prod * prev) from final row
            matr[i] -= prod * matr[k];
        }
    
        // normalise the row 
        matr[i] = getNormalised(matr[i]);
    }
    
    // return the new orthonormal matrix
    return matr;
}

QMatrix getRandomDiagonalUnitary(int numQb) {
    DEMAND( numQb >= 1 );

    QMatrix matr = getZeroMatrix(1 << numQb);
    for (size_t i=0; i<matr.size(); i++)
        matr[i][i] = expI(getRandomReal(0,4*M_PI));

    return matr;
}

QMatrix getRandomUnitary(int numQb) {
    DEMAND( numQb >= 1 );

    // create Z ~ random complex matrix (distribution not too important)
    size_t dim = 1 << numQb;
    QMatrix matrZ = getRandomQMatrix(dim);
    QMatrix matrZT = getTranspose(matrZ);

    // create Z = Q R (via QR decomposition) ...
    QMatrix matrQT = getOrthonormalisedRows(matrZ);
    QMatrix matrQ = getTranspose(matrQT);
    QMatrix matrR = getZeroMatrix(dim);

    // ... where R_rc = (columm c of Z) . (column r of Q) = (row c of ZT) . (row r of QT)
    for (size_t r=0; r<dim; r++)
        for (size_t c=r; c<dim; c++)
            matrR[r][c] = matrZT[c] * matrQT[r];

    // create D = normalised diagonal of R
    QMatrix matrD = getZeroMatrix(dim);
    for (size_t i=0; i<dim; i++)
        matrD[i][i] = matrR[i][i] / abs(matrR[i][i]);

    // create U = Q D
    QMatrix matrU = matrQ * matrD;

    // in the rare scenario the result is not sufficiently precisely unitary,
    // replace it with a trivially unitary diagonal matrix
    QMatrix daggerProd = matrU * getConjugateTranspose(matrU);
    QMatrix iden = getIdentityMatrix(dim);
    if( ! areEqual(daggerProd, iden) )
        matrU = getRandomDiagonalUnitary(numQb);

    return matrU;
}

vector<QMatrix> getRandomKrausMap(int numQb, int numOps) {
    DEMAND( numOps >= 1 );
    DEMAND( numOps <= 4*numQb*numQb );

    // generate random unitaries
    vector<QMatrix> ops;
    for (int i=0; i<numOps; i++)
        ops.push_back(getRandomUnitary(numQb));
        
    // generate random weights
    vector<qreal> weights(numOps);
    for (int i=0; i<numOps; i++)
        weights[i] = getRandomReal(0, 1);
        
    // normalise random weights
    qreal weightSum = 0;
    for (int i=0; i<numOps; i++)
        weightSum += weights[i];
    for (int i=0; i<numOps; i++)
        weights[i] = sqrt((qreal) weights[i]/weightSum);
        
    // normalise ops
    for (int i=0; i<numOps; i++)
        ops[i] *= weights[i];
        
    // check what we produced was a valid Kraus map
    QMatrix iden = getIdentityMatrix(1 << numQb);
    QMatrix prodSum = getZeroMatrix(1 << numQb);
    for (int i=0; i<numOps; i++)
        prodSum += getConjugateTranspose(ops[i]) * ops[i];

    // in the rare scenario it is insufficiently numerically precise,
    // replace the map with trivially precise diagonals
    if( ! areEqual(prodSum, iden) )
        for (int i=0; i<numOps; i++)
            ops[i] = weights[i] * getRandomDiagonalUnitary(numQb);

    return ops;
}

vector<QVector> getRandomOrthonormalVectors(int numQb, int numStates) {
    DEMAND( numQb >= 1 );
    DEMAND( numStates >= 1);
    
    // set of orthonormal vectors
    vector<QVector> vecs;
    
    for (int n=0; n<numStates; n++) {
        
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

QMatrix getMixedDensityMatrix(vector<qreal> probs, vector<QVector> states) {
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
        for (size_t y=0; y<dim; y++)
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

long long int getIndexOfTargetValues(long long int ref, int* targs, int numTargs, long long int targVal) {
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
    long long int inDim = (long long int) in.size();
    long long int targDim = (1LL << numTargs);
    
    for (long long int j=0; j<inDim; j++) {
        
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
    QVector &state, int *targs, int numTargs, QMatrix op
) {
    // for state-vectors, the op is always just left-multiplied
    applyReferenceOp(state, targs, numTargs, op);
}
void applyReferenceMatrix(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    // for density matrices, op is left-multiplied only
    int numQubits = calcLog2(state.size());
    QMatrix leftOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    state = leftOp * state;
}
void applyReferenceMatrix(
    QMatrix &state, int *targs, int numTargs, QMatrix op
) {
    applyReferenceMatrix(state, NULL, 0, targs, numTargs, op);
}

bool areEqual(Qureg qureg1, Qureg qureg2, qreal precision) {
    DEMAND( qureg1.isDensityMatrix == qureg2.isDensityMatrix );
    DEMAND( qureg1.numAmps == qureg2.numAmps );
        
    copyStateFromGPU(qureg1);
    copyStateFromGPU(qureg2);
    syncQuESTEnv();
    
    // loop terminates when areEqual = 0
    int ampsAgree = 1;
    for (long long int i=0; ampsAgree && i<qureg1.numAmpsPerNode; i++)
        ampsAgree = absComp(qureg1.cpuAmps[i] - qureg2.cpuAmps[i]) < precision;
            
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAmpsAgree = ampsAgree;
#if COMPILE_MPI
    MPI_Allreduce(&ampsAgree, &allAmpsAgree, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif

    return allAmpsAgree;
}
bool areEqual(Qureg qureg1, Qureg qureg2) {
    return areEqual(qureg1, qureg2, REAL_EPS);
}

bool areEqual(Qureg qureg, QVector vec, qreal precision) {
    DEMAND( !qureg.isDensityMatrix );
    DEMAND( (int) vec.size() == qureg.numAmps );
    
    copyStateFromGPU(qureg);
    syncQuESTEnv();
    
    // the starting index in vec of this node's qureg partition.
    long long int startInd = qureg.rank * qureg.numAmpsPerNode;
            
    int ampsAgree = 1;
    for (long long int i=0; i<qureg.numAmpsPerNode; i++) {
        qcomp dif = (qureg.cpuAmps[i] - vec[startInd+i]);

        if (absComp(dif) > precision) {
            ampsAgree = 0;

            // debug
            char buff[200];
            snprintf(buff, 200, "Disagreement at %lld of (%s) + i(%s):\n\t%s + i(%s) VS %s + i(%s)\n",
                startInd+i,
                QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER, 
                QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER);
            printf(buff,
                real(dif), imag(dif),
                real(qureg.cpuAmps[i]), imag(qureg.cpuAmps[i]),
                real(vec[startInd+i]), imag(vec[startInd+i]));
            
            break;
        }
    }
            
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAmpsAgree = ampsAgree;
#if COMPILE_MPI
    MPI_Allreduce(&ampsAgree, &allAmpsAgree, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif
    
    return allAmpsAgree;
}
bool areEqual(Qureg qureg, QVector vec) {
    return areEqual(qureg, vec, REAL_EPS);
}

bool areEqual(Qureg qureg, QMatrix matr, qreal precision) {
    DEMAND( qureg.isDensityMatrix );
    DEMAND( (long long int) (matr.size()*matr.size()) == qureg.numAmps );
    
    // ensure local qureg amps is up to date
    copyStateFromGPU(qureg);
    syncQuESTEnv();
    
    // the starting index in vec of this node's qureg partition.
    long long int startInd = qureg.rank * qureg.numAmpsPerNode;
    long long int globalInd, row, col, i;
    int ampsAgree;
    
    // compare each of this node's amplitude to the corresponding matr sub-matrix
    for (i=0; i<qureg.numAmpsPerNode; i++) {
        globalInd = startInd + i;
        row = globalInd % matr.size();
        col = globalInd / matr.size();

        qreal realDif = absReal(real(qureg.cpuAmps[i]) - real(matr[row][col]));
        qreal imagDif = absReal(imag(qureg.cpuAmps[i]) - imag(matr[row][col]));
        ampsAgree = (realDif < precision && imagDif < precision);
        
        // DEBUG
        if (!ampsAgree) {
            char buff[200];
            snprintf(buff, 200, "[msg from utilities.cpp] node %d has a disagreement at %lld of (%s) + i(%s):\n\t[qureg] %s + i(%s) VS [ref] %s + i(%s)\n",
                qureg.rank, startInd+i,
                QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER, 
                QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER, QREAL_FORMAT_SPECIFIER);
            printf(buff,
                realDif, imagDif,
                real(qureg.cpuAmps[i]), imag(qureg.cpuAmps[i]),
                real(matr[row][col]),   imag(matr[row][col]));
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
#if COMPILE_MPI
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
#define macro_copyQMatrixToDeprecatedComplexMatrix(dest, src) { \
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
    macro_copyQMatrixToDeprecatedComplexMatrix(cm, qm);
    return cm;
}
ComplexMatrix4 toComplexMatrix4(QMatrix qm) {
    DEMAND( qm.size() == 4 );
    ComplexMatrix4 cm;
    macro_copyQMatrixToDeprecatedComplexMatrix(cm, qm);
    return cm;
}

/** Copies ComplexMatrix structures into a QMatrix */
#define macro_copyComplexMatrix(dest, src, dim) \
    for (size_t i=0; i<dim; i++) \
        for (size_t j=0; j<dim; j++) \
            dest[i][j] = src[i][j];

void toComplexMatrixN(QMatrix qm, ComplexMatrixN cm) {
    DEMAND( qm.size() == (1u<<cm.numQubits) );
    macro_copyComplexMatrix(cm.cpuElems, qm, qm.size());
    syncCompMatr(cm);
}

QMatrix toQMatrix(CompMatr1 src) {
    QMatrix dest = getZeroMatrix(2);
    macro_copyComplexMatrix(dest, src.elems, dest.size());
    return dest;
}
QMatrix toQMatrix(CompMatr2 src) {
    QMatrix dest = getZeroMatrix(4);
    macro_copyComplexMatrix(dest, src.elems, dest.size());
    return dest;
}
QMatrix toQMatrix(CompMatr src) {
    QMatrix dest = getZeroMatrix(1 << src.numQubits);
    macro_copyComplexMatrix(dest, src.cpuElems, dest.size());
    return dest;
}

QMatrix toQMatrix(Qureg qureg) {
    DEMAND( qureg.isDensityMatrix );
#if COMPILE_MPI
    DEMAND( qureg.numAmps < MPI_MAX_AMPS_IN_MSG );
#endif
    
    // ensure local qureg amps are up to date
    copyStateFromGPU(qureg);
    syncQuESTEnv();

    // collect all amps between all nodes
    qcomp* allAmps = qureg.cpuAmps;
    
    // in distributed mode, give every node the full state vector
#if COMPILE_MPI
    if (qureg.isDistributed) {
        allAmps = (qcomp*) malloc(qureg.numAmps * sizeof *allAmps);
        MPI_Allgather(
            qureg.cpuAmps, qureg.numAmpsPerNode, MPI_QCOMP,
            allAmps, qureg.numAmpsPerNode, MPI_QCOMP, MPI_COMM_WORLD);
    }
#endif
        
    // copy full state vector into a QVector
    long long int dim = (1LL << qureg.numQubits);
    QMatrix matr = getZeroMatrix(dim);
    for (long long int n=0; n<qureg.numAmps; n++)
        matr[n%dim][n/dim] = allAmps[n];
    
    // clean up if we malloc'd the distributed array
    if (qureg.isDistributed)
        free(allAmps);
    return matr;
}

QVector toQVector(Qureg qureg) {
    DEMAND( !qureg.isDensityMatrix );
#if COMPILE_MPI
    DEMAND( qureg.numAmps < MPI_MAX_AMPS_IN_MSG );
#endif
    
    // ensure local qureg amps are up to date
    copyStateFromGPU(qureg);
    syncQuESTEnv();
    
    qcomp* allAmps = qureg.cpuAmps;
    
    // in distributed mode, give every node the full state vector
#if COMPILE_MPI
    if (qureg.isDistributed) {
        allAmps = (qcomp*) malloc(qureg.numAmps * sizeof *allAmps);

        MPI_Allgather(
            qureg.cpuAmps, qureg.numAmpsPerNode, MPI_QCOMP,
            allAmps, qureg.numAmpsPerNode, MPI_QCOMP, MPI_COMM_WORLD);
    }
#endif
    
    // copy full state vector into a QVector
    QVector vec = QVector(qureg.numAmps);
    for (long long int i=0; i<qureg.numAmps; i++)
        vec[i] = allAmps[i];
            
    // clean up if we malloc'd distrib array
    if (qureg.isDistributed)
        free(allAmps);

    return vec;
}

QVector toQVector(DiagMatr matr) {

    return vector<qcomp>(matr.cpuElems, matr.cpuElems + matr.numElems);
}

QVector toQVector(FullStateDiagMatr matr) {

#if COMPILE_MPI
    DEMAND( matr.numElems < MPI_MAX_AMPS_IN_MSG );
#endif

    vector<qcomp> vec(matr.numElems);

    // in distributed mode, give every node the full diagonal operator
    if (matr.isDistributed) {
        #if COMPILE_MPI
            MPI_Allgather(
                matr.cpuElems, matr.numElemsPerNode, MPI_QCOMP,
                vec.data(),    matr.numElemsPerNode, MPI_QCOMP, MPI_COMM_WORLD);
        #endif
    } else {
        vec.assign(matr.cpuElems, matr.cpuElems + matr.numElems);
    }

    return vec;
}

QMatrix toQMatrix(FullStateDiagMatr in) {
    QVector vec = toQVector(in);
    QMatrix mat = getZeroMatrix(in.numElems);
    for (size_t i=0; i<mat.size(); i++)
        mat[i][i] = vec[i];
    return mat;
}

QMatrix toQMatrix(DiagMatr in) {
    QMatrix mat = getZeroMatrix(in.numElems);
    for (size_t i=0; i<mat.size(); i++)
        mat[i][i] = in.cpuElems[i];
    return mat;
}

void toQureg(Qureg qureg, QVector vec) {
    DEMAND( !qureg.isDensityMatrix );
    DEMAND( qureg.numAmps == (long long int) vec.size() );
    
    syncQuESTEnv();
    
    for (int i=0; i<qureg.numAmpsPerNode; i++) {
        int ind = qureg.rank*qureg.numAmpsPerNode + i;
        qureg.cpuAmps[i] = vec[ind];
    }
    copyStateToGPU(qureg);
}
void toQureg(Qureg qureg, QMatrix mat) {
    DEMAND( qureg.isDensityMatrix );
    DEMAND( (1LL << qureg.numQubits) == (long long int) mat.size() );
    
    syncQuESTEnv();
    
    int len = (1 << qureg.numQubits);
    for (int i=0; i<qureg.numAmpsPerNode; i++) {
        int ind = qureg.rank*qureg.numAmpsPerNode + i;
        qureg.cpuAmps[i] = mat[ind%len][ind/len];
    }
    copyStateToGPU(qureg);
}

PauliStr getRandomPauliStr(int numQubits) {

    std::string paulis = "";
    for (int i=0; i<numQubits; i++)
        paulis += "IXYZ"[getRandomInt(0,4)];

    return getPauliStr(paulis);
}
PauliStr getRandomDiagPauliStr(int numQubits) {

    std::string paulis = "";
    for (int i=0; i<numQubits; i++)
        paulis += "IX"[getRandomInt(0,2)];

    return getPauliStr(paulis);
}

void setRandomPauliSum(qreal* coeffs, pauliOpType* codes, int numQubits, int numTerms) {
    int i=0;
    for (int n=0; n<numTerms; n++) {
        coeffs[n] = getRandomReal(-5, 5);
        for (int q=0; q<numQubits; q++)
            codes[i++] = (pauliOpType) getRandomInt(0,4);
    }
}

void setRandomPauliSum(PauliHamil hamil, int numQubits) {

    for (int n=0; n<hamil.numTerms; n++) {
        hamil.coeffs[n] = getRandomReal(-5, 5);
        hamil.strings[n] = getRandomPauliStr(numQubits);
    }
}

void setRandomDiagPauliHamil(PauliHamil hamil, int numQubits) {
    for (int n=0; n<hamil.numTerms; n++) {
        hamil.coeffs[n] = getRandomReal(-5, 5);
        hamil.strings[n] = getRandomDiagPauliStr(numQubits);
    }
}

void setRandomTargets(int* targs, int numTargs, int numQb) {
    DEMAND( numQb >= 1 );
    DEMAND( numTargs >= 1);
    DEMAND( numTargs <= numQb );

    // create an ordered list of all possible qubits
    vector<int> allQb(numQb);
    for (int q=0; q<numQb; q++)
        allQb[q] = q;

    // shuffle all qubits (must be consistent on each node)
    std::shuffle(&allQb[0], &allQb[numQb], randomGenerator);

    // select numTargs of all qubits
    for (int i=0; i<numTargs; i++)
        targs[i] = allQb[i];
}
void setRandomTargets(vector<int> &targs, int numQb) {

    setRandomTargets(targs.data(), targs.size(), numQb);
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

// QMatrix toQMatrix(PauliHamil hamil) {
//     return toQMatrix(hamil.termCoeffs, hamil.pauliCodes, hamil.numQubits, hamil.numSumTerms);
// }


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

// void setDiagMatrixOverrides(QMatrix &matr, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
//     DEMAND( (encoding == UNSIGNED || encoding == TWOS_COMPLEMENT) );
//     DEMAND( numRegs > 0 );
//     DEMAND( numOverrides >= 0 );
    
//     int totalQb = 0;
//     for (int r=0; r<numRegs; r++) {
//         DEMAND( numQubitsPerReg[r] > 0 );
//         totalQb += numQubitsPerReg[r];
//     }
//     DEMAND( matr.size() == (1 << totalQb) );
    
//     // record whether a diagonal index has been already overriden
//     vector<int> hasBeenOverriden(1 << totalQb);
//     for (int i=0; i<(1 << totalQb); i++)
//         hasBeenOverriden[i] = 0;
    
//     int flatInd = 0;
//     for (int v=0; v<numOverrides; v++) {
//         int matrInd = 0;
//         int numQubitsLeft = 0;
        
//         for (int r=0; r<numRegs; r++) {
            
//             if (encoding == UNSIGNED)
//                 matrInd += overrideInds[flatInd] * (1 << numQubitsLeft);
//             else if (encoding == TWOS_COMPLEMENT)
//                 matrInd += getUnsigned(overrideInds[flatInd], numQubitsPerReg[r]) * (1 << numQubitsLeft);
                
//             numQubitsLeft += numQubitsPerReg[r];
//             flatInd += 1;
//         }
        
//         if (!hasBeenOverriden[matrInd]) {
//             matr[matrInd][matrInd] = expI(overridePhases[v]);
//             hasBeenOverriden[matrInd] = 1;
//         }
//     }
// }

static int fn_unique_suffix_id = 0;

void setUniqueFilename(char* outFn, int maxlen, char* prefix) {
    snprintf(outFn, maxlen, "%s_%d.txt", prefix, fn_unique_suffix_id++);
}

void writeToFileSynch(char* fn, const string& contents) {
    
    // master node writes
    if (getQuESTEnv().rank == 0) {
        FILE* file = fopen(fn, "w");
        fputs(contents.c_str(), file);
        fclose(file);
    }
    
    // other nodes wait
    syncQuESTEnv();
}

void deleteFilesWithPrefixSynch(char* prefix) {
    
    // master node deletes all files
    if (getQuESTEnv().rank == 0) {
        char cmd[200];
        snprintf(cmd, 200, "exec rm %s*", prefix);
        system(cmd);
    }
    
    // other nodes wait 
    syncQuESTEnv();
}

/// @private
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
        if (std::next_permutation(sublist, sublist+sublen))
            return true;

        // else generate the next combination
        if (std::next_permutation(featured.begin(), featured.end())) {
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
            Catch::Detail::make_unique<SubListGenerator>(list, len, sublen));
}
Catch::Generators::GeneratorWrapper<int*> sublists(
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps, const int* exclude, int numExclude
) {    
    return Catch::Generators::GeneratorWrapper<int*>(
        Catch::Detail::make_unique<SubListGenerator>(std::move(gen), numSamps, exclude, numExclude));
}
Catch::Generators::GeneratorWrapper<int*> sublists(
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps, int excluded
) {
    int exclude[] = {excluded};  
    return Catch::Generators::GeneratorWrapper<int*>(
        Catch::Detail::make_unique<SubListGenerator>(std::move(gen), numSamps, exclude, 1));
}
Catch::Generators::GeneratorWrapper<int*> sublists(
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps
) {
    int exclude[] = {-1}; // non-empty to satisfy MSVC
    return Catch::Generators::GeneratorWrapper<int*>(
            Catch::Detail::make_unique<SubListGenerator>(std::move(gen), numSamps, exclude, 0));
}

/// @private
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
        Catch::Detail::make_unique<SequenceGenerator<int>>(1, numBits));
}
Catch::Generators::GeneratorWrapper<int*> sequences(int base, int numDigits) {    
    return Catch::Generators::GeneratorWrapper<int*>(
        Catch::Detail::make_unique<SequenceGenerator<int>>(base-1, numDigits));
}
Catch::Generators::GeneratorWrapper<pauliOpType*> pauliseqs(int numPaulis) {    
    return Catch::Generators::GeneratorWrapper<pauliOpType*>(
        Catch::Detail::make_unique<SequenceGenerator<pauliOpType>>(PAULI_Z, numPaulis));
}
