/** @file
 * Simple implementations of matrix operations used by QuEST_unit_tests
 *
 * @author Tyson Jones
 */

#include "QuEST.h"
#include "QuEST_test_utils.hpp"
#include "catch.hpp"

/** produces a dim-by-dim square complex matrix, initialised to zero 
 */
QMatrix getZeroMatrix(size_t dim) {
    REQUIRE( dim > 1 );
    QMatrix matr = QMatrix(dim);
    for (size_t i=0; i<dim; i++)
        matr[i].resize(dim);
    return matr;
}

/** produces a dim-by-dim identity matrix 
 */
QMatrix getIdentityMatrix(size_t dim) {
    REQUIRE( dim > 1 );
    QMatrix matr = getZeroMatrix(dim);
    for (size_t i=0; i<dim; i++)
        matr[i][i] = 1;
    return matr;
}

/** returns a (otimes) b, where a and b are square but possibly different-sized 
 */
QMatrix getKroneckerProduct(QMatrix a, QMatrix b) {
    QMatrix prod = getZeroMatrix(a.size() * b.size());
    for (size_t r=0; r<b.size(); r++)
        for (size_t c=0; c<b.size(); c++)
            for (size_t i=0; i<a.size(); i++)
                for (size_t j=0; j<a.size(); j++)
                    prod[r+b.size()*i][c+b.size()*j] = a[i][j] * b[r][c];
    return prod;
}

/** returns a square matrix, the product of a and b 
 */
QMatrix getMatrixProduct(QMatrix a, QMatrix b) {
    REQUIRE( a.size() == b.size() );
    QMatrix prod = getZeroMatrix(a.size());
    for (size_t r=0; r<a.size(); r++)
        for (size_t c=0; c<a.size(); c++)
            for (size_t k=0; k<a.size(); k++)
                prod[r][c] += a[r][k] * b[k][c];
    return prod;
}

/** returns the conjugate transpose of the complex square matrix a 
 */
QMatrix getConjugateTranspose(QMatrix a) {
    QMatrix b = a;
    for (size_t r=0; r<a.size(); r++)
        for (size_t c=0; c<a.size(); c++)
            b[r][c] = conj(a[c][r]);
    return b;
}

/** modifies dest by overwriting its submatrix (from top-left corner 
 * (r, c) to bottom-right corner (r+dest.size(), c+dest.size()) with the 
 * complete elements of sub 
 */
void setSubMatrix(QMatrix &dest, QMatrix sub, size_t r, size_t c) {
    REQUIRE( sub.size() + r <= dest.size() );
    REQUIRE( sub.size() + c <= dest.size() );
    for (size_t i=0; i<sub.size(); i++)
        for (size_t j=0; j<sub.size(); j++)
            dest[r+i][c+j] = sub[i][j];
}

/** returns the 2^numQb-by-2^numQb unitary matrix which swaps qubits qb1 and qb2.
 * If qb1==qb2, returns the identity matrix.
 */
QMatrix getSwapMatrix(int qb1, int qb2, int numQb) {
    REQUIRE( numQb > 1 );
    REQUIRE( (qb1 >= 0 && qb1 < numQb) );
    REQUIRE( (qb2 >= 0 && qb2 < numQb) );
    
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

/** iterates list1 (of length len1) and replaces element oldEl with newEl, which is 
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

/** takes a 2^numTargs-by-2^numTargs matrix op and a returns a 2^numQubits-by-2^numQubits
 * matrix where op is controlled on the given ctrls qubits. The union of {ctrls}
 * and {targs} must be unique, and every element must be 0 or positive. 
 * The passed {ctrls} and {targs} arrays are unmodified.
 * This funciton works by first swapping {targs} and {ctrls} (via swap unitaries) 
 * to be strictly increasing {0,1,...}, building controlled(op), tensoring it to 
 * the full Hilbert space, and then 'unswapping'. The returned matrix has form:
 * swap1 ... swapN . c(op) . swapN ... swap1
 */
QMatrix getFullOperatorMatrix(
    int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
) {        
    REQUIRE( numCtrls >= 0 );
    REQUIRE( numTargs >= 0 );
    REQUIRE( numQubits >= (numCtrls+numTargs) );
    REQUIRE( op.size() == (1 << numTargs) );
    
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
            swaps = getMatrixProduct(matr, swaps);
            unswaps = getMatrixProduct(unswaps, matr);
            
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
            swaps = getMatrixProduct(matr, swaps);
            unswaps = getMatrixProduct(unswaps, matr);
            
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
    fullOp = getMatrixProduct(fullOp, swaps);
    fullOp = getMatrixProduct(unswaps, fullOp);
    
    // restore {ctrls and targs}
    for (int i=0; i<numCtrls; i++)
        ctrls[i] = ctrlsCopy[i];
    for (int i=0; i<numTargs; i++)
        targs[i] = targsCopy[i];

    return fullOp;
}

/** returns the product of complex matrix m onto complex vector v 
 */
QVector getMatrixVectorProduct(QMatrix m, QVector v) {
    REQUIRE( m.size() == v.size() );
    QVector prod = QVector(v.size());
    for (size_t r=0; r<v.size(); r++)
        for (size_t c=0; c<v.size(); c++)
            prod[r] += m[r][c] * v[c];
    return prod;
}

/** overwrites state to be the result of applying the unitary matrix op (with the
 * specified control and target qubits) to statevector state, i.e. fullOp(op) * |state>
 */
void applyQUnitary(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
) {
    REQUIRE( state.size() == (1<<numQubits) );
    QMatrix fullOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    state = getMatrixVectorProduct(fullOp, state);
}
void applyQUnitary(
    QVector &state, int *targs, int numTargs, QMatrix op, int numQubits
) {
    applyQUnitary(state, NULL, 0, targs, numTargs, op, numQubits);
}
void applyQUnitary(
    QVector &state, int ctrl, int targ, QMatrix op, int numQubits
) {
    int ctrls[1] = {ctrl};
    int targs[1] = {targ};
    applyQUnitary(state, ctrls, 1, targs, 1, op, numQubits);
}
void applyQUnitary(
    QVector &state, int targ, QMatrix op, int numQubits
) {
    int targs[1] = {targ};
    applyQUnitary(state, NULL, 0, targs, 1, op, numQubits);
}

/** overwrites state to be the result of applying the unitary matrix op (with the
 * specified control and target qubits) to density matrix state, i.e. 
 *  fullOp(op) * state * fullOp(op)^dagger
 */
void applyQUnitary(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits
) {
    REQUIRE( state.size() == (1<<numQubits) );
    QMatrix leftOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    QMatrix rightOp = getConjugateTranspose(leftOp);
    state = getMatrixProduct(getMatrixProduct(leftOp, state), rightOp);
}
void applyQUnitary(
    QMatrix &state, int *targs, int numTargs, QMatrix op, int numQubits
) {
    applyQUnitary(state, NULL, 0, targs, numTargs, op, numQubits);
}
void applyQUnitary(
    QMatrix &state, int ctrl, int targ, QMatrix op, int numQubits
) {
    int ctrls[1] = {ctrl};
    int targs[1] = {targ};
    applyQUnitary(state, ctrls, 1, targs, 1, op, numQubits);
}
void applyQUnitary(
    QMatrix &state, int targ, QMatrix op, int numQubits
) {
    int targs[1] = {targ};
    applyQUnitary(state, NULL, 0, targs, 1, op, numQubits);
}



/** hardware-agnostic comparison of the given vector and pure qureg, to within 
 * the QuEST_PREC-specific REAL_EPS (defined in QuEST_precision) precision.
 * In GPU mode, this involves a GPU to CPU memory copy overhead.
 * In distributed mode, this involves a all-to-all single-int broadcast
 */
bool areEqual(QVector vec, Qureg qureg) {
    REQUIRE( !qureg.isDensityMatrix );
    REQUIRE( vec.size() == qureg.numAmpsTotal );
    
    copyStateFromGPU(qureg);
    
    // the starting index in vec of this node's qureg partition.
    long long int startInd = qureg.chunkId * qureg.numAmpsPerChunk;
    
    // loop terminates when areEqual = 0
    int areEqual = 1;
    for (long long int i=0; areEqual && i<qureg.numAmpsPerChunk; i++)
        areEqual = (
               absReal(qureg.stateVec.real[i] - real(vec[startInd+i])) < REAL_EPS
            && absReal(qureg.stateVec.imag[i] - imag(vec[startInd+i])) < REAL_EPS);
    
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAreEqual = areEqual;
#ifdef MPI_VERSION
    if (qureg.numChunks > 1)
        MPI_Allreduce(&areEqual, &allAreEqual, 1, MPI_INTEGER, MPI_LAND, MPI_COMM_WORLD);
#endif
        
    return allAreEqual;
}

/** hardware-agnostic comparison of the given matrix and mixed qureg, to within 
 * the QuEST_PREC-specific REAL_EPS (defined in QuEST_precision) precision.
 * In GPU mode, this involves a GPU to CPU memory copy overhead.
 * In distributed mode, this involves a all-to-all single-int broadcast
 */
bool areEqual(QMatrix matr, Qureg qureg) {
    REQUIRE( qureg.isDensityMatrix );
    REQUIRE( matr.size()*matr.size() == qureg.numAmpsTotal );
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    
    // the starting index in vec of this node's qureg partition.
    long long int startInd = qureg.chunkId * qureg.numAmpsPerChunk;
    long long int globalInd, row, col;
    
    // loop terminates when areEqual = 0
    int areEqual = 1;
    for (long long int i=0; areEqual && i<qureg.numAmpsPerChunk; i++) {
        globalInd = startInd + i;
        row = globalInd % matr.size();
        col = globalInd / matr.size();
        areEqual = (
               absReal(qureg.stateVec.real[i] - real(matr[row][col])) < REAL_EPS
            && absReal(qureg.stateVec.imag[i] - imag(matr[row][col])) < REAL_EPS);
    }
    
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAreEqual = areEqual;
#ifdef MPI_VERSION
    if (qureg.numChunks > 1)
        MPI_Allreduce(&areEqual, &allAreEqual, 1, MPI_INTEGER, MPI_LAND, MPI_COMM_WORLD);
#endif
        
    return allAreEqual;
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
    REQUIRE( src.real != NULL );
    REQUIRE( src.imag != NULL );
    QMatrix dest = getZeroMatrix(1 << src.numQubits);
    macro_copyComplexMatrix(dest, src);
    return dest;
}

/** Encodes a Complex pair into the matrix (a=alpha, b=beta)
 * {{a, -conj(b)}},
 * {{b,  conj(a)}}
 */
QMatrix toQMatrix(Complex alpha, Complex beta) {
    qcomp a = qcomp(alpha.real, alpha.imag);
    qcomp b = qcomp(beta.real, beta.imag);
    QMatrix matr{
        {a, -conj(b)},
        {b,  conj(a)}};
    return matr;
}

/** Encodes a density-matrix qureg into a QMatrix.
 * In GPU mode, this involves a GPU to CPU memory copy overhead.
 * In distributed mode, this involves a all-to-all full-statevector broadcast
 */
QMatrix toQMatrix(Qureg qureg) {
    REQUIRE( qureg.isDensityMatrix );
    REQUIRE( qureg.numAmpsTotal < MPI_MAX_AMPS_IN_MSG );
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    
    qreal* fullRe;
    qreal* fullIm;
    
    // in distributed mode, give every node the full state vector
#ifdef MPI_VERSION
    fullRe = malloc(qureg.numAmpsTotal * sizeof *fullRe);
    fullIm = malloc(qureg.numAmpsTotal * sizeof *fullIm);
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
    
    // clean up
#ifdef MPI_VERSION
    free(fullRe);
    free(fullIm);
#endif
    return matr;
}

/** Encodes a state-vector qureg into a QVector.
 * In GPU mode, this involves a GPU to CPU memory copy overhead.
 * In distributed mode, this involves a all-to-all full-statevector broadcast
 */
QVector toQVector(Qureg qureg) {
    REQUIRE( !qureg.isDensityMatrix );
    REQUIRE( qureg.numAmpsTotal < MPI_MAX_AMPS_IN_MSG );
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    
    qreal* fullRe;
    qreal* fullIm;
    
    // in distributed mode, give every node the full state vector
#ifdef MPI_VERSION
    fullRe = malloc(qureg.numAmpsTotal * sizeof *fullRe);
    fullIm = malloc(qureg.numAmpsTotal * sizeof *fullIm);
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
    
    // clean up
#ifdef MPI_VERSION
    free(fullRe);
    free(fullIm);
#endif
    return vec;
}

