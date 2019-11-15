/** @file
 * Unoptimised, analytic implementations of matrix operations used by QuEST_unit_tests
 *
 * @author Tyson Jones
 */

#include "QuEST.h"
#include "QuEST_test_utils.hpp"
#include "catch.hpp"
#include <random>
#include <algorithm>



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

/** returns log2 of numbers which must be gauranteed to be 2^n */
unsigned int calcLog2(unsigned int res) {
    unsigned int n = 0;
    while (res >>= 1)
        n++;
    return n;
}

/** overwrites state to be the result of applying the unitary matrix op (with the
 * specified control and target qubits) to statevector state, i.e. fullOp(op) * |state>
 */
void applyUnitaryOp(
    QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    int numQubits = calcLog2(state.size());
    QMatrix fullOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    state = getMatrixVectorProduct(fullOp, state);
}
void applyUnitaryOp(
    QVector &state, int* ctrls, int numCtrls, int target, QMatrix op
) {
    int targs[1] = {target};
    applyUnitaryOp(state, ctrls, numCtrls, targs, 1, op);
}
void applyUnitaryOp(
    QVector &state, int *targs, int numTargs, QMatrix op
) {
    applyUnitaryOp(state, NULL, 0, targs, numTargs, op);
}
void applyUnitaryOp(
    QVector &state, int ctrl, int targ, QMatrix op
) {
    int ctrls[1] = {ctrl};
    int targs[1] = {targ};
    applyUnitaryOp(state, ctrls, 1, targs, 1, op);
}
void applyUnitaryOp(
    QVector &state, int targ, QMatrix op
) {
    int targs[1] = {targ};
    applyUnitaryOp(state, NULL, 0, targs, 1, op);
}

/** overwrites state to be the result of applying the unitary matrix op (with the
 * specified control and target qubits) to density matrix state, i.e. 
 *  fullOp(op) * state * fullOp(op)^dagger
 */
void applyUnitaryOp(
    QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op
) {
    int numQubits = calcLog2(state.size());
    QMatrix leftOp = getFullOperatorMatrix(ctrls, numCtrls, targs, numTargs, op, numQubits);
    QMatrix rightOp = getConjugateTranspose(leftOp);
    state = getMatrixProduct(getMatrixProduct(leftOp, state), rightOp);
}
void applyUnitaryOp(
    QMatrix &state, int* ctrls, int numCtrls, int target, QMatrix op
) {
    int targs[1] = {target};
    applyUnitaryOp(state, ctrls, numCtrls, targs, 1, op);
}
void applyUnitaryOp(
    QMatrix &state, int *targs, int numTargs, QMatrix op
) {
    applyUnitaryOp(state, NULL, 0, targs, numTargs, op);
}
void applyUnitaryOp(
    QMatrix &state, int ctrl, int targ, QMatrix op
) {
    int ctrls[1] = {ctrl};
    int targs[1] = {targ};
    applyUnitaryOp(state, ctrls, 1, targs, 1, op);
}
void applyUnitaryOp(
    QMatrix &state, int targ, QMatrix op
) {
    int targs[1] = {targ};
    applyUnitaryOp(state, NULL, 0, targs, 1, op);
}

/** hardware-agnostic comparison of the given vector and pure qureg, to within 
 * the QuEST_PREC-specific REAL_EPS (defined in QuEST_precision) precision.
 * In GPU mode, this involves a GPU to CPU memory copy overhead.
 * In distributed mode, this involves a all-to-all single-int broadcast
 */
bool areEqual(Qureg qureg, QVector vec) {
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
 * the QuEST_PREC-specific 1E4*REAL_EPS (defined in QuEST_precision) precision.
 * We check against 1E4*REAL_EPS, since effecting operations on unitaries involve 
 * twice as many operations are on a statevector, and the state is stored in square
 * as many elements.
 * In GPU mode, this involves a GPU to CPU memory copy overhead.
 * In distributed mode, this involves a all-to-all single-int broadcast
 */
bool areEqual(Qureg qureg, QMatrix matr) {
    REQUIRE( qureg.isDensityMatrix );
    REQUIRE( matr.size()*matr.size() == qureg.numAmpsTotal );
    
    // ensure local qureg.stateVec is up to date
    copyStateFromGPU(qureg);
    
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
        ampsAgree = (realDif < 1E4*REAL_EPS && imagDif < 1E4*REAL_EPS);

        // break loop as soon as amplitudes disagree
        if (!ampsAgree)
            break;
            
        /* TODO:
         * of the nodes which disagree, the lowest-rank should send its 
         * disagreeing (i, row, col, stateVec[i]) to rank 0 which should 
         * report it immediately (before the impending REQUIRE failure)
         * using FAIL_CHECK, so users can determine nature of disagreement 
         * (e.g. numerical precision).
         * Note FAIL_CHECK accepts << like cout, e.g.
         * FAIL_CHECK( "Amp at (" << row << ", " << col ") disagreed" );
         */
    }
    
    // if one node's partition wasn't equal, all-nodes must report not-equal
    int allAmpsAgree = ampsAgree;
#ifdef MPI_VERSION
    if (qureg.numChunks > 1)
        MPI_Allreduce(&ampsAgree, &allAmpsAgree, 1, MPI_INTEGER, MPI_LAND, MPI_COMM_WORLD);
#endif
        
    return allAmpsAgree;
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
    REQUIRE( qm.size() == 2 );
    ComplexMatrix2 cm;
    macro_copyQMatrix(cm, qm);
    return cm;
}
ComplexMatrix4 toComplexMatrix4(QMatrix qm) {
    REQUIRE( qm.size() == 4 );
    ComplexMatrix4 cm;
    macro_copyQMatrix(cm, qm);
    return cm;
}
void toComplexMatrixN(QMatrix qm, ComplexMatrixN cm) {
    REQUIRE( qm.size() == (1<<cm.numQubits) );
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

/** Generates every fixed-length sublist of the constructor-given list, in increasing 
 * lexographic order. That is, generates every combination of the given list 
 * and every permutation of each. If the sublist length is the full list length, 
 * this generator produces every permutation correctly. Note that the (same) pointer 
 * returned by get()must not be modified between invocations of next(). QuEST's 
 * internal functions will indeed modify but restore the qubit index lists given 
 * to them, which is ok. Assumes the constructor-given list contains no duplicate, 
 * otherwise the generated sublists may be duplicated.
 */
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
        
        REQUIRE( numSamps <= numElems );
                
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
        int numSamps, int* exclude, int numExclude
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
        
        REQUIRE( numSamps <= len );

        // prepare sublist
        sublen = numSamps;
        createSublist();
        prepareSublist();
    }
    
    int* const& get() const {
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
    Catch::Generators::GeneratorWrapper<int>&& gen, int numSamps, int* exclude, int numExclude
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