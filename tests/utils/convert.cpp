/** @file
 * Testing utilities for converting QuEST API structures
 * (like Qureg, CompMatr, PauliStr) to/from testing types 
 * (like qvector and qmatrix).
 *
 * @author Tyson Jones
 */

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "linalg.hpp"
#include "macros.hpp"
#include "lists.hpp"

#include "quest/include/quest.h"

#include <type_traits>
using std::is_same_v;

#include <vector>
using std::vector;



/*
 * TO QUREG
 */


void setQuregToReference(Qureg qureg, qvector vector) {
    DEMAND( !qureg.isDensityMatrix );
    DEMAND( qureg.numAmps == (qindex) vector.size() );

    setQuregAmps(qureg, 0, vector.data(), vector.size());
}


void setQuregToReference(Qureg qureg, qmatrix matrix) {
    DEMAND( qureg.isDensityMatrix );
    DEMAND( getPow2(qureg.numQubits) == (qindex) matrix.size() );

    qindex numRows = 1;
    qindex numCols = getPow2(qureg.numQubits);

    for (size_t r=0; r<matrix.size(); r++) {
        qcomp* arr[] = {matrix[r].data()};
        setDensityQuregAmps(qureg, r, 0, arr, numRows, numCols);
    }
}



/*
 * FROM QUREG
 */


qvector getVector(Qureg qureg) {
    DEMAND( !qureg.isDensityMatrix );

    qvector out = getZeroVector(qureg.numAmps);
    getQuregAmps(out.data(), qureg, 0, qureg.numAmps);
    return out;
}


qmatrix getMatrix(Qureg qureg) {
    DEMAND( qureg.isDensityMatrix );

    qindex numRows = 1;
    qindex numCols = getPow2(qureg.numQubits);
    qmatrix out = getZeroMatrix(numCols);

    for (size_t r=0; r<out.size(); r++) {
        qcomp* arr[] = {out[r].data()};
        getDensityQuregAmps(arr, qureg, r, 0, numRows, numCols);

    }

    return out;
}



/*
 * FROM API MATRIX
 */


template <typename T> 
qcomp getElem(T m, size_t r, size_t c) {

    if constexpr (is_same_v<T, CompMatr> || is_same_v<T, SuperOp>)
        return m.cpuElems[r][c];
    
    if constexpr (is_same_v<T, CompMatr1> || is_same_v<T, CompMatr2>)
        return m.elems[r][c];

    if constexpr (is_same_v<T, DiagMatr>)
        return (r==c)? m.cpuElems[r] : 0;

    if constexpr (is_same_v<T, DiagMatr1> || is_same_v<T, DiagMatr2>)
        return (r==c)? m.elems[r] : 0;
}


template <typename T> 
qmatrix getMatrixInner(T m) {

    qindex dim = (is_same_v<T, SuperOp>)?
        getPow2(2*m.numQubits):
        getPow2(  m.numQubits);

    qmatrix out = getZeroMatrix(dim);

    for (size_t r=0; r<out.size(); r++)
        for (size_t c=0; c<out.size(); c++)
            out[r][c] = getElem<T>(m, r, c);

    return out;
}


qmatrix getMatrix(CompMatr1 m) { return getMatrixInner(m); }
qmatrix getMatrix(CompMatr2 m) { return getMatrixInner(m); }
qmatrix getMatrix(CompMatr  m) { return getMatrixInner(m); }
qmatrix getMatrix(DiagMatr1 m) { return getMatrixInner(m); }
qmatrix getMatrix(DiagMatr2 m) { return getMatrixInner(m); }
qmatrix getMatrix(DiagMatr  m) { return getMatrixInner(m); }
qmatrix getMatrix(SuperOp   m) { return getMatrixInner(m); }



/*
 * FROM PAULI STRING
 */


extern int paulis_getPauliAt(PauliStr str, int ind);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStrSum sum);


qmatrix getMatrix(PauliStr str, vector<int> targs) {
    DEMAND( targs.size() >= 1 );

    qmatrix out = getIdentityMatrix(1);

    for (auto t : targs) {
        int ind = paulis_getPauliAt(str, t);
        qmatrix matr = getPauliMatrix(ind);
        out = getKroneckerProduct(matr, out);
    }

    return out;
}


qmatrix getMatrix(PauliStr str, int numQubits) {
    DEMAND( numQubits >= paulis_getIndOfLefmostNonIdentityPauli(str) );

    return getMatrix(str, getRange(numQubits));
}


qmatrix getMatrix(PauliStrSum sum, int numQubits) {
    DEMAND( sum.numTerms > 0 );
    DEMAND( numQubits >= paulis_getIndOfLefmostNonIdentityPauli(sum) );

    qmatrix out = getZeroMatrix(getPow2(numQubits));

    for (qindex i=0; i<sum.numTerms; i++)
        out += sum.coeffs[i] * getMatrix(sum.strings[i], numQubits);

    return out;
}
