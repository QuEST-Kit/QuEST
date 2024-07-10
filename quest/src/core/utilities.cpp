/** @file
 * Miscellaneous utility functions needed internally.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/utilities.hpp"

#include <complex>



/*
 * MATRIX CONJUGATION
 */

// type T can be qcomp** or qcomp*[]
template <typename T>
void setDenseElemsConj(T elems, qindex dim) {
    for (qindex i=0; i<dim; i++)
        for (qindex j=0; j<dim; j++)
           elems[i][j] = conj(elems[i][j]);
}

// diagonals don't need templating because arrays decay to pointers, yay!
void setDiagElemsConj(qcomp* elems, qindex dim) {
    for (qindex i=0; i<dim; i++)
        elems[i] = conj(elems[i]);
}

CompMatr1 util_getConj(CompMatr1 matrix) {
    CompMatr1 conj = matrix;
    setDenseElemsConj(conj.elems, matrix.numRows);
    return conj;
}
CompMatr2 util_getConj(CompMatr2 matrix) {
    CompMatr2 conj = matrix;
    setDenseElemsConj(conj.elems, matrix.numRows);
    return conj;
}

DiagMatr1 util_getConj(DiagMatr1 matrix) {
    DiagMatr1 conj = matrix;
    setDiagElemsConj(conj.elems, matrix.numElems);
    return conj;
}
DiagMatr2 util_getConj(DiagMatr2 matrix) {
    DiagMatr2 conj = matrix;
    setDiagElemsConj(conj.elems, matrix.numElems);
    return conj;
}

void util_setConj(CompMatr matrix) {
    setDenseElemsConj(matrix.cpuElems, matrix.numRows);
}
void util_setConj(DiagMatr matrix) {
    setDiagElemsConj(matrix.cpuElems, matrix.numElems);
}



/*
 * MATRIX UNITARITY
 */

// type T can be qcomp** or qcomp*[]
template <typename T>
bool isUnitary(T elems, qindex dim) {

    qreal epsSq = VALIDATION_EPSILON * VALIDATION_EPSILON;

    // check m * dagger(m) == identity
    for (qindex r=0; r<dim; r++) {
        for (qindex c=0; c<dim; c++) {

            // compute m[r,...] * dagger(m)[...,c]
            qcomp elem = 0;
            for (qindex i=0; i<dim; i++)
                elem += elems[r][i] * conj(elems[c][i]);

            // check if further than epsilon from identity[r,c]
            qcomp dif = elem - qcomp(r == c, 0);
            qreal dist = real(dif)*real(dif) + imag(dif)*imag(dif);
            if (dist > epsSq)
                return false;
        }
    }

    return true;
}

// diagonal version doesn't need templating because array decays to pointer, yay!
bool isUnitary(qcomp* diags, qindex dim) {

    // check every element has unit magnitude
    for (qindex i=0; i<dim; i++) {
        qreal mag = std::abs(diags[i]);
        qreal dif = std::abs(1 - mag);

        if (dif > VALIDATION_EPSILON)
            return false;
    }

    return true;
}

bool util_isUnitary(CompMatr1 matrix) {
    return isUnitary(matrix.elems, matrix.numRows);
}
bool util_isUnitary(CompMatr2 matrix) {
    return isUnitary(matrix.elems, matrix.numRows);
}
bool util_isUnitary(CompMatr matrix) {
    return isUnitary(matrix.cpuElems, matrix.numRows);
}

bool util_isUnitary(DiagMatr1 matrix) {
    return isUnitary(matrix.elems, matrix.numElems);
}
bool util_isUnitary(DiagMatr2 matrix) {
    return isUnitary(matrix.elems, matrix.numElems);
}
bool util_isUnitary(DiagMatr matrix) {
    return isUnitary(matrix.cpuElems, matrix.numElems);
}



/*
 * QUBIT SHIFTING
 */

int util_getShifted(int qubit, Qureg qureg) {
    assert_shiftedQuregIsDensMatr(qureg);
    
    return qubit + qureg.numQubits;
}