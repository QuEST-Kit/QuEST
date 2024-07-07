/** @file
 * Miscellaneous utility functions needed internally.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/structures.h"

#include "quest/src/core/errors.hpp"



/*
 * MATRIX CONJUGATION
 */

template <typename T>
void setElemsConj(T& matrix) {
    for (qindex i=0; i<matrix.numRows; i++)
        for (qindex j=0; j<matrix.numRows; j++)
            matrix.elems[i][j] = conj(matrix.elems[i][j]);
}

void util_setConj(CompMatr matrix) {
    setElemsConj(matrix);
}

CompMatr1 util_getConj(CompMatr1 matrix) {
    CompMatr1 conj = matrix;
    setElemsConj(conj);
    return conj;
}

CompMatr2 util_getConj(CompMatr2 matrix) {
    CompMatr2 conj = matrix;
    setElemsConj(conj);
    return conj;
}



/*
 * MATRIX UNITARITY
 */

bool isUnitary(qcomp** matrix, qindex dim) {

    // matrix is simply qcomp** but inherits the consts
    // from CompMatr.elems, which protects changing
    // the memory addresses but permits element change

    qreal epsSq = VALIDATION_EPSILON * VALIDATION_EPSILON;

    // check m * dagger(m) == identity
    for (qindex r=0; r<dim; r++) {
        for (qindex c=0; c<dim; c++) {

            // compute m[r,...] * dagger(m)[...,c]
            qcomp elem = 0;
            for (qindex i=0; i<dim; i++)
                elem += matrix[r][i] * conj(matrix[c][i]);

            // check if further than epsilon from identity[r,c]
            qcomp dif = elem - qcomp(r == c, 0);
            qreal dist = real(dif)*real(dif) + imag(dif)*imag(dif);
            if (dist > epsSq)
                return false;
        }
    }

    return true;
}

bool util_isUnitary(CompMatr1 matrix) {

    qcomp* ptrs[] = {
        matrix.elems[0], 
        matrix.elems[1]
    };

    return isUnitary(ptrs, matrix.numRows);
}

bool util_isUnitary(CompMatr2 matrix) {

    qcomp* ptrs[] = {
        matrix.elems[0], 
        matrix.elems[1], 
        matrix.elems[2], 
        matrix.elems[3]
    };

    return isUnitary(ptrs, matrix.numRows);
}

bool util_isUnitary(CompMatr matrix) {

    return isUnitary(matrix.elems, matrix.numRows);
}



/*
 * QUBIT SHIFTING
 */

int util_getShifted(int qubit, Qureg qureg) {
    assert_shiftedQuregIsDensMatr(qureg);
    
    return qubit + qureg.numQubits;
}