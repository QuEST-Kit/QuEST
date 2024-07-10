/** @file
 * Miscellaneous utility functions needed internally.
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/structures.h"

#include "quest/src/core/errors.hpp"

#include <type_traits>



/*
 * MATRIX TYPING
 *
 * defined here in the header since templated, and which use compile-time inspection
 */

// T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr
template<class T>
constexpr bool isDenseMatrixType() {

    // CompMatr are "dense", storing all 2D elements
    if constexpr (
        std::is_same_v<T, CompMatr1> ||
        std::is_same_v<T, CompMatr2> ||
        std::is_same_v<T, CompMatr>
    )
        return true;

    // DiagMatr are "sparse", storing only the diagonals
    if constexpr (
        std::is_same_v<T, DiagMatr1> ||
        std::is_same_v<T, DiagMatr2> ||
        std::is_same_v<T, DiagMatr>
    )
        return false;

    // this line is unreachable but throwing errors in a template expansion is ludicrous;
    // above type checks are explicit in case we add more matrix types later
    return false;
}

// T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr
template<class T>
constexpr bool isFixedSizeMatrixType() {

    return (
        std::is_same_v<T, CompMatr1> ||
        std::is_same_v<T, CompMatr2> ||
        std::is_same_v<T, DiagMatr1> ||
        std::is_same_v<T, DiagMatr2>
    );
}

// T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr
template<class T>
constexpr qindex getMatrixDim(T matr) {
    
    if constexpr (isDenseMatrixType<T>())
        return matr.numRows;
    else
        return matr.numElems;
}



/*
 * MATRIX CONJUGATION
 */

CompMatr1 util_getConj(CompMatr1 matrix);
CompMatr2 util_getConj(CompMatr2 matrix);
DiagMatr1 util_getConj(DiagMatr1 matrix);
DiagMatr2 util_getConj(DiagMatr2 matrix);

void util_setConj(CompMatr matrix);
void util_setConj(DiagMatr matrix);



/*
 * MATRIX UNITARITY
 */

bool util_isUnitary(CompMatr1 matrix);
bool util_isUnitary(CompMatr2 matrix);
bool util_isUnitary(CompMatr matrix);
bool util_isUnitary(DiagMatr1 matrix);
bool util_isUnitary(DiagMatr2 matrix);
bool util_isUnitary(DiagMatr matrix);



/*
 * QUBIT SHIFTING
 */

int util_getShifted(int qubit, Qureg qureg);



#endif // UTILITIES_HPP