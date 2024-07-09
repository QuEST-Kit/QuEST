/** @file
 * Miscellaneous utility functions needed internally.
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "types.h"
#include "qureg.h"
#include "structures.h"

#include "../core/errors.hpp"



/*
 * MATRIX CONJUGATION
 */

void util_setConj(CompMatrN matrix);

CompMatr1 util_getConj(CompMatr1 matrix);

CompMatr2 util_getConj(CompMatr2 matrix);



/*
 * MATRIX UNITARITY
 */

bool util_isUnitary(CompMatr1 matrix);

bool util_isUnitary(CompMatr2 matrix);

bool util_isUnitary(CompMatrN matrix);



/*
 * QUBIT SHIFTING
 */

int util_getShifted(int qubit, Qureg qureg);



#endif // UTILITIES_HPP