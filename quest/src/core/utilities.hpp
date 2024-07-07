/** @file
 * Miscellaneous utility functions needed internally.
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/structures.h"

#include "quest/src/core/errors.hpp"



/*
 * MATRIX CONJUGATION
 */

void util_setConj(CompMatr matrix);

CompMatr1 util_getConj(CompMatr1 matrix);

CompMatr2 util_getConj(CompMatr2 matrix);



/*
 * MATRIX UNITARITY
 */

bool util_isUnitary(CompMatr1 matrix);

bool util_isUnitary(CompMatr2 matrix);

bool util_isUnitary(CompMatr matrix);



/*
 * QUBIT SHIFTING
 */

int util_getShifted(int qubit, Qureg qureg);



#endif // UTILITIES_HPP