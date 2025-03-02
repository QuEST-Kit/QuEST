#ifndef MACROS_HPP
#define MACROS_HPP

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>


// TODO:
// should this even be a macro? eh
#define NUM_UNIT_QUREG_QUBITS 6



// TODO:
// move this (and above # quregs) out of macros.hpp and into 
// a new e.g. 'config.cpp' so we don't have to recompile
// all tests when changing precision and qureg sizes, DAYUM

#if FLOAT_PRECISION == 1
    #define TEST_EPSILON 1E-4
#elif FLOAT_PRECISION == 2
    #define TEST_EPSILON 1E-10
#elif FLOAT_PRECISION == 4
    #define TEST_EPSILON 1E-12
#endif



/*
 * preconditions to the internal unit testing functions are checked using 
 * DEMAND rather than Catch2's REQUIRE, so that they are not counted in the 
 * total unit testing statistics (e.g. number of checks passed). 
 */
#define DEMAND( cond ) do { if (!(cond)) { FAIL( ); } } while (0)


// section labels

#define LABEL_CORRECTNESS "correctness"
#define LABEL_VALIDATION "validation"
#define LABEL_STATEVEC "statevector"
#define LABEL_DENSMATR "densitymatrix"


#endif // MACROS_HPP
