/** @file
 * Macros used by the tests and testing utilities.
 *
 * @author Tyson Jones
 */

#ifndef MACROS_HPP
#define MACROS_HPP

#include <catch2/catch_test_macros.hpp>


// TODO:
// this shouldn't be a macro; Catch2 can handle fine when
// we make this a runtime variable (changing some of our
// GENERATE to GENERATE_COPY). As a macro, changing the
// qureg size requires recompiling the entire tests, which
// is certainly annoying for development!

#define NUM_UNIT_QUREG_QUBITS 6



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
