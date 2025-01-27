#ifndef MACROS_HPP
#define MACROS_HPP

#include "catch.hpp"

// TODO: 
// use CMake to obtain catch2 and replace #include "catch.hpp"
// with "#include <catch2/catch_test_macros.hpp>""



// TODO:
// adjust TEST_EPSILON according to QuEST's FLOAT_PRECISION

#define TEST_EPSILON 1E-10


/*
 * preconditions to the internal unit testing functions are checked using 
 * DEMAND rather than Catch2's REQUIRE, so that they are not counted in the 
 * total unit testing statistics (e.g. number of checks passed).
 */
#define DEMAND( cond ) if (!(cond)) FAIL( );


#endif // MACROS_HPP