/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsmacros Macros
 * @ingroup testutils
 * @brief
 * Macros used by the tests and testing utilities.
 * @{
 */

#ifndef MACROS_HPP
#define MACROS_HPP

#include <catch2/catch_test_macros.hpp>


/*
 * macros which affect the speed and rigour of the unit tests, useful
 * for accelerating tests on particular platforms (e.g. paid github runners).
 * The default values are those which perform the most rigorous tests at
 * the slowest speed, so adjusting these macros accelerates tests.
 */

// 0 = perform all, and a sensible value to accelerate tests is 50
#ifndef TEST_MAX_NUM_QUBIT_PERMUTATIONS
#define TEST_MAX_NUM_QUBIT_PERMUTATIONS 0
#endif

// 0 = use all available deployments at once, 1 = try all combinations in-turn
#ifndef TEST_ALL_DEPLOYMENTS
#define TEST_ALL_DEPLOYMENTS 1
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

/** @} (end defgroup) */
