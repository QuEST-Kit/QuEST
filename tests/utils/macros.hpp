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


/**
 * macros which affect the speed and rigour of the unit tests, useful
 * for accelerating tests on particular platforms (e.g. paid github runners).
 * The default values are those which perform the most rigorous tests at
 * the slowest speed, so adjusting these macros accelerates tests. 
 *
 * @todo
 * These are clunky preprocessors (invoking full recompilation when changed),
 * rather than runtime arguments, because of the nuisance of passing such
 * args to cmake. It can be done however using environment variables; see
 * https://stackoverflow.com/questions/28812533/
 */

// 0 = perform all, and a sensible value to accelerate tests is 50
#ifndef TEST_MAX_NUM_QUBIT_PERMUTATIONS
#define TEST_MAX_NUM_QUBIT_PERMUTATIONS 0
#endif

// 0 = perform all (very slow), while 4 limits to superops = 8-qubit matrices
#ifndef TEST_MAX_NUM_SUPEROP_TARGETS
#define TEST_MAX_NUM_SUPEROP_TARGETS 4
#endif

// 0 = use all available deployments at once, 1 = try all combinations in-turn
#ifndef TEST_ALL_DEPLOYMENTS
#define TEST_ALL_DEPLOYMENTS 1
#endif

// number of times to repeat each "[mixed]" test (minimum 1)
#ifndef TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS
#define TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS 10
#endif

// spoofing above macros as consts to doc
#if 0

    /// @macrodoc
    const int TEST_MAX_NUM_QUBIT_PERMUTATIONS = 0;

    /// @macrodoc
    const int TEST_MAX_NUM_SUPEROP_TARGETS = 4;

    /// @macrodoc
    const int TEST_ALL_DEPLOYMENTS = 1;

    /// @macrodoc
    const int TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS = 10;

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
#define LABEL_C_INTERFACE "C interface"
#define LABEL_CPP_INTERFACE "C++ interface"

#define LABEL_DELIMITER ", "

#define LABEL_UNIT_TAG "[unit]"
#define LABEL_MIXED_DEPLOY_TAG "[mixed]"
#define LABEL_INTEGRATION_TAG "[integration]"


// detect LLVM address sanitizer (on GCC and Clang only)
#if defined(__SANITIZE_ADDRESS__)
    #define SANITIZER_IS_ACTIVE
#elif defined(__has_feature)
    #if __has_feature(address_sanitizer)
        #define SANITIZER_IS_ACTIVE
    #endif
#endif


#endif // MACROS_HPP

/** @} (end defgroup) */
