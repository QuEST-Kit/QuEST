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
