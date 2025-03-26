/** @file
 * Unit tests of the types module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unittypes Types
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[types]"



/** 
 * TESTS
 * 
 * @ingroup unittypes
 * @{
 */


TEST_CASE( "getQcomp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qreal re = getRandomReal(-10, 10);
        qreal im = getRandomReal(-10, 10);

        qcomp comp = getQcomp(re, im);
        REQUIRE( std::real(comp) == re );
        REQUIRE( std::imag(comp) == im );
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


TEST_CASE( "complex arithmetic", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qcomp x;
        x  =   1 + 2_i;
        x +=   3 - 4_i;
        x -= - 5 + 6_i;
        x *= - 7 - 8_i;
        x /=   9 + 10_i;

        qcomp ref = getQcomp(-1303/181., 1126/181.);
        REQUIRE_AGREE( x, ref );
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void reportStr(const char* label);
void reportStr(std::string str);

void reportScalar(const char* label, qcomp num);
void reportScalar(const char* label, qreal num);
void reportScalar(std::string label, qcomp num);
void reportScalar(std::string label, qreal num);