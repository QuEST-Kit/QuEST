/** @file
 * Unit tests of the environment module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitenv Environment
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[environment]"



/** 
 * TESTS
 * 
 * @ingroup unitenv
 * @{
 */


TEST_CASE( "initQuESTEnv", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // cannot be meaningfully tested since env already active
        SUCCEED( );
    }

    SECTION( LABEL_VALIDATION ) {

        REQUIRE_THROWS_WITH( initQuESTEnv(), ContainsSubstring( "already been initialised") );
    }
}


TEST_CASE( "initCustomQuESTEnv", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // cannot be meaningfully tested since env already active
        SUCCEED( );
    }

    SECTION( LABEL_VALIDATION ) {

        REQUIRE_THROWS_WITH( initCustomQuESTEnv(0,0,0), ContainsSubstring( "already been initialised") );

        // cannot check arguments since env-already-initialised
        // validation is performed first
    }
}


TEST_CASE( "finalizeQuESTEnv", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // cannot be meaningfully tested since calling
        // mid-tests will break subsequent testing
        SUCCEED( );
    }

    SECTION( LABEL_VALIDATION ) {

        // cannot be validated since calling it once
        // is always valid but would break subsequent tests
        SUCCEED( );
    }
}


TEST_CASE( "syncQuESTEnv", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // always legal to call
        REQUIRE_NOTHROW( syncQuESTEnv() );
    }

    SECTION( LABEL_VALIDATION ) {

        // cannot test validation (that env is already
        // created) since we cannot call before initQuESTEnv
        SUCCEED( );
    }
}


TEST_CASE( "isQuESTEnvInit", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // cannot test for other outcome
        REQUIRE( isQuESTEnvInit() == 1 );
    }

    SECTION( LABEL_VALIDATION ) {

        // performs no validation
        SUCCEED( );
    }
}


TEST_CASE( "getQuESTEnv", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        QuESTEnv env = getQuESTEnv();

        REQUIRE( (env.isMultithreaded  == 0 || env.isMultithreaded  == 1) );
        REQUIRE( (env.isGpuAccelerated == 0 || env.isGpuAccelerated == 1) );
        REQUIRE( (env.isDistributed    == 0 || env.isDistributed    == 1) );
        
        REQUIRE( env.rank     >= 0 );
        REQUIRE( env.numNodes >= 0 );
        
        if (!env.isDistributed) {
            REQUIRE( env.rank     == 0 );
            REQUIRE( env.numNodes == 1 );
        }

        bool isNumNodesPow2 = ((env.numNodes & (env.numNodes - 1)) == 0);
        REQUIRE( isNumNodesPow2 );
    }

    SECTION( LABEL_VALIDATION ) {

        // performs no validation
        SUCCEED( );
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void reportQuESTEnv();
