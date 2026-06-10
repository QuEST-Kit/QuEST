/** @file
 * Unit tests of the environment module.
 *
 * @author Oliver Brown
 * @author Tyson Jones
 * 
 * @defgroup unitexperi Experimental
 * @ingroup unittests
 */

#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/config.hpp"

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[experimental]"



/** 
 * TESTS
 * 
 * @ingroup unitexperi
 * @{
 */


TEST_CASE( "setQuESTNumGpuThreadsPerBlock", TEST_CATEGORY ) {

    // remember the default number for later restoration (hence static)
    static int initNumTPB = getQuESTNumGpuThreadsPerBlock();

    SECTION( LABEL_CORRECTNESS ) {

        // begin at 64 (AMD min, larger than NVIDIA min of 32),
        // stop at 1024 (should be less than dev-specific max)
        int inNumTPB = GENERATE( 64, 128, 256, 512, 1024 ); 
        setQuESTNumGpuThreadsPerBlock(inNumTPB);

        int outNumTPB = getQuESTNumGpuThreadsPerBlock();
        REQUIRE( inNumTPB == outNumTPB );
        
        // BEWARE that we do not here test whether all QuEST
        // operators succeed with the various numTBP; that must
        // be ad hoc asssesed via updating the numTBP env-var
        // before launching the entirety of the tests
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "Negative" ) {

            int badNumTPB = GENERATE( 0, -1, -9999 );
            REQUIRE_THROWS_WITH( setQuESTNumGpuThreadsPerBlock(badNumTPB), ContainsSubstring( "must be positive" ) );
        }

        SECTION( "Indivisible by warp size" ) {

            // If HIP status was attached to QuESTEnv, we could do:
            //     QuESTEnv env = getQuESTEnv();
            //     int warpSize = (env.isGpuAccelerated && env.isHipCompiled)? 64 : 32;
            // Since this currently isn't the case, we assume a warp size of 32,
            // which will mean when this test is run on AMD GPUs, the below tested
            // badNumTBP won't be as interestingly/rigorously spread
            int warpSize = 32;

            int badNumTPB = GENERATE_COPY( warpSize - 1, warpSize + 1, warpSize + warpSize/2, 3*warpSize + warpSize/2 );

            REQUIRE_THROWS_WITH( setQuESTNumGpuThreadsPerBlock(badNumTPB), ContainsSubstring( "does not divide evenly into the warp size" ) );
        }

        SECTION( "Exceeds device maximum" ) {

            int badNumTPB = 999999; // exceeds expected 1024 max

            // Cannot be tested (since validation not imposed) when GPU is not actively used
            if (getQuESTEnv().isGpuAccelerated)
                REQUIRE_THROWS_WITH( setQuESTNumGpuThreadsPerBlock(badNumTPB), ContainsSubstring( "Exceeds the hardware-imposed maximum" ) );

            SUCCEED( );
        }
    }

    // restore numTBP, so as not to interfere with other tests
    setQuESTNumGpuThreadsPerBlock(initNumTPB);
}


TEST_CASE( "getQuESTNumGpuThreadsPerBlock", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // check initial value matches either the env-var (if set),
        // or the fixed default in the codebase (hardcoded in test utils)
        int defaultNum = getDefaultNumGpuThreadsPerBlock(); // test util via env-var
        int reportedNum = getQuESTNumGpuThreadsPerBlock();  // QuEST API

        REQUIRE( defaultNum == reportedNum );

        // further testing of this function appears in setQuESTNumGpuThreadsPerBlock()
    }

    SECTION( LABEL_VALIDATION ) {

        // there is none (except untestable env is init!)
        SUCCEED( );
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

// nothing! :^)
