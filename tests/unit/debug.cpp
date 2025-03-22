/** @file
 * Unit tests of the debug module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitdebug Debug
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/macros.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[debug]"



/** 
 * TESTS
 * 
 * @ingroup unitdebug
 * @{
 */


TEST_CASE( "setMaxNumReportedSigFigs", TEST_CATEGORY) {

    SECTION( LABEL_CORRECTNESS ) {

        qcomp scalar = getQcomp(0.12345, 0.12345);

        std::vector<std::string> refs = {
            "0.1+0.1i\n",
            "0.12+0.12i\n",
            "0.123+0.123i\n",
            "0.1235+0.1235i\n", // rounded
            "0.12345+0.12345i\n"
        };

        for (size_t numSigFigs=1; numSigFigs<=refs.size(); numSigFigs++) {

            setMaxNumReportedSigFigs(numSigFigs);

            // redirect stdout to buffer
            std::stringstream buffer;
            std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());

            reportQcomp(scalar);

            // restore stdout
            std::cout.rdbuf(old);
            std::string ref = refs[numSigFigs-1];

            CAPTURE( numSigFigs, ref );
            REQUIRE( buffer.str() == ref );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "number" ) {

            int num = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( setMaxNumReportedSigFigs(num), ContainsSubstring("Cannot be less than one") );
        }
    }   
}
 

/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void setSeeds(unsigned* seeds, int numSeeds);
void setSeedsToDefault();

void getSeeds(unsigned* seeds);
int getNumSeeds();

// void invalidQuESTInputError(const char* msg, const char* func);

void setValidationOn();
void setValidationOff();

void setValidationEpsilonToDefault();
void setValidationEpsilon(qreal eps);
qreal getValidationEpsilon();

void setMaxNumReportedItems(qindex numRows, qindex numCols);

qindex getGpuCacheSize();
void clearGpuCache();

void getEnvironmentString(char str[200]);
