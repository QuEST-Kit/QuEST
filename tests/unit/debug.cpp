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



/*
 * UTILITIES
 */

#define TEST_CATEGORY "[unit][debug]"



/** 
 * TESTS
 * 
 * @ingroup unitdebug
 * @{
 */
 
TEST_CASE( "placeholder", TEST_CATEGORY) {
 
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
void setMaxNumReportedSigFigs(int numSigFigs);

qindex getGpuCacheSize();
void clearGpuCache();

void getEnvironmentString(char str[200]);
