/** @file
 * Unit tests of the environment module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitenv Environment unit tests
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

#define TEST_CATEGORY "[unit][environment]"



/** 
 * TESTS
 * 
 * @ingroup unitenv
 * @{
 */

TEST_CASE( "placeholder", TEST_CATEGORY) {

}

/** @} (end defgroup) */



/*
 * @todo
 * UNTESTED FUNCTIONS BELOW
 */

void initQuESTEnv();

void initCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread);

void finalizeQuESTEnv();

void syncQuESTEnv();

void reportQuESTEnv();

int isQuESTEnvInit();

QuESTEnv getQuESTEnv();