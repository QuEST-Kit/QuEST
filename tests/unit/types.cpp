/** @file
 * Unit tests of the types module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unittypes Types unit tests
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

#define TEST_CATEGORY "[unit][types]"



/** 
 * TESTS
 * 
 * @ingroup unittypes
 * @{
 */

TEST_CASE( "placeholder", TEST_CATEGORY) {

}

/** @} (end defgroup) */



/*
 * @todo
 * UNTESTED FUNCTIONS BELOW
 */


static inline qcomp getQcomp(qreal re, qreal im);

void reportQcomp(qcomp num);

// qcomp literals

// qcomp arithmatic overloads
