/** @file
 * Unit tests of the channels module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitchannels Channels
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>



/*
 * UTILITIES
 */

#define TEST_CATEGORY "[unit][channels]"



/** 
 * TESTS
 * 
 * @ingroup unitchannels
 * @{
 */
 
TEST_CASE( "placeholder", TEST_CATEGORY) {
 
}
 
/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */


KrausMap createKrausMap(int numQubits, int numOperators);

void syncKrausMap(KrausMap map);

void destroyKrausMap(KrausMap map);

void reportKrausMap(KrausMap map);

void setKrausMap(KrausMap map, qcomp*** matrices);

void setKrausMap(KrausMap map, std::vector<std::vector<std::vector<qcomp>>> matrices);

void setInlineKrausMap(KrausMap map, int numQb, int numOps, std::vector<std::vector<std::vector<qcomp>>> matrices);

KrausMap createInlineKrausMap(int numQubits, int numOperators, std::vector<std::vector<std::vector<qcomp>>> matrices);


SuperOp createSuperOp(int numQubits);

void syncSuperOp(SuperOp op);

void destroySuperOp(SuperOp op);

void reportSuperOp(SuperOp op);

void setSuperOp(SuperOp op, qcomp** matrix);

void setSuperOp(SuperOp op, std::vector<std::vector<qcomp>> matrix);

void setInlineSuperOp(SuperOp op, int numQb, std::vector<std::vector<qcomp>> matrix);

SuperOp createInlineSuperOp(int numQubits, std::vector<std::vector<qcomp>> matrix);
