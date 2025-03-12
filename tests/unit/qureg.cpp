/** @file
 * Unit tests of the qureg module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitqureg Qureg
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

#define TEST_CATEGORY "[unit][qureg]"



/** 
 * TESTS
 * 
 * @ingroup unitqureg
 * @{
 */

TEST_CASE( "placeholder8", TEST_CATEGORY) {

}

/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */


Qureg createQureg(int numQubits);
Qureg createDensityQureg(int numQubits);

Qureg createForcedQureg(int numQubits);
Qureg createForcedDensityQureg(int numQubits);

Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread);

Qureg createCloneQureg(Qureg qureg);

void destroyQureg(Qureg qureg);

void reportQuregParams(Qureg qureg);
void reportQureg(Qureg qureg);

void syncQuregToGpu  (Qureg qureg);
void syncQuregFromGpu(Qureg qureg);

void syncSubQuregToGpu  (Qureg qureg, qindex localStartInd, qindex numLocalAmps);
void syncSubQuregFromGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps);

void getQuregAmps(qcomp* outAmps, Qureg qureg, qindex startInd, qindex numAmps);
void getDensityQuregAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


qcomp getQuregAmp(Qureg qureg, qindex index);

qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column);
