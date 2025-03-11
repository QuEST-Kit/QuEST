/** @file
 * Unit tests of the initialisations module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitinit Initialisation unit tests
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

#define TEST_CATEGORY "[unit][initialisations]"



/** 
 * TESTS
 * 
 * @ingroup unitinit
 * @{
 */

TEST_CASE( "placeholder", TEST_CATEGORY) {

}

/** @} (end defgroup) */



/*
 * TODO:
 * UNTESTED FUNCTIONS BELOW
 */

void initBlankState(Qureg qureg);

void initZeroState(Qureg qureg);

void initPlusState(Qureg qureg);

void initPureState(Qureg qureg, Qureg pure);

void initClassicalState(Qureg qureg, qindex stateInd);

void initDebugState(Qureg qureg);

void initArbitraryPureState(Qureg qureg, qcomp* amps);

void initRandomPureState(Qureg qureg);

void initRandomMixedState(Qureg qureg, qindex numPureStates);


void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);

void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);

void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);

void setQuregToClone(Qureg targetQureg, Qureg copyQureg);

void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2);

qreal setQuregToRenormalized(Qureg qureg);

void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum);


void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);

void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits);
