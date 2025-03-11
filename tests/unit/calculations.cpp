/** @file
 * Unit tests of the calculations module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitcalcs Calculations unit tests
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/measure.hpp"

#include <vector>
#include <type_traits>

using std::vector;
using namespace Catch::Matchers;



/*
 * UTILITIES
 */


#define TEST_CATEGORY "[unit][calculations]"


void TEST_ON_CACHED_QUREGS(quregCache quregs, auto refState, auto apiFunc, auto refFunc) {

    for (auto& [label, qureg]: quregs) {

        DYNAMIC_SECTION( label ) {

            // we initialise to a valid pure or mixed state, rather than the 
            // debug state, since many API functions validate the output domain
            // which is often invalid for unnormalised states
            setToRandomState(refState);
            setQureg(qureg, refState);

            // results can be a scalar, a vector, or even a Qureg & qmatrix/qvector
            auto apiResult = apiFunc(qureg);
            auto refResult = refFunc(refState);

            REQUIRE_AGREE( apiResult, refResult );

            // free API result if it is a newly created Qureg
            if constexpr (std::is_same_v<decltype(apiResult),Qureg>)
                destroyQureg(apiResult);
        }
    }
}

void TEST_ON_CACHED_STATEVECS_AND_DENSMATRS(auto apiFunc, auto refFunc) {
    SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), getRefStatevec(), apiFunc, refFunc); }
    SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), getRefDensmatr(), apiFunc, refFunc); }
}

#define TEST_ALL_QUREGS( quregVar, apiExpr, stateVar, refExpr ) \
    TEST_ON_CACHED_STATEVECS_AND_DENSMATRS( \
        [&](Qureg quregVar) { return apiExpr; }, \
        [&](auto& stateVar) { return refExpr; } ) \



/** 
 * TESTS
 * 
 * @ingroup unitcalcs
 * @{
 */


TEST_CASE( "calcExpecPauliStr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_COPY( sublists(range(0,numQubits), numTargs) );

        GENERATE( range(0,10) );
        PauliStr str = getRandomPauliStr(targets);

        TEST_ALL_QUREGS(
            qureg, calcExpecPauliStr(qureg, str),  
            state, real(getReferenceExpectationValue(state, str))
        );
    }

    // @todo input validation
}



TEST_CASE( "calcExpecPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTerms = GENERATE_COPY( 1, numQubits, getPow2(2*numQubits) );

        GENERATE( range(0,100) );
        PauliStrSum sum = createRandomPauliStrSum(numQubits, numTerms);

        TEST_ALL_QUREGS(
            qureg, calcExpecPauliStrSum(qureg, sum),  
            state, real(getReferenceExpectationValue(state, sum))
        );

        destroyPauliStrSum(sum);
    }

    // @todo input validation
}



TEST_CASE( "calcExpecNonHermitianPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,100) );
        int numQubits = getNumCachedQubits();
        int numTerms = GENERATE_COPY( 1, numQubits, getPow2(2*numQubits) );

        PauliStrSum sum = createRandomNonHermitianPauliStrSum(numQubits, numTerms);

        TEST_ALL_QUREGS(
            qureg, calcExpecNonHermitianPauliStrSum(qureg, sum),  
            state, getReferenceExpectationValue(state, sum)
        );

        destroyPauliStrSum(sum);
    }

    // @todo input validation
}



TEST_CASE( "calcProbOfBasisState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qindex index = GENERATE( range(0, 1 << getNumCachedQubits()) );

        TEST_ALL_QUREGS(
            qureg, calcProbOfBasisState(qureg, index),  
            state, getReferenceProbability(state, index)
        );
    }

    // @todo input validation
}



TEST_CASE( "calcProbOfQubitOutcome", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int target = GENERATE( range(0, getNumCachedQubits()) );
        int outcome = GENERATE( 0, 1 );

        TEST_ALL_QUREGS(
            qureg, calcProbOfQubitOutcome(qureg, target, outcome),  
            state, getReferenceProbability(state, {target}, {outcome})
        );
    }

    // @todo input validation
}



TEST_CASE( "calcProbOfMultiQubitOutcome", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTargs = GENERATE_COPY( range(1, numQubits+1) );
        auto targets = GENERATE_COPY( sublists(range(0,numQubits), numTargs) );
        auto outcomes = getRandomInts(0, 1+1, numTargs);

        TEST_ALL_QUREGS(
            qureg, calcProbOfMultiQubitOutcome(qureg, targets.data(), outcomes.data(), numTargs),  
            state, getReferenceProbability(state, targets, outcomes)
        );
    }

    // @todo input validation
}



TEST_CASE( "calcProbsOfAllMultiQubitOutcomes", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_COPY( sublists(range(0,numQubits), numTargs) );

        auto apiFunc = [&](Qureg qureg) { 
            vector<qreal> out(getPow2(numTargs));
            calcProbsOfAllMultiQubitOutcomes(out.data(), qureg, targets.data(), numTargs);
            return out;
        };
        auto refFunc = [&](auto& state) {
            return getAllReferenceProbabilities(state, targets);
        };

        TEST_ON_CACHED_STATEVECS_AND_DENSMATRS( apiFunc, refFunc );
    }

    // @todo input validation
}



TEST_CASE( "calcTotalProb", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,100) );

        TEST_ALL_QUREGS(
            qureg, calcTotalProb(qureg),  
            state, getReferenceProbability(state)
        );
    }

    // @todo input validation
}



TEST_CASE( "calcPurity", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,100) );

        TEST_ALL_QUREGS(
            qureg, calcPurity(qureg),  
            state, getReferencePurity(state)
        );
    }

    // @todo input validation
}



/// @private
int getMaxNumTracedQubits(int numQubits) {

    // cannot reduce all qubits, nor so many that the final qureg is
    // illegally distributed (which we enforce even upon the testsed
    // non-distributed Quregs for simplicity)

    return std::min(
        numQubits - 1,
        numQubits - getLog2(getQuESTEnv().numNodes));
}


TEST_CASE( "calcPartialTrace", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( LABEL_DENSMATR ) {

            int numQubits = getNumCachedQubits();
            int maxNumTargs = getMaxNumTracedQubits(numQubits);
            int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            auto targets = GENERATE_COPY( sublists(range(0,numQubits), numTargs) );

            auto apiFunc = [&](Qureg qureg) { return calcPartialTrace(qureg, targets.data(), numTargs); };
            auto refFunc = [&](qmatrix ref) { return getPartialTrace(ref, targets); };
            
            TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), getRefDensmatr(), apiFunc, refFunc);
        }
    }

    // @todo input validation
}



TEST_CASE( "calcReducedDensityMatrix", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( LABEL_DENSMATR ) {

            int numQubits = getNumCachedQubits();
            int maxNumTraced = getMaxNumTracedQubits(numQubits);
            int minNumRetained = numQubits - maxNumTraced;
            int maxNumRetained = numQubits - 1;
            int numRetained = GENERATE_COPY( range(minNumRetained, maxNumRetained+1) );
            auto retains = GENERATE_COPY( sublists(range(0,numQubits), numRetained) );
            auto targets = getComplement(getRange(numQubits), retains);

            auto apiFunc = [&](Qureg qureg) { return calcReducedDensityMatrix(qureg, retains.data(), numRetained); };
            auto refFunc = [&](qmatrix ref) { return getPartialTrace(ref, targets); };
            
            TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), getRefDensmatr(), apiFunc, refFunc);
        }
    }

    // @todo input validation
}


/** @} (end defgroup) */



/*
 * @todo
 * UNTESTED FUNCTIONS BELOW
 */

// these require we deploy input objects (Qureg,FullStateDiagMatr) differently
// (as is respectively permitted) to thoroughly test all QuEST control flows

qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);

qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);

qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);

qcomp calcExpecNonHermitianFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


qreal calcFidelity(Qureg qureg, Qureg other);

qreal calcDistance(Qureg qureg1, Qureg qureg2);

qcomp calcInnerProduct(Qureg qureg1, Qureg qureg2);
