/** @file
 * Unit tests of the trotterisation module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unittrotter Trotterisation
 * @ingroup unittests
 */

#include "quest.h"

#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/cache.hpp"
#include "tests/utils/random.hpp"



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[trotterisation]"


void TEST_ON_CACHED_QUREGS(quregCache quregs, auto& refFunc, auto& regularFunc, auto& randFunc) {

    for (auto& [label, refQureg]: quregs) {

        DYNAMIC_SECTION( label ) {
            initDebugState(refQureg);

            Qureg regularQureg = createCloneQureg(refQureg);
            Qureg randQureg = createCloneQureg(refQureg);

            refFunc(refQureg);
            regularFunc(regularQureg);
            randFunc(randQureg);

            double regularDistance = calcDistance(regularQureg, refQureg);
            double randDistance = calcDistance(randQureg, refQureg);

            REQUIRE( randDistance < regularDistance );

            destroyQureg(regularQureg);
            destroyQureg(randQureg);
        }
    }
}

/**
 * @todo
 * Basic validation for randomisation, should be expanded and merged
 * once the Trotterisation function tests have been implemented.
 */

TEST_CASE( "randomisedTrotter", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTerms = 25;
        int reps = 50;
        double time = 1.0;

        int refOrder = 4;
        int order = GENERATE_COPY(1, 2);

        GENERATE( range(0, 10) );
        PauliStrSum sum = createRandomPauliStrSum(numQubits, numTerms);

        auto refFunc = [&](Qureg qureg) { applyTrotterizedUnitaryTimeEvolution(qureg, sum, time, refOrder, reps, false); };
        auto regularFunc = [&](Qureg qureg) { applyTrotterizedUnitaryTimeEvolution(qureg, sum, time, order, reps, false); };
        auto randFunc = [&](Qureg qureg) { applyTrotterizedUnitaryTimeEvolution(qureg, sum, time, order, reps, true); };
        
        TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), refFunc, regularFunc, randFunc);
        TEST_ON_CACHED_QUREGS(getCachedStatevecs(), refFunc, regularFunc, randFunc);

        destroyPauliStrSum(sum);

    }

}


/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void applyTrotterizedNonUnitaryPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps, bool permutePaulis);

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedControlledPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedUnitaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal time, int order, int reps, bool permutePaulis);

void applyTrotterizedImaginaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal tau, int order, int reps, bool permutePaulis);

void applyTrotterizedNoisyTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal* damps, PauliStr* jumps, int numJumps, qreal time, int order, int reps, bool permutePaulis);
