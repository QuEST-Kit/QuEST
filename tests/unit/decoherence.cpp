/** @file
 * Unit tests of the decoherence module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitdeco Decoherence
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"

#include <vector>
#include <algorithm>

using std::vector;



/*
 * UTILITIES
 */


#define TEST_CATEGORY "[unit][decoherence]"


void TEST_ON_CACHED_QUREGS(auto apiFunc, vector<int> targs, vector<qmatrix> kraus) {

    // all tests use a fixed-size density matrix
    qmatrix reference = getZeroMatrix(getPow2(getNumCachedQubits()));

    for (auto& [label, qureg]: getCachedDensmatrs()) {

        DYNAMIC_SECTION( label ) {

            // no need to validate whether qureg successfully
            // enters the debug state here, because the below
            // serial setToDebugState() is guaranteed to succeed
            initDebugState(qureg);
            setToDebugState(reference);

            apiFunc(qureg);
            applyReferenceOperator(reference, targs, kraus);

            REQUIRE_AGREE( qureg, reference );
        }
    }
}



/** 
 * TESTS
 * 
 * @ingroup unitdeco
 * @{
 */


TEST_CASE( "mixDephasing", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int targ = GENERATE_COPY( range(0,numQubits) );
        qreal prob = getRandomReal(0, 1/2.);

        vector<qmatrix> kraus = {
            sqrt(1-prob) * getPauliMatrix(0),
            sqrt(prob)   * getPauliMatrix(3)
        };

        auto func = [&](Qureg qureg) { mixDephasing(qureg, targ, prob); };

        CAPTURE( targ, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, {targ}, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixDepolarising", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int targ = GENERATE_COPY( range(0,numQubits) );
        qreal prob = getRandomReal(0, 3/4.);

        vector<qmatrix> kraus = {
            sqrt(1-prob) * getPauliMatrix(0),
            sqrt(prob/3) * getPauliMatrix(1), 
            sqrt(prob/3) * getPauliMatrix(2), 
            sqrt(prob/3) * getPauliMatrix(3),
        };

        auto func = [&](Qureg qureg) { mixDepolarising(qureg, targ, prob); };

        CAPTURE( targ, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, {targ}, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixDamping", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int targ = GENERATE_COPY( range(0,numQubits) );
        qreal prob = getRandomReal(0, 1);

        vector<qmatrix> kraus = {
            {{1,0},{0,sqrt(1-prob)}},
            {{0,sqrt(prob)}, {0,0}}
        };

        auto func = [&](Qureg qureg) { mixDamping(qureg, targ, prob); };

        CAPTURE( targ, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, {targ}, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixPaulis", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int targ = GENERATE_COPY( range(0,numQubits) );

        qreal pX = getRandomReal(0, 1);
        qreal pY = getRandomReal(0, 1);
        qreal pZ = getRandomReal(0, 1);

        // we require pX+pY+pZ <= 1
        qreal norm = pX + pY + pZ;
        pX /= norm;
        pY /= norm;
        pZ /= norm;

        // and max(pX,pY,pZ) <= 1-pX-pY-pZ, which we'll
        // lazily achieve with iteration (truly stinky)
        qreal pI = 1 - pX - pY - pZ;
        while (std::max({pX,pY,pZ}) > pI) {
            pX /= 1.1;
            pY /= 1.1;
            pZ /= 1.1;
            pI = 1 - pX - pY - pZ;
        }

        vector<qmatrix> kraus = {
            sqrt(pI) * getPauliMatrix(0),
            sqrt(pX) * getPauliMatrix(1),
            sqrt(pY) * getPauliMatrix(2),
            sqrt(pZ) * getPauliMatrix(3)
        };

        auto func = [&](Qureg qureg) { mixPaulis(qureg, targ, pX, pY, pZ); };

        CAPTURE( targ, pX, pY, pZ );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, {targ}, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixTwoQubitDephasing", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        auto targs = GENERATE_TARGS( getNumCachedQubits(), 2 );
        qreal prob = getRandomReal(0, 3/4.);

        qmatrix i = getPauliMatrix(0);
        qmatrix z = getPauliMatrix(3);

        vector<qmatrix> kraus = {
            sqrt(1-prob) * getKroneckerProduct(i, i),
            sqrt(prob/3) * getKroneckerProduct(i, z),
            sqrt(prob/3) * getKroneckerProduct(z, i),
            sqrt(prob/3) * getKroneckerProduct(z, z)
        };

        auto func = [&](Qureg qureg) { mixTwoQubitDephasing(qureg, targs[0], targs[1], prob); };

        CAPTURE( targs, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, targs, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixTwoQubitDepolarising", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        auto targs = GENERATE_TARGS( getNumCachedQubits(), 2 );
        qreal prob = getRandomReal(0, 15/16.);

        vector<qmatrix> kraus = { sqrt(1-16*prob/15) * getIdentityMatrix(4) };
        for (int a=0; a<4; a++)
            for (int b=0; b<4; b++)
                kraus.push_back( sqrt(prob/15) * 
                    getKroneckerProduct(getPauliMatrix(a), getPauliMatrix(b)));

        auto func = [&](Qureg qureg) { mixTwoQubitDepolarising(qureg, targs[0], targs[1], prob); };

        CAPTURE( targs, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, targs, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixKrausMap", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int maxFlag = TEST_MAX_NUM_SUPEROP_TARGETS;
        int numQubits = getNumCachedQubits();
        int maxNumTargs = (maxFlag != 0 && numQubits > maxFlag)?
            maxFlag : numQubits;

        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
        int numKraus = GENERATE( 1, 2, 10 );
        auto targs = GENERATE_TARGS( numQubits, numTargs );
        auto matrices = getRandomKrausMap(numTargs, numKraus);

        KrausMap map = createKrausMap(numTargs, numKraus);
        setKrausMap(map, matrices);
        auto func = [&](Qureg qureg) { mixKrausMap(qureg, targs.data(), numTargs, map); };

        CAPTURE( maxNumTargs, targs, numKraus );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(func, targs, matrices); }

        destroyKrausMap(map);
    }

    /// @todo input validation
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

// these require we deploy the Quregs differently
// to thoroughly test all QuEST control flows

void mixQureg(Qureg qureg, Qureg other, qreal prob);