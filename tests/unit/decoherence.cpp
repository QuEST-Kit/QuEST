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
    qmatrix reference = getRefDensmatr();

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


void TEST_ON_MIXED_CACHED_QUREGS(auto altQuregCache, auto apiFunc, auto refAlt, auto refFunc) {

    // test all combinations of deployments (where cacheA != cacheB)

    for (auto& [labelA, quregOut]: getCachedDensmatrs()) {
        for (auto& [labelB, quregAlt]: altQuregCache) {

            // skip illegal (local densitymatrix, distributed statevector) combo
            if (!quregOut.isDistributed && quregAlt.isDistributed && !quregAlt.isDensityMatrix)
                continue;

            // skip illegal (local densitymatrix, distributed densitymatrix) combo
            if (quregAlt.isDensityMatrix && quregOut.isDistributed != quregAlt.isDistributed)
                continue;

            DYNAMIC_SECTION( labelA + LABEL_DELIMITER + labelB ) {

                // randomise the output density matrix (both qureg and reference)
                qmatrix refOut = getRandomDensityMatrix(getNumCachedQubits());
                setQuregToReference(quregOut, refOut);

                // randomise the alternate state (may be statevector or density matrix)
                setToRandomState(refAlt);
                setQuregToReference(quregAlt, refAlt);

                apiFunc(quregOut, quregAlt);
                refFunc(refOut, refAlt);
                REQUIRE_AGREE( quregOut, refOut );
            }
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
            std::sqrt(1-prob) * getPauliMatrix(0),
            std::sqrt(prob)   * getPauliMatrix(3)
        };

        auto apiFunc = [&](Qureg qureg) { mixDephasing(qureg, targ, prob); };

        CAPTURE( targ, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, {targ}, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixDepolarising", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int targ = GENERATE_COPY( range(0,numQubits) );
        qreal prob = getRandomReal(0, 3/4.);

        vector<qmatrix> kraus = {
            std::sqrt(1-prob) * getPauliMatrix(0),
            std::sqrt(prob/3) * getPauliMatrix(1), 
            std::sqrt(prob/3) * getPauliMatrix(2), 
            std::sqrt(prob/3) * getPauliMatrix(3),
        };

        auto apiFunc = [&](Qureg qureg) { mixDepolarising(qureg, targ, prob); };

        CAPTURE( targ, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, {targ}, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixDamping", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int targ = GENERATE_COPY( range(0,numQubits) );
        qreal prob = getRandomReal(0, 1);

        vector<qmatrix> kraus = {
            {{1,0},{0,std::sqrt(1-prob)}},
            {{0,std::sqrt(prob)}, {0,0}}
        };

        auto apiFunc = [&](Qureg qureg) { mixDamping(qureg, targ, prob); };

        CAPTURE( targ, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, {targ}, kraus); }
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
            std::sqrt(pI) * getPauliMatrix(0),
            std::sqrt(pX) * getPauliMatrix(1),
            std::sqrt(pY) * getPauliMatrix(2),
            std::sqrt(pZ) * getPauliMatrix(3)
        };

        auto apiFunc = [&](Qureg qureg) { mixPaulis(qureg, targ, pX, pY, pZ); };

        CAPTURE( targ, pX, pY, pZ );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, {targ}, kraus); }
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
            std::sqrt(1-prob) * getKroneckerProduct(i, i),
            std::sqrt(prob/3) * getKroneckerProduct(i, z),
            std::sqrt(prob/3) * getKroneckerProduct(z, i),
            std::sqrt(prob/3) * getKroneckerProduct(z, z)
        };

        auto apiFunc = [&](Qureg qureg) { mixTwoQubitDephasing(qureg, targs[0], targs[1], prob); };

        CAPTURE( targs, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, targs, kraus); }
    }

    /// @todo input validation
}


TEST_CASE( "mixTwoQubitDepolarising", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        auto targs = GENERATE_TARGS( getNumCachedQubits(), 2 );
        qreal prob = getRandomReal(0, 15/16.);

        vector<qmatrix> kraus = { std::sqrt(1-16*prob/15) * getIdentityMatrix(4) };
        for (int a=0; a<4; a++)
            for (int b=0; b<4; b++)
                kraus.push_back( std::sqrt(prob/15) * 
                    getKroneckerProduct(getPauliMatrix(a), getPauliMatrix(b)));

        auto apiFunc = [&](Qureg qureg) { mixTwoQubitDepolarising(qureg, targs[0], targs[1], prob); };

        CAPTURE( targs, prob );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, targs, kraus); }
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
        auto apiFunc = [&](Qureg qureg) { mixKrausMap(qureg, targs.data(), numTargs, map); };

        CAPTURE( maxNumTargs, targs, numKraus );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, targs, matrices); }

        destroyKrausMap(map);
    }

    /// @todo input validation
}


TEST_CASE( "mixSuperOp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int maxFlag = TEST_MAX_NUM_SUPEROP_TARGETS;
        int maxNumTargs = (maxFlag != 0 && numQubits > maxFlag)?
            maxFlag : numQubits;

        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
        auto targs = GENERATE_TARGS( numQubits, numTargs );
        auto matrices = getRandomKrausMap(numTargs, getRandomInt(1,4+1));

        SuperOp superOp = createSuperOp(numTargs);
        setSuperOp(superOp, getSuperOperator(matrices));
        auto apiFunc = [&](Qureg qureg) { mixSuperOp(qureg, targs.data(), numTargs, superOp); };

        CAPTURE( targs );
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(apiFunc, targs, matrices); }

        destroySuperOp(superOp);
    }

    /// @todo input validation
}


TEST_CASE( "mixQureg", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qreal prob = getRandomReal(0, 1);
        auto apiFunc = [&](Qureg a, Qureg b) { mixQureg(a, b, prob); };

        CAPTURE( prob );
        
        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_STATEVEC ) { 

            auto refFunc = [&](qmatrix& a, qvector b) { a = (1-prob)*a + prob*getOuterProduct(b,b); };
            
            TEST_ON_MIXED_CACHED_QUREGS( getAltCachedStatevecs(), apiFunc, getRefStatevec(), refFunc); 
        }

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_DENSMATR ) {

            auto refFunc = [&](qmatrix& a, qmatrix b) { a = (1-prob)*a + prob*b; };

            TEST_ON_MIXED_CACHED_QUREGS( getAltCachedDensmatrs(), apiFunc, getRefDensmatr(), refFunc); 
        }
    }

    /// @todo input validation
}


/** @} (end defgroup) */
