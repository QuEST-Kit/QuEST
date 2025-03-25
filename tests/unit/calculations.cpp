/** @file
 * Unit tests of the calculations module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitcalcs Calculations
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
#include <algorithm>
#include <type_traits>

using std::vector;
using namespace Catch::Matchers;



/*
 * UTILITIES
 */


#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[calculations]"


void TEST_ON_CACHED_QUREGS(quregCache cache, auto refState, auto apiFunc, auto refFunc) {

    for (auto& [label, qureg]: cache) {

        DYNAMIC_SECTION( label ) {

            // we initialise to a valid pure or mixed state, rather than the 
            // debug state, since many API functions validate the output domain
            // which is often invalid for unnormalised states
            setToRandomState(refState);
            setQuregToReference(qureg, refState);

            /// @todo 
            /// validate that above is correct, i.e. not the all-zero
            /// state, which would permit flawed tests to succeed

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


void TEST_ON_MIXED_CACHED_QUREGS(quregCache cacheA, quregCache cacheB, auto refA, auto refB, auto apiFunc, auto refFunc) {

    // test all combinations of deployments (where cacheA != cacheB)

    for (auto& [labelA, quregA]: cacheA) {
        for (auto& [labelB, quregB]: cacheB) {

            // skip illegal (local densitymatrix, distributed statevector) combo
            if (quregA.isDensityMatrix != quregB.isDensityMatrix &&  // (sv,dm) or (dm,sv)
                quregA.isDistributed   != quregB.isDistributed   &&  // differing distributions
                quregA.isDistributed   == quregB.isDensityMatrix)    // dm.dist=0
                continue;

            // skip illegal same-size different-GPU-accel combo
            if (quregA.isDensityMatrix  == quregB.isDensityMatrix && // (sv,sv) or (dm,dm)
                quregA.isGpuAccelerated != quregB.isGpuAccelerated)  // differing GPU-accel
                continue;

            DYNAMIC_SECTION( labelA + LABEL_DELIMITER + labelB ) {

                // initialise to random states (rather than debug) since
                // tested function will likely validate output domain
                setToRandomState(refA); setQuregToReference(quregA, refA);
                setToRandomState(refB); setQuregToReference(quregB, refB);

                // results will always be a scalar
                auto apiResult = apiFunc(quregA, quregB);
                auto refResult = refFunc(refA, refB);
                REQUIRE_AGREE( apiResult, refResult );
            }
        }
    }
}


void TEST_ON_CACHED_QUREG_AND_MATRIX(quregCache quregs, matrixCache matrices, auto apiFunc, auto refState, auto refMatr, auto refFunc) {

    // test all combinations of deployments (where cacheA != cacheB)

    for (auto& [labelA, qureg]: quregs) {
        for (auto& [labelB, matrix]: matrices) {

            // all combinations are legal; none are skipped!

            DYNAMIC_SECTION( labelA + LABEL_DELIMITER + labelB ) {

                // set API matrix to pre-initialised reference matrix
                setFullStateDiagMatr(matrix, 0, getDiagonals(refMatr));

                // initialise to random states (rather than debug) since
                // tested function will likely validate output domain
                setToRandomState(refState);
                setQuregToReference(qureg, refState);

                // results will always be a scalar
                auto apiResult = apiFunc(qureg, matrix);
                auto refResult = refFunc(refState, refMatr);
                REQUIRE_AGREE( apiResult, refResult );
            }
        }
    }
}



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
        auto targets = GENERATE_TARGS(numQubits, numTargs);

        PauliStr str = getRandomPauliStr(targets);

        TEST_ALL_QUREGS(
            qureg, calcExpecPauliStr(qureg, str),  
            state, std::real(getReferenceExpectationValue(state, str))
        );
    }

    /// @todo input validation
}



TEST_CASE( "calcExpecPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTerms = GENERATE_COPY( 1, numQubits, getPow2(2*numQubits) );

        GENERATE( range(0,10) );
        PauliStrSum sum = createRandomPauliStrSum(numQubits, numTerms);

        TEST_ALL_QUREGS(
            qureg, calcExpecPauliStrSum(qureg, sum),  
            state, std::real(getReferenceExpectationValue(state, sum))
        );

        destroyPauliStrSum(sum);
    }

    /// @todo input validation
}



TEST_CASE( "calcExpecNonHermitianPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        int numQubits = getNumCachedQubits();
        int numTerms = GENERATE_COPY( 1, numQubits, getPow2(2*numQubits) );

        PauliStrSum sum = createRandomNonHermitianPauliStrSum(numQubits, numTerms);

        TEST_ALL_QUREGS(
            qureg, calcExpecNonHermitianPauliStrSum(qureg, sum),  
            state, getReferenceExpectationValue(state, sum)
        );

        destroyPauliStrSum(sum);
    }

    /// @todo input validation
}



TEST_CASE( "calcProbOfBasisState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qindex index = GENERATE( range(0, 1 << getNumCachedQubits()) );

        TEST_ALL_QUREGS(
            qureg, calcProbOfBasisState(qureg, index),  
            state, getReferenceProbability(state, index)
        );
    }

    /// @todo input validation
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

    /// @todo input validation
}



TEST_CASE( "calcProbOfMultiQubitOutcome", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTargs = GENERATE_COPY( range(1, numQubits+1) );
        auto targets = GENERATE_TARGS(numQubits, numTargs);
        auto outcomes = getRandomOutcomes(numTargs);

        TEST_ALL_QUREGS(
            qureg, calcProbOfMultiQubitOutcome(qureg, targets.data(), outcomes.data(), numTargs),  
            state, getReferenceProbability(state, targets, outcomes)
        );
    }

    /// @todo input validation
}



TEST_CASE( "calcProbsOfAllMultiQubitOutcomes", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_TARGS(numQubits, numTargs);

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

    /// @todo input validation
}



TEST_CASE( "calcTotalProb", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,100) );

        TEST_ALL_QUREGS(
            qureg, calcTotalProb(qureg),  
            state, getReferenceProbability(state)
        );
    }

    /// @todo input validation
}



TEST_CASE( "calcPurity", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,100) );

        TEST_ALL_QUREGS(
            qureg, calcPurity(qureg),  
            state, getReferencePurity(state)
        );
    }

    /// @todo input validation
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

    // no statevector equivalent

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( LABEL_DENSMATR ) {

            int numQubits = getNumCachedQubits();
            int maxNumTargs = getMaxNumTracedQubits(numQubits);
            int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            auto targets = GENERATE_TARGS(numQubits, numTargs);

            auto apiFunc = [&](Qureg qureg) { return calcPartialTrace(qureg, targets.data(), numTargs); };
            auto refFunc = [&](qmatrix ref) { return getPartialTrace(ref, targets); };
            
            TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), getRefDensmatr(), apiFunc, refFunc);
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcReducedDensityMatrix", TEST_CATEGORY ) {

    // no statevector equivalent

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( LABEL_DENSMATR ) {

            int numQubits = getNumCachedQubits();
            int maxNumTraced = getMaxNumTracedQubits(numQubits);
            int minNumRetained = numQubits - maxNumTraced;
            int maxNumRetained = numQubits - 1;
            int numRetained = GENERATE_COPY( range(minNumRetained, maxNumRetained+1) );
            auto retains = GENERATE_TARGS(numQubits, numRetained);
            auto targets = getComplement(getRange(numQubits), retains);

            auto apiFunc = [&](Qureg qureg) { return calcReducedDensityMatrix(qureg, retains.data(), numRetained); };
            auto refFunc = [&](qmatrix ref) { return getPartialTrace(ref, targets); };
            
            TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), getRefDensmatr(), apiFunc, refFunc);
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcInnerProduct", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qvector refSV = getRefStatevec();
        qmatrix refDM = getRefDensmatr();
        auto apiFunc = calcInnerProduct;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC LABEL_DELIMITER LABEL_STATEVEC ) {

            // <a|b>
            auto refFunc = [&](qvector a, qvector b) { return getInnerProduct (a,b); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedStatevecs(), getAltCachedStatevecs(), refSV, refSV, apiFunc, refFunc );
        }

        SECTION( LABEL_STATEVEC LABEL_DELIMITER LABEL_DENSMATR ) {

            // <A|B|A>
            auto refFunc = [&](qvector a, qmatrix b) { return getInnerProduct(a, b * a); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedStatevecs(), getAltCachedDensmatrs(), refSV, refDM, apiFunc, refFunc );
        }

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_STATEVEC ) {

            // <B|A^dagger|B>
            
            auto refFunc = [&](qmatrix a, qvector b) { return getInnerProduct(b, getConjugateTranspose(a) * b); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedDensmatrs(), getAltCachedStatevecs(), refDM, refSV, apiFunc, refFunc );
        }

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_DENSMATR ) {

            // Tr(a^dagger b)
            auto refFunc = [&](qmatrix a, qmatrix b) { return getTrace(getConjugateTranspose(a) * b); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedDensmatrs(), getAltCachedDensmatrs(), refDM, refDM, apiFunc, refFunc );
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcFidelity", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qvector refSV = getRefStatevec();
        qmatrix refDM = getRefDensmatr();
        auto apiFunc = calcFidelity;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC LABEL_DELIMITER LABEL_STATEVEC ) {

            // |<a|b>|^2
            auto refFunc = [&](qvector a, qvector b) { return std::norm(getInnerProduct(a,b)); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedStatevecs(), getAltCachedStatevecs(), refSV, refSV, apiFunc, refFunc );
        }

        SECTION( LABEL_STATEVEC LABEL_DELIMITER LABEL_DENSMATR ) {

            // Re[<A|B|A>]
            auto refFunc = [&](qvector a, qmatrix b) { return std::real(getInnerProduct(a, b * a)); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedStatevecs(), getAltCachedDensmatrs(), refSV, refDM, apiFunc, refFunc );
        }

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_STATEVEC ) {

            // Re[<B|A|B>]
            auto refFunc = [&](qmatrix a, qvector b) { return std::real(getInnerProduct(b, a * b)); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedDensmatrs(), getAltCachedStatevecs(), refDM, refSV, apiFunc, refFunc );
        }

        // (densitymatrix, densitymatrix) not supported
    }

    /// @todo input validation
}



TEST_CASE( "calcDistance", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qvector refSV = getRefStatevec();
        qmatrix refDM = getRefDensmatr();
        auto apiFunc = calcDistance;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC LABEL_DELIMITER LABEL_STATEVEC ) {

            // sqrt(2 - 2 |<A|B>|)
            auto refFunc = [&](qvector a, qvector b) { return std::sqrt(2 - 2 * std::abs(getInnerProduct(a,b))); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedStatevecs(), getAltCachedStatevecs(), refSV, refSV, apiFunc, refFunc);
        }

        SECTION( LABEL_STATEVEC LABEL_DELIMITER LABEL_DENSMATR ) {

            // sqrt(1 - <psi|rho|psi>)
            auto refFunc = [&](qvector a, qmatrix b) { return std::sqrt(1 - std::real(getInnerProduct(a, b * a))); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedStatevecs(), getAltCachedDensmatrs(), refSV, refDM, apiFunc, refFunc );
        }

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_STATEVEC ) {

            // sqrt(1 - <psi|rho|psi>)
            auto refFunc = [&](qmatrix a, qvector b) { return std::sqrt(1 - std::real(getInnerProduct(b, a * b))); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedDensmatrs(), getAltCachedStatevecs(), refDM, refSV, apiFunc, refFunc );
        }

        SECTION( LABEL_DENSMATR LABEL_DELIMITER LABEL_DENSMATR ) {

            // sqrt(Tr((A-B)(A-B)^dagger)
            auto refFunc = [&](qmatrix a, qmatrix b) { return std::sqrt(std::real(getTrace((a-b)*getConjugateTranspose(a-b)))); };

            TEST_ON_MIXED_CACHED_QUREGS( getCachedDensmatrs(), getAltCachedDensmatrs(), refDM, refDM, apiFunc, refFunc );
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcExpecFullStateDiagMatr", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalHermitian(getNumCachedQubits());
        auto apiFunc = calcExpecFullStateDiagMatr;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector state, qmatrix matr) { return getReferenceExpectationValue(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedStatevecs(), getCachedFullStateDiagMatrs(), apiFunc, getRefStatevec(), refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix state, qmatrix matr) { return getReferenceExpectationValue(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedDensmatrs(), getCachedFullStateDiagMatrs(), apiFunc, getRefDensmatr(), refMatr, refFunc);
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcExpecNonHermitianFullStateDiagMatr", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalMatrix(getPow2(getNumCachedQubits()));
        auto apiFunc = calcExpecNonHermitianFullStateDiagMatr;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector state, qmatrix matr) { return getReferenceExpectationValue(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedStatevecs(), getCachedFullStateDiagMatrs(), apiFunc, getRefStatevec(), refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix state, qmatrix matr) { return getReferenceExpectationValue(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedDensmatrs(), getCachedFullStateDiagMatrs(), apiFunc, getRefDensmatr(), refMatr, refFunc);
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcExpecFullStateDiagMatrPower", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalHermitian(getNumCachedQubits());

        // integer and non-integer exponents have different requirements
        bool isIntegerExp = GENERATE( true, false );
        qreal exponent = (isIntegerExp)?
            getRandomInt(-3, 3):
            getRandomReal(-3, 3);
        
        // when exponent is non-integer, the matrix cannot contain negative 
        // numbers which otherwise become complex, making the matrix non-
        // Hermitian and triggering validation
        if (!isIntegerExp)
            for (size_t i=0; i<refMatr.size(); i++)
                refMatr[i][i] = std::abs(std::real(refMatr[i][i]));

        // when exponent is negative, the matrix cannot contain near-zero
        // magnitude numbers which create divergences and trigger validation.
        // Negative exponents are particularly unstable (they can produce 
        // very large matrix elements which ergo have less post-decimal
        // precision and sum catastrophically), so we scale up matrix.
        if (exponent < 0)
            for (size_t i=0; i<refMatr.size(); i++)
                refMatr[i][i] *= 100;

        auto apiFunc = [&](Qureg qureg, FullStateDiagMatr matr) { 
            return calcExpecFullStateDiagMatrPower(qureg, matr, exponent);
        };

        CAPTURE( exponent );

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, qcomp(exponent,0));
                return getReferenceExpectationValue(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedStatevecs(), getCachedFullStateDiagMatrs(), apiFunc, getRefStatevec(), refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, qcomp(exponent,0));
                return getReferenceExpectationValue(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedDensmatrs(), getCachedFullStateDiagMatrs(), apiFunc, getRefDensmatr(), refMatr, refFunc);
        }
    }

    /// @todo input validation
}



TEST_CASE( "calcExpecNonHermitianFullStateDiagMatrPower", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalMatrix(getPow2(getNumCachedQubits()));
        qcomp exponent = getRandomComplex();

        auto apiFunc = [&](Qureg qureg, FullStateDiagMatr matr) { 
            return calcExpecNonHermitianFullStateDiagMatrPower(qureg, matr, exponent);
        };

        CAPTURE( exponent );

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, exponent);
                return getReferenceExpectationValue(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedStatevecs(), getCachedFullStateDiagMatrs(), apiFunc, getRefStatevec(), refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, exponent);
                return getReferenceExpectationValue(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( getCachedDensmatrs(), getCachedFullStateDiagMatrs(), apiFunc, getRefDensmatr(), refMatr, refFunc);
        }
    }

    /// @todo input validation
}



/** @} (end defgroup) */
