/** @file
 * Unit tests of the initialisations module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitinit Initialisation
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
#include "tests/utils/measure.hpp"
#include "tests/utils/random.hpp"



/*
 * UTILITIES
 */


#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[initialisations]"


void TEST_ON_CACHED_QUREGS(quregCache quregs, auto testFunc) {

    for (auto& [label, qureg]: quregs) {

        DYNAMIC_SECTION( label ) {

            testFunc(qureg);
        }
    }
}


void TEST_ON_CACHED_QUREGS(quregCache quregs, auto apiFunc, auto refState) {

    // assumes refState is already initialised

    auto testFunc = [&](Qureg qureg) {
        apiFunc(qureg);
        REQUIRE_AGREE( qureg, refState );
    };

    TEST_ON_CACHED_QUREGS(quregs, testFunc);
}



/** 
 * TESTS
 * 
 * @ingroup unitinit
 * @{
 */


TEST_CASE( "initBlankState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), initBlankState, getRefStatevec()); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), initBlankState, getRefDensmatr()); }
    }

    /// @todo input validation
}


TEST_CASE( "initZeroState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qvector refVec = getRefStatevec(); refVec[0]    = 1; // |0> = {1, 0...}
        qmatrix refMat = getRefDensmatr(); refMat[0][0] = 1; // |0><0| = {{1,0...},{0...}...}

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), initZeroState, refVec); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), initZeroState, refMat); }
    }

    /// @todo input validation
}


TEST_CASE( "initPlusState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        qreal vecElem = 1. / std::sqrt(getPow2(numQubits));
        qreal matElem = 1. / getPow2(numQubits);

        qvector refVec = getConstantVector(getPow2(numQubits), vecElem); // |+> = 1/sqrt(2^N) {1, ...}
        qmatrix refMat = getConstantMatrix(getPow2(numQubits), matElem); // |+><+| = 1/2^N {{1, ...}, ...}

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), initPlusState, refVec); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), initPlusState, refMat); }
    }

    /// @todo input validation
}


TEST_CASE( "initClassicalState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numInds = (int) getPow2(numQubits);
        int stateInd = GENERATE_COPY( range(0,numInds) );

        qvector refVec = getRefStatevec(); refVec[stateInd]           = 1; // |i> = {0, ..., 1, 0, ...}
        qmatrix refMat = getRefDensmatr(); refMat[stateInd][stateInd] = 1; // |i><i| 

        auto apiFunc = [&](Qureg qureg) { initClassicalState(qureg, stateInd); };

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), apiFunc, refVec); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), apiFunc, refMat); }
    }

    /// @todo input validation
}


TEST_CASE( "initDebugState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qvector refVec = getRefStatevec(); setToDebugState(refVec); // |debug>
        qmatrix refMat = getRefDensmatr(); setToDebugState(refMat); // ||debug>>

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), initDebugState, refVec); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), initDebugState, refMat); }
    }

    /// @todo input validation
}


TEST_CASE( "initRandomPureState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // this test does not use reference states
        GENERATE( range(0,10) );

        auto testFunc = [&](Qureg qureg) {

            initRandomPureState(qureg);
            
            /// @todo
            /// these not-all-same-amp checks can be made much more rigorous,
            /// by e.g. asserting distinct nodes haven't generated all the same
            /// amplitudes (we currently observe this by eye)
            syncQuregFromGpu(qureg);
            REQUIRE( qureg.cpuAmps[0] != qureg.cpuAmps[1] );

            qreal prob = calcTotalProb(qureg);
            REQUIRE_AGREE( prob, 1 );

            qreal purity = calcPurity(qureg);
            REQUIRE_AGREE( purity, 1 );
        };

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "initRandomMixedState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // this test does not use reference states

        GENERATE( range(0,10) );
        int numPureStates = GENERATE( 1, 2, 10 );

        auto testFunc = [&](Qureg qureg) {

            initRandomMixedState(qureg, numPureStates);
            
            /// @todo
            /// these not-all-same-amp checks can be made much more rigorous,
            /// by e.g. asserting distinct nodes haven't generated all the same
            /// amplitudes (we currently observe this by eye)
            syncQuregFromGpu(qureg);
            REQUIRE( qureg.cpuAmps[0] != qureg.cpuAmps[1] ); // performed on all nodes

            qreal prob = calcTotalProb(qureg);
            REQUIRE_AGREE( prob, 1 );

            qreal purity = calcPurity(qureg);
            if (numPureStates == 1)
                REQUIRE_AGREE( purity, 1 );
            else
                REQUIRE( purity < 1 );
        };

        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "initArbitraryPureState", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // works for unnormalised states
        qvector refVec = getRandomVector(getPow2(getNumCachedQubits()));
        qmatrix refMat = getOuterProduct(refVec, refVec);

        auto apiFunc = [&](Qureg qureg) { initArbitraryPureState(qureg, refVec.data()); };

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), apiFunc, refVec); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), apiFunc, refMat); }
    }

    /// @todo input validation
}




TEST_CASE( "setQuregAmps", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numTotalAmps = getPow2(getNumCachedQubits());
        int numSetAmps = GENERATE_COPY( range(0,numTotalAmps+1) ); 
        int startInd = GENERATE_COPY( range(0,numTotalAmps-numSetAmps) );
        qvector amps = getRandomVector(numSetAmps);

        auto testFunc = [&](Qureg qureg) {

            // initialise qureg randomly
            qvector refVec = getRandomVector(numTotalAmps);
            setQuregToReference(qureg, refVec);

            // modify only subset of refVec amps and qureg...
            setSubVector(refVec, amps, startInd);
            setQuregAmps(qureg, startInd, amps.data(), numSetAmps);
        
            // so that we simultaneously check targeted amps 
            // are modified while non-targeted are not
            REQUIRE_AGREE( qureg, refVec );
        };

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "setDensityQuregFlatAmps", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numTotalRows = getPow2(getNumCachedQubits());
        int numTotalAmps = numTotalRows * numTotalRows;

        // systematic iteration is WAY too slow
        GENERATE( range(0,1000) );
        int numSetAmps = getRandomInt(0, numTotalAmps + 1);
        int startInd = getRandomInt(0, numTotalAmps - numSetAmps);
        qvector amps = getRandomVector(numSetAmps);

        auto testFunc = [&](Qureg qureg) {

            // initialise qureg randomly
            qmatrix refMat = getRandomMatrix(numTotalRows);
            setQuregToReference(qureg, refMat);

            // overwrite a contiguous region of row-major refMat, column-wise
            refMat = getTranspose(refMat);
            setSubMatrix(refMat, amps, startInd);
            refMat = getTranspose(refMat);

            // modify the same contiguous region of column-major qureg
            setDensityQuregFlatAmps(qureg, startInd, amps.data(), numSetAmps);
        
            // check that both targeted and non-targeted amps agree
            REQUIRE_AGREE( qureg, refMat );
        };

        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "setDensityQuregAmps", TEST_CATEGORY ) {

    //void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);


    SECTION( LABEL_CORRECTNESS ) {

        int numTotalRowsCols = getPow2(getNumCachedQubits());

        // systematic iteration is WAY too slow
        GENERATE( range(0,1000) );
        int numSetRows = getRandomInt(1, numTotalRowsCols+1);
        int numSetCols = getRandomInt(1, numTotalRowsCols+1);
        int startRow = getRandomInt(0, numTotalRowsCols - numSetRows);
        int startCol = getRandomInt(0, numTotalRowsCols - numSetCols);

        // caution that amps is 'qmatrix' despite not being square
        qmatrix amps = getRandomNonSquareMatrix(numSetRows, numSetCols);

        auto testFunc = [&](Qureg qureg) {

            // initialise qureg randomly
            qmatrix refMat = getRandomMatrix(numTotalRowsCols);
            setQuregToReference(qureg, refMat);

            // API needs nested pointers
            std::vector<qcomp*> rowPtrs(numSetRows);
            for (size_t r=0; r<numSetRows; r++)
                rowPtrs[r] = amps[r].data();

            // overwrite a sub-matrix of refMat and Qureg
            setSubMatrix(refMat, amps, startRow, startCol);
            setDensityQuregAmps(qureg, startRow, startCol, rowPtrs.data(), numSetRows, numSetCols);

            // check that both targeted and non-targeted amps agree
            REQUIRE_AGREE( qureg, refMat );
        };

        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), testFunc); }
    
    }

    /// @todo input validation
}


TEST_CASE( "setQuregToRenormalized", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        qindex dim = getPow2(getNumCachedQubits());
        qvector refVec = getRandomVector(dim);
        qmatrix refMat = getRandomMatrix(dim);

        // eliminate random chance of tr(refMat)=0, triggering validation
        if (doScalarsAgree(getTrace(refMat), 0))
            refMat[0][0] += 1/(qreal) dim;

        // [=] stores current (pre-normalised) reference objects
        auto funcVec = [=](Qureg qureg) { setQuregToReference(qureg, refVec); setQuregToRenormalized(qureg); };
        auto funcMat = [=](Qureg qureg) { setQuregToReference(qureg, refMat); setQuregToRenormalized(qureg); };

        // setQuregToRenormalized() makes statevectors become valid
        refVec = getNormalised(refVec);

        // but it only divides density matrices by the sum of the real-elems of their diagonals
        refMat /= getReferenceProbability(refMat);

        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(getCachedStatevecs(), funcVec, refVec); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), funcMat, refMat); }
    }

    /// @todo input validation
}


TEST_CASE( "setQuregToPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        int numQubits = getNumCachedQubits();
        int numTerms = GENERATE_COPY( 1, numQubits, getPow2(2*numQubits) );
        PauliStrSum sum = createRandomPauliStrSum(numQubits, numTerms);
        qmatrix refMat = getMatrix(sum, numQubits);

        auto apiFunc = [&](Qureg qureg) { setQuregToPauliStrSum(qureg, sum); };

        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), apiFunc, refMat); }
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

void initPureState(Qureg qureg, Qureg pure);

void setQuregToClone(Qureg targetQureg, Qureg copyQureg);

void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2);

void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);

void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits);
