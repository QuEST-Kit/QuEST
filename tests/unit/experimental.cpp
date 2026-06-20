/** @file
 * Unit tests of the environment module.
 *
 * @author Oliver Brown
 * @author Tyson Jones
 * 
 * @defgroup unitexperi Experimental
 * @ingroup unittests
 */

#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/config.hpp"
#include "tests/utils/cache.hpp"

#include <filesystem>

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */


#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[experimental]"


void TEST_ON_CACHED_QUREGS(quregCache quregs, auto testFunc) {

    for (auto& [label, qureg]: quregs) {

        DYNAMIC_SECTION( label ) {

            testFunc(qureg);
        }
    }
}



/** 
 * TESTS
 * 
 * @ingroup unitexperi
 * @{
 */


TEST_CASE( "setQuESTNumGpuThreadsPerBlock", TEST_CATEGORY ) {

    // remember the default number for later restoration (hence static)
    static int initNumTPB = getQuESTNumGpuThreadsPerBlock();

    SECTION( LABEL_CORRECTNESS ) {

        // begin at 64 (AMD min, larger than NVIDIA min of 32),
        // stop at 1024 (should be less than dev-specific max)
        int inNumTPB = GENERATE( 64, 128, 256, 512, 1024 ); 
        setQuESTNumGpuThreadsPerBlock(inNumTPB);

        int outNumTPB = getQuESTNumGpuThreadsPerBlock();
        REQUIRE( inNumTPB == outNumTPB );
        
        // BEWARE that we do not here test whether all QuEST
        // operators succeed with the various numTBP; that must
        // be ad hoc asssesed via updating the numTBP env-var
        // before launching the entirety of the tests
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "Negative" ) {

            int badNumTPB = GENERATE( 0, -1, -9999 );
            REQUIRE_THROWS_WITH( setQuESTNumGpuThreadsPerBlock(badNumTPB), ContainsSubstring( "must be positive" ) );
        }

        SECTION( "Indivisible by warp size" ) {

            // If HIP status was attached to QuESTEnv, we could do:
            //     QuESTEnv env = getQuESTEnv();
            //     int warpSize = (env.isGpuAccelerated && env.isHipCompiled)? 64 : 32;
            // Since this currently isn't the case, we assume a warp size of 32,
            // which will mean when this test is run on AMD GPUs, the below tested
            // badNumTBP won't be as interestingly/rigorously spread
            int warpSize = 32;

            int badNumTPB = GENERATE_COPY( warpSize - 1, warpSize + 1, warpSize + warpSize/2, 3*warpSize + warpSize/2 );

            REQUIRE_THROWS_WITH( setQuESTNumGpuThreadsPerBlock(badNumTPB), ContainsSubstring( "does not divide evenly into the warp size" ) );
        }

        SECTION( "Exceeds device maximum" ) {

            int badNumTPB = 999999; // exceeds expected 1024 max

            // Cannot be tested (since validation not imposed) when GPU is not actively used
            if (getQuESTEnv().isGpuAccelerated)
                REQUIRE_THROWS_WITH( setQuESTNumGpuThreadsPerBlock(badNumTPB), ContainsSubstring( "Exceeds the hardware-imposed maximum" ) );

            SUCCEED( );
        }
    }

    // restore numTBP, so as not to interfere with other tests
    setQuESTNumGpuThreadsPerBlock(initNumTPB);
}


TEST_CASE( "getQuESTNumGpuThreadsPerBlock", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // check initial value matches either the env-var (if set),
        // or the fixed default in the codebase (hardcoded in test utils)
        int defaultNum = getDefaultNumGpuThreadsPerBlock(); // test util via env-var
        int reportedNum = getQuESTNumGpuThreadsPerBlock();  // QuEST API

        REQUIRE( defaultNum == reportedNum );

        // further testing of this function appears in setQuESTNumGpuThreadsPerBlock()
    }

    SECTION( LABEL_VALIDATION ) {

        // there is none (except untestable env is init!)
        SUCCEED( );
    }
}


TEST_CASE( "saveQuregToFile", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {
        
        const char* outFn = "test_checkpoint.bp";

        auto testFunc = [&](Qureg qureg) {
            initRandomPureState(qureg);
            REQUIRE_NOTHROW( saveQuregToFile(qureg, outFn) );

            // note that we are NOT validating the contents was correct;
            // that will be performed by the createQuregFromFile() test
        };

        // skip correctness tests if ADIOS2 not compiled
        SECTION( LABEL_STATEVEC ) { if (QUEST_COMPILE_ADIOS2) TEST_ON_CACHED_QUREGS(getCachedStatevecs(), testFunc); SUCCEED( ); }
        SECTION( LABEL_DENSMATR ) { if (QUEST_COMPILE_ADIOS2) TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), testFunc); SUCCEED( ); }

        // single process deletes checkpoint file (assumes a shared filesystem; if not, who cares about the scraps?)
        syncQuESTEnv();
        if (getQuESTEnv().rank == 0)
            std::filesystem::remove_all(outFn);
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedStatevec();

        SECTION( "adios2 not compiled" ) {

            if (!QUEST_COMPILE_ADIOS2)
                REQUIRE_THROWS_WITH( saveQuregToFile(qureg, "dummy.bp"), ContainsSubstring("blah") );

            SUCCEED( );
        }

        SECTION( "qureg uninitialised" ) {

            if (QUEST_COMPILE_ADIOS2) {
                Qureg badQureg;
                badQureg.numQubits = -123;
                REQUIRE_THROWS_WITH( saveQuregToFile(badQureg, "dummy.bp"), ContainsSubstring("Received an invalid Qureg") );
            }

            SUCCEED( );
        }

        SECTION( "bad name" ) {

            if (QUEST_COMPILE_ADIOS2) {
                auto badFn = GENERATE( "" ); // surprisingly hard to find cross-OS illegal names!
                REQUIRE_THROWS_WITH( saveQuregToFile(qureg, badFn), ContainsSubstring("could not be opened") );
            }

            SUCCEED( );
        }
    }
}


    // TODO:
    // - fix this guard! Just runtime skip 
    // - fix tests; don't use custom comparison, use existing utils
    // - negative test of when PRECISION CHANGES
    //   (can we invoke a QuEST subprocess to WRITE to file?!?! Probs not )
    // - extend tests to CHANGE DEPLOYMENT of the Qureg pre and post restoration!
    // - note we cannot actually make negative test changes of precision!
    // - separate test into two functions, for each API func
    

#ifdef QUEST_COMPILE_ADIOS2

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <string>

namespace {

    const char* SV_FILE = "test_checkpoint_statevector.bp";
    const char* DM_FILE = "test_checkpoint_densitymatrix.bp";

    qreal maxStatevectorAmpDiff(Qureg a, Qureg b) {
        qreal m = 0;
        for (qindex i = 0; i < a.numAmps; i++)
            m = std::max(m, std::abs(getQuregAmp(a, i) - getQuregAmp(b, i)));
        return m;
    }

    qreal maxDensityMatrixAmpDiff(Qureg a, Qureg b) {
        qreal m = 0;
        qindex dim = (qindex) 1 << a.numQubits;
        for (qindex r = 0; r < dim; r++)
            for (qindex c = 0; c < dim; c++)
                m = std::max(m, std::abs(getDensityQuregAmp(a, r, c) - getDensityQuregAmp(b, r, c)));
        return m;
    }

    // distributed-safe cleanup: a barrier guarantees every node has finished
    // reading the shared file, only rank 0 deletes it (concurrent removal races),
    // and a second barrier stops the next write racing a half-removed directory.
    void removeCheckpointFile(const char* fn) {
        syncQuESTEnv();
        if (getQuESTEnv().rank == 0)
            std::filesystem::remove_all(fn);
        syncQuESTEnv();
    }
}

TEST_CASE( "saveQuregToFile and createQuregFromFile", TEST_CATEGORY ) {

    // TODO / DEBUG / BEWARE!
    // These tests are insufficient! They only ever test createQuregFromFile()
    // (and ergo validate saveQuregToFile() worked properly) for non-distributed
    // Quregs! This is because createQuregFromFile() uses the distribution of the
    // autodeployer, which for our tiny unit-test Quregs, will always default to
    // non-distributed. Distributed Qureg restoration is totally untested!

    SECTION( LABEL_CORRECTNESS ) {

        // We will iterate the cached Quregs so the save path is exercised under every
        // deployment combination (serial, OMP, MPI, GPU and their mixtures). However,
        // the restored Qureg uses a distribution chosen by the auto-deployer, which is
        // not permitted to differ from the checkpointed distribution; we skip those!
        Qureg svDummy = createQureg(getNumCachedQubits());
        Qureg dmDummy = createDensityQureg(getNumCachedQubits());
        int legalSvNumNodes = svDummy.numNodes;
        int legalDmNumNodes = dmDummy.numNodes;
        destroyQureg(svDummy);
        destroyQureg(dmDummy);

        SECTION( LABEL_STATEVEC ) {

            for (auto& [label, q] : getCachedStatevecs()) {
                DYNAMIC_SECTION( label ) {

                    // always test writing succeeds
                    initRandomPureState(q);
                    REQUIRE_NOTHROW( saveQuregToFile(q, SV_FILE) );

                    // skip restoration when new Qureg distribution would disagree with old
                    if (q.numNodes != legalSvNumNodes)
                        continue;

                    Qureg r = createQuregFromFile(SV_FILE);

                    CHECK( r.numQubits       == q.numQubits );
                    CHECK( r.isDensityMatrix == q.isDensityMatrix );
                    CHECK( maxStatevectorAmpDiff(q, r) < 1e-12 );

                    destroyQureg(r);
                    removeCheckpointFile(SV_FILE);
                }
            }
        }

        SECTION( LABEL_DENSMATR ) {

            for (auto& [label, q] : getCachedDensmatrs()) {
                DYNAMIC_SECTION( label ) {

                    // always test writing succeeds
                    initRandomMixedState(q, /*numPureStates=*/10);
                    REQUIRE_NOTHROW( saveQuregToFile(q, DM_FILE) );

                    // skip cached quregs with illegal distributions
                    if (q.numNodes != legalDmNumNodes)
                        continue;

                    Qureg r = createQuregFromFile(DM_FILE);

                    CHECK( r.numQubits       == q.numQubits );
                    CHECK( r.isDensityMatrix == q.isDensityMatrix );
                    CHECK( maxDensityMatrixAmpDiff(q, r) < 1e-12 );

                    destroyQureg(r);
                    removeCheckpointFile(DM_FILE);
                }
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        // The only checkpointing-specific validation - calling the API when QuEST
        // was compiled without checkpointing - is unreachable here, since this
        // file only compiles under QUEST_COMPILE_ADIOS2. ADIOS2's own
        // runtime errors (e.g. a missing file) are not QuEST validation errors.
        SUCCEED( );
    }
}

#endif // QUEST_COMPILE_ADIOS2



/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

// nothing! :^)
