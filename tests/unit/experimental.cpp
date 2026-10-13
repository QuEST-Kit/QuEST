/** @file
 * Unit tests of the environment module.
 *
 * @author Oliver Brown
 * @author Tyson Jones
 * @author Ashmit JaiSarita Gupta (checkpoint test prototype)
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
#include "tests/utils/compare.hpp"

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

            int badNumTPB = 102400; // exceeds expected 1024 max

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

        // Single process deletes checkpoint file (assumes a shared filesystem; if not, who cares about the scraps?)
        // Note these syncs are ESSENTIAL for correct behaviour, else root can begin deletion while a subsequent node
        // proceeds to the below validation and re-creates some files within the same direc, causing MPI hangs. Ouch!
        syncQuESTEnv();
        if (getQuESTEnv().rank == 0)
            std::filesystem::remove_all(outFn);
        syncQuESTEnv();
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedStatevec();

        SECTION( "adios2 not compiled" ) {

            if (!QUEST_COMPILE_ADIOS2)
                REQUIRE_THROWS_WITH( saveQuregToFile(qureg, "dummy.bp"), ContainsSubstring("compiled with ADIOS2") );

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
                // surprisingly hard to find cross-OS illegal names!
                #if defined(_MSC_VER)
                    auto badFn = GENERATE( ":", "?", "*" );
                #else
                    auto badFn = GENERATE( "", "\0" );
                #endif
                REQUIRE_THROWS_WITH( saveQuregToFile(qureg, badFn), ContainsSubstring("could not be opened") );
            }

            SUCCEED( );
        }
    }
}


TEST_CASE( "createQuregFromFile", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {
        
        const char* checkpointFn = "test_checkpoint.bp";

        // We will iterate the cached Quregs so the save path is exercised under every
        // deployment combination (serial, OMP, MPI, GPU and their mixtures). However,
        // the restored Qureg uses a distribution chosen by the auto-deployer, which is
        // not permitted to differ from the checkpointed distribution. We know, given
        // the unit test Quregs are so small, that distribution is NEVER automatically
        // enabled; so we will forbid testing with distributed Quregs
        int legalNumNodes = 1;

        auto testFunc = [&](Qureg qureg) {

            initRandomPureState(qureg);
            REQUIRE_NOTHROW( saveQuregToFile(qureg, checkpointFn) );

            // skip restoration when new Qureg distribution would disagree with old
            if (qureg.numNodes != legalNumNodes)
                return;

            Qureg newQureg = createQuregFromFile(checkpointFn);
            REQUIRE_AGREE(qureg, newQureg);

            destroyQureg(newQureg);
        };

        // skip correctness tests if ADIOS2 not compiled
        SECTION( LABEL_STATEVEC ) { if (QUEST_COMPILE_ADIOS2) TEST_ON_CACHED_QUREGS(getCachedStatevecs(), testFunc); SUCCEED( ); }
        SECTION( LABEL_DENSMATR ) { if (QUEST_COMPILE_ADIOS2) TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), testFunc); SUCCEED( ); }

        CAPTURE( checkpointFn );

        // Single process deletes checkpoint file (assumes a shared filesystem; if not, who cares about the scraps?).
        // Note these syncs are ESSENTIAL for correct behaviour, else root can begin deletion while a subsequent node
        // proceeds to the below validation and re-creates some files within the same direc, causing MPI hangs. Ouch!
        syncQuESTEnv();
        if (getQuESTEnv().rank == 0)
            std::filesystem::remove_all(checkpointFn);
        syncQuESTEnv();
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "adios2 not compiled" ) {

            if (!QUEST_COMPILE_ADIOS2)
                REQUIRE_THROWS_WITH( createQuregFromFile("dummy.bp"), ContainsSubstring("compiled with ADIOS2") );

            SUCCEED( );
        }

        SECTION( "bad name" ) {

            if (QUEST_COMPILE_ADIOS2)
                REQUIRE_THROWS_WITH( createQuregFromFile("BAD_FILENAME"), ContainsSubstring("could not be opened") );

            SUCCEED( );
        }

        SECTION( "differing distributions" ) {

            // Distributions can only differ when QuEST is distributed over more than 1 node
            if (QUEST_COMPILE_ADIOS2 && getQuESTEnv().numNodes > 1) {

                // Create a new distributed qureg; we know createQuregFromFile() will create
                // non-distributed, since unit-test-size Quregs auto-deploy to non-distributed
                Qureg quregDistrib = createCustomQureg(getNumCachedQubits(), 0, /*useDistrib=*/1, 0, 0);

                CAPTURE( quregDistrib.numNodes );

                // Write qureg to file, then deliberately fail to restore it
                const char* fn = "test_checkpoint.bp";
                saveQuregToFile(quregDistrib, fn);
                REQUIRE_THROWS_WITH( createQuregFromFile(fn), ContainsSubstring("distributions must match") );

                // cleanup
                destroyQureg(quregDistrib);
                syncQuESTEnv();
                if (getQuESTEnv().rank == 0)
                    std::filesystem::remove_all(fn);
                syncQuESTEnv();
            }

            SUCCEED( );
        }

        // We do not presently test the below validations, since it will require
        // externally generating and saving ADIOS2 files; quite a pain!
        // SECTION( "differing precision" ) { }
        // SECTION( "overflow" ) { }
        // SECTION( "insufficient RAM" ) { }
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

// nothing! :^)
