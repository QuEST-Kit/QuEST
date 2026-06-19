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

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[experimental]"



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


    // TODO:
    // - fix this guard! Just runtime skip 
    // - fix test
    

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

    SECTION( LABEL_CORRECTNESS ) {

        // iterate the cached Quregs so the save path is exercised under every
        // deployment combination (serial, OMP, MPI, GPU and their mixtures);
        // each restored Qureg chooses its own deployment independently
        SECTION( LABEL_STATEVEC ) {

            for (auto& [label, q] : getCachedStatevecs()) {
                DYNAMIC_SECTION( label ) {

                    initRandomPureState(q);

                    saveQuregToFile(q, SV_FILE);
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

                    initRandomPureState(q); // works even for density matrices

                    saveQuregToFile(q, DM_FILE);
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
