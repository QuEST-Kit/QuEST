/** @file
 * Unit tests of Qureg checkpointing (saveQuregToFile / createQuregFromFile).
 *
 * These tests are only compiled when QuEST is built with the CMake option
 * -DQUEST_ENABLE_CHECKPOINTING=ON (which additionally requires the ADIOS2 library).
 *
 * @author Ashmit JaiSarita Gupta
 *
 * @defgroup unitcheckpoint Checkpointing
 * @ingroup unittests
 */

#include "quest.h"

#if QUEST_COMPILE_CHECKPOINTING

#include <catch2/catch_test_macros.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <string>



/*
 * file constants and helpers
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[checkpoint]"

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



/** TESTS
 *
 * @ingroup unitcheckpoint
 * @{
 */

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
        // file only compiles under QUEST_COMPILE_CHECKPOINTING. ADIOS2's own
        // runtime errors (e.g. a missing file) are not QuEST validation errors.
        SUCCEED( );
    }
}

/** @} (end defgroup) */

#endif // QUEST_COMPILE_CHECKPOINTING
