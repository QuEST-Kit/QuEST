/** @file
 * Unit tests of Qureg checkpointing (saveQuregToFile / createQuregFromFile).
 *
 * These tests are only compiled when QuEST is built with the CMake option
 * -DENABLE_CHECKPOINTING=ON (which additionally requires the ADIOS2 library).
 *
 * @author Ashmit JaiSarita Gupta
 *
 * @defgroup unitcheckpoint Checkpointing
 * @ingroup unittests
 */

#include "quest.h"

#if QUEST_COMPILE_CHECKPOINTING

#include <catch2/catch_test_macros.hpp>

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
}

TEST_CASE( "saveQuregToFile and createQuregFromFile", "[checkpoint]" ) {

    SECTION( "statevector round-trip preserves dimension and amplitudes" ) {

        Qureg q = createQureg(6);
        initRandomPureState(q);

        saveQuregToFile(q, SV_FILE);
        Qureg r = createQuregFromFile(SV_FILE);

        CHECK( r.numQubits       == q.numQubits );
        CHECK( r.isDensityMatrix == q.isDensityMatrix );
        CHECK( maxStatevectorAmpDiff(q, r) < 1e-12 );

        destroyQureg(q);
        destroyQureg(r);

        // In distributed runs every node opened the same shared file, so only one
        // may delete it; a barrier first guarantees all nodes have finished
        // reading, and a barrier after keeps the next section's collective write
        // from racing a half-removed directory.
        syncQuESTEnv();
        if (getQuESTEnv().rank == 0)
            std::filesystem::remove_all(SV_FILE);
        syncQuESTEnv();
    }

    SECTION( "density-matrix round-trip preserves dimension and amplitudes" ) {

        Qureg q = createDensityQureg(4);
        initZeroState(q);
        for (int t = 0; t < q.numQubits; t++)
            applyHadamard(q, t);
        applyT(q, 0);
        applyControlledPauliX(q, 0, 1);

        saveQuregToFile(q, DM_FILE);
        Qureg r = createQuregFromFile(DM_FILE);

        CHECK( r.numQubits       == q.numQubits );
        CHECK( r.isDensityMatrix == q.isDensityMatrix );
        CHECK( maxDensityMatrixAmpDiff(q, r) < 1e-12 );

        destroyQureg(q);
        destroyQureg(r);

        // see the statevector section: one node deletes, barriers bracket cleanup
        syncQuESTEnv();
        if (getQuESTEnv().rank == 0)
            std::filesystem::remove_all(DM_FILE);
        syncQuESTEnv();
    }
}

#endif // QUEST_COMPILE_CHECKPOINTING
