/** @file
 * Unit tests of the qureg module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitqureg Qureg
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/random.hpp"

#include <string>
#include <algorithm>

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[qureg]"



/** 
 * TESTS
 * 
 * @ingroup unitqureg
 * @{
 */


TEST_CASE( "createQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1, 10) );
        CAPTURE( numQubits );

        Qureg qureg = createQureg(numQubits);
        QuESTEnv env = getQuESTEnv();

        // check fixed fields
        REQUIRE( qureg.numQubits == numQubits );
        REQUIRE( qureg.isDensityMatrix == 0 );
        REQUIRE( qureg.numAmps == getPow2(numQubits) );
        REQUIRE( qureg.logNumAmps == getLog2(qureg.numAmps) );
        REQUIRE( qureg.logNumAmpsPerNode == getLog2(qureg.numAmpsPerNode) );
        REQUIRE( (qureg.isMultithreaded  == 0 || qureg.isMultithreaded  == 1) );
        REQUIRE( (qureg.isGpuAccelerated == 0 || qureg.isGpuAccelerated == 1) );
        REQUIRE( (qureg.isDistributed    == 0 || qureg.isDistributed    == 1) );
        
        // check deployment-specific fields
        if (qureg.isDistributed) {
            REQUIRE( qureg.rank == env.rank );
            REQUIRE( qureg.numNodes == env.numNodes );
            REQUIRE( qureg.logNumNodes == getLog2(env.numNodes) );
            REQUIRE( qureg.numAmpsPerNode == qureg.numAmps / env.numNodes );
        } else {
            REQUIRE( qureg.rank == 0 );
            REQUIRE( qureg.numNodes == 1 );
            REQUIRE( qureg.logNumNodes == 0 );
            REQUIRE( qureg.numAmpsPerNode == qureg.numAmps );
        }

        // check memory allocs
        REQUIRE( qureg.cpuAmps != nullptr );
        if (qureg.isGpuAccelerated)
            REQUIRE( qureg.gpuAmps != nullptr );
        if (qureg.isDistributed)
            REQUIRE( qureg.cpuCommBuffer != nullptr );
        if (qureg.isGpuAccelerated && qureg.isDistributed)
            REQUIRE( qureg.gpuCommBuffer != nullptr );

        /// @todo check autodeployments

        // check begins in zero state
        qvector ref = getZeroVector(qureg.numAmps); ref[0] = 1; // |0>
        REQUIRE_AGREE( qureg, ref );

        destroyQureg(qureg);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createQureg(numQubits), ContainsSubstring("must contain one or more qubits") );
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createQureg(100), ContainsSubstring("maximum which can be addressed by the qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createQureg(62), ContainsSubstring("memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createQureg(50), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createDensityQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1, 7) ); // 7qb densematr = 14qb statevec
        CAPTURE( numQubits );

        Qureg qureg = createDensityQureg(numQubits);
        QuESTEnv env = getQuESTEnv();

        // check fixed fields
        REQUIRE( qureg.numQubits == numQubits );
        REQUIRE( qureg.isDensityMatrix == 1 );
        REQUIRE( qureg.numAmps == getPow2(2*numQubits) );
        REQUIRE( qureg.logNumAmps == getLog2(qureg.numAmps) );
        REQUIRE( qureg.logNumAmpsPerNode == getLog2(qureg.numAmpsPerNode) );
        REQUIRE( qureg.logNumNodes == getLog2(qureg.numNodes) );
        REQUIRE( qureg.logNumColsPerNode == getLog2(getPow2(numQubits) / qureg.numNodes) );
        REQUIRE( qureg.numAmpsPerNode == qureg.numAmps / qureg.numNodes );
        REQUIRE( (qureg.isMultithreaded  == 0 || qureg.isMultithreaded  == 1) );
        REQUIRE( (qureg.isGpuAccelerated == 0 || qureg.isGpuAccelerated == 1) );
        REQUIRE( (qureg.isDistributed    == 0 || qureg.isDistributed    == 1) );
        
        // check deployment-specific fields
        if (qureg.isDistributed) {
            REQUIRE( qureg.rank == env.rank );
            REQUIRE( qureg.numNodes == env.numNodes );
        } else {
            REQUIRE( qureg.rank == 0 );
            REQUIRE( qureg.numNodes == 1 );
        }

        // check memory allocs
        REQUIRE( qureg.cpuAmps != nullptr );
        if (qureg.isGpuAccelerated)
            REQUIRE( qureg.gpuAmps != nullptr );
        if (qureg.isDistributed)
            REQUIRE( qureg.cpuCommBuffer != nullptr );
        if (qureg.isGpuAccelerated && qureg.isDistributed)
            REQUIRE( qureg.gpuCommBuffer != nullptr );

        /// @todo check autodeployments

        // check begins in zero state
        qmatrix ref = getZeroMatrix(getPow2(numQubits)); ref[0][0] = 1; // |0><0|
        REQUIRE_AGREE( qureg, ref );

        destroyQureg(qureg);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createDensityQureg(numQubits), ContainsSubstring("must contain one or more qubits") );
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createDensityQureg(50), ContainsSubstring("qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createDensityQureg(31), ContainsSubstring("memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createDensityQureg(25), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createForcedQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        QuESTEnv env = getQuESTEnv();

        int minNumQubits = std::max({1, env.isDistributed? getLog2(env.numNodes) : 1});
        int maxNumQubits = std::min({minNumQubits + 10, 14});
        int numQubits = GENERATE_COPY( range(minNumQubits, maxNumQubits) );
        CAPTURE( numQubits );

        Qureg qureg = createForcedQureg(numQubits);

        // check fixed fields
        REQUIRE( qureg.numQubits == numQubits );
        REQUIRE( qureg.isDensityMatrix == 0 );
        REQUIRE( qureg.numAmps == getPow2(numQubits) );
        REQUIRE( qureg.logNumAmps == getLog2(qureg.numAmps) );
        REQUIRE( qureg.logNumAmpsPerNode == getLog2(qureg.numAmpsPerNode) );

        // check deployments
        REQUIRE( qureg.numNodes == env.numNodes );
        REQUIRE( qureg.logNumNodes == getLog2(env.numNodes) );
        REQUIRE( qureg.isMultithreaded  == env.isMultithreaded);
        REQUIRE( qureg.isGpuAccelerated == env.isGpuAccelerated);
        REQUIRE( (qureg.isDistributed   == env.isDistributed || env.numNodes == 1) ); // permit auto-disable MPI
        
        // check deployment-specific fields
        if (qureg.isDistributed) {
            REQUIRE( qureg.rank == env.rank );
            REQUIRE( qureg.numAmpsPerNode == qureg.numAmps / env.numNodes );
        } else {
            REQUIRE( qureg.rank == 0 );
            REQUIRE( qureg.numAmpsPerNode == qureg.numAmps );
        }

        // check memory allocs
        REQUIRE( qureg.cpuAmps != nullptr );
        if (qureg.isGpuAccelerated)
            REQUIRE( qureg.gpuAmps != nullptr );
        if (qureg.isDistributed)
            REQUIRE( qureg.cpuCommBuffer != nullptr );
        if (qureg.isGpuAccelerated && qureg.isDistributed)
            REQUIRE( qureg.gpuCommBuffer != nullptr );

        // check begins in zero state
        qvector ref = getZeroVector(getPow2(numQubits)); ref[0] = 1; // |0>
        REQUIRE_AGREE( qureg, ref );

        destroyQureg(qureg);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            REQUIRE_THROWS_WITH( createForcedQureg(-1), ContainsSubstring("must contain one or more qubits") );
            REQUIRE_THROWS_WITH( createForcedQureg(+0), ContainsSubstring("must contain one or more qubits") );

            int numNodes = getQuESTEnv().numNodes;
            if (numNodes > 2) { // nodes=2 => min=1
                int minNumQubits = getLog2(numNodes);
                REQUIRE_THROWS_WITH( createForcedQureg(minNumQubits-1), ContainsSubstring("each node would contain fewer than one amplitude") );
            }
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createForcedQureg(100), ContainsSubstring("maximum which can be addressed by the qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createForcedQureg(62), ContainsSubstring("memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createForcedQureg(50), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createForcedDensityQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        QuESTEnv env = getQuESTEnv();

        int minNumQubits = std::max({1, env.isDistributed? getLog2(env.numNodes) : 1});
        int maxNumQubits = std::min({minNumQubits + 5, 8}); // 8qb densmatr = 16qb statevec
        int numQubits = GENERATE_COPY( range(minNumQubits, maxNumQubits) );
        CAPTURE( numQubits );

        Qureg qureg = createForcedDensityQureg(numQubits);

        // check fixed fields
        REQUIRE( qureg.numQubits == numQubits );
        REQUIRE( qureg.isDensityMatrix == 1 );
        REQUIRE( qureg.numAmps == getPow2(2*numQubits) );
        REQUIRE( qureg.logNumAmps == getLog2(qureg.numAmps) );
        REQUIRE( qureg.logNumAmpsPerNode == getLog2(qureg.numAmpsPerNode) );
        REQUIRE( qureg.logNumColsPerNode == getLog2(getPow2(numQubits) / qureg.numNodes) );

        // check deployments
        REQUIRE( qureg.numNodes == env.numNodes );
        REQUIRE( qureg.logNumNodes == getLog2(env.numNodes) );
        REQUIRE( qureg.isMultithreaded  == env.isMultithreaded);
        REQUIRE( qureg.isGpuAccelerated == env.isGpuAccelerated);
        REQUIRE( (qureg.isDistributed   == env.isDistributed || env.numNodes == 1) ); // permit auto-disable MPI
        
        // check deployment-specific fields
        if (qureg.isDistributed) {
            REQUIRE( qureg.rank == env.rank );
            REQUIRE( qureg.numAmpsPerNode == qureg.numAmps / env.numNodes );
        } else {
            REQUIRE( qureg.rank == 0 );
            REQUIRE( qureg.numAmpsPerNode == qureg.numAmps );
        }

        // check memory allocs
        REQUIRE( qureg.cpuAmps != nullptr );
        if (qureg.isGpuAccelerated)
            REQUIRE( qureg.gpuAmps != nullptr );
        if (qureg.isDistributed)
            REQUIRE( qureg.cpuCommBuffer != nullptr );
        if (qureg.isGpuAccelerated && qureg.isDistributed)
            REQUIRE( qureg.gpuCommBuffer != nullptr );

        // check begins in zero state
        qmatrix ref = getZeroMatrix(getPow2(numQubits)); ref[0][0] = 1; // |0><0|
        REQUIRE_AGREE( qureg, ref );

        destroyQureg(qureg);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            REQUIRE_THROWS_WITH( createForcedDensityQureg(-1), ContainsSubstring("must contain one or more qubits") );
            REQUIRE_THROWS_WITH( createForcedDensityQureg(+0), ContainsSubstring("must contain one or more qubits") );

            int numNodes = getQuESTEnv().numNodes;
            if (numNodes > 2) { // nodes=2 => min=1
                int minNumQubits = getLog2(numNodes);
                REQUIRE_THROWS_WITH( createForcedDensityQureg(minNumQubits-1), ContainsSubstring("each node would contain fewer than a column") );
            }
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createForcedDensityQureg(50), ContainsSubstring("qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createForcedDensityQureg(31), ContainsSubstring("memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createForcedDensityQureg(25), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createCustomQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        QuESTEnv env = getQuESTEnv();

        for (auto [deployLabel, mpi, gpu, omp] : getSupportedDeployments()) {

            int isDenseMatr = GENERATE( 0, 1 );

            std::string quregLabel = (isDenseMatr)? LABEL_STATEVEC : LABEL_DENSMATR; 
            std::string secLabel = quregLabel + LABEL_DELIMITER + deployLabel;

            SECTION( secLabel ) {

                int minNumQubits = std::max({1, mpi? getLog2(env.numNodes) : 1});
                int maxNumQubits = std::min({minNumQubits + 3, isDenseMatr? 6 : 12});
                int numQubits = GENERATE_COPY( range(minNumQubits, maxNumQubits) );
                CAPTURE( numQubits );

                Qureg qureg = createCustomQureg(numQubits, isDenseMatr, mpi, gpu, omp);

                // check fixed fields
                REQUIRE( qureg.numQubits == numQubits );
                REQUIRE( qureg.isDensityMatrix == isDenseMatr );
                REQUIRE( qureg.numAmps == getPow2(numQubits * (isDenseMatr? 2:1)) );
                REQUIRE( qureg.logNumNodes == getLog2(qureg.numNodes) );
                REQUIRE( qureg.logNumAmps == getLog2(qureg.numAmps) );
                REQUIRE( qureg.logNumAmpsPerNode == getLog2(qureg.numAmpsPerNode) );
                REQUIRE( qureg.logNumColsPerNode == (isDenseMatr? getLog2(getPow2(numQubits) / qureg.numNodes) : 0) );

                // check deployments
                REQUIRE( qureg.isMultithreaded  == omp );
                REQUIRE( qureg.isGpuAccelerated == gpu );
                REQUIRE( (qureg.isDistributed   == mpi || env.numNodes == 1) ); // permit auto-disable MPI
                
                // check deployment-specific fields
                if (mpi) {
                    REQUIRE( qureg.rank == env.rank );
                    REQUIRE( qureg.numNodes == env.numNodes );
                    REQUIRE( qureg.numAmpsPerNode == qureg.numAmps / env.numNodes );
                } else {
                    REQUIRE( qureg.rank == 0 );
                    REQUIRE( qureg.numNodes == 1 );
                    REQUIRE( qureg.numAmpsPerNode == qureg.numAmps );
                }

                // check memory allocs
                REQUIRE( qureg.cpuAmps != nullptr );
                if (qureg.isGpuAccelerated)
                    REQUIRE( qureg.gpuAmps != nullptr );
                if (qureg.isDistributed)
                    REQUIRE( qureg.cpuCommBuffer != nullptr );
                if (qureg.isGpuAccelerated && qureg.isDistributed)
                    REQUIRE( qureg.gpuCommBuffer != nullptr );

                // check begins in zero state
                if (isDenseMatr) {
                    qmatrix ref = getZeroMatrix(getPow2(numQubits)); ref[0][0] = 1; // |0><0|
                    REQUIRE_AGREE( qureg, ref );
                } else {
                    qvector ref = getZeroVector(getPow2(numQubits)); ref[0] = 1; // |0>
                    REQUIRE_AGREE( qureg, ref );
                }

                destroyQureg(qureg);
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "invalid isDensityMatrix flag" ) {

            int isDensMatr = GENERATE( -1, 2 );

            REQUIRE_THROWS_WITH( createCustomQureg(10, isDensMatr, 0,0,0), ContainsSubstring("isDensityMatrix") );
        }

        SECTION( "unsupported deployment" ) {

            QuESTEnv env = getQuESTEnv();
            const int numQubits = 10;
            const int isDensMatr = 0;

            if (!env.isDistributed)
                REQUIRE_THROWS_WITH( createCustomQureg(numQubits,isDensMatr, 1,0,0), ContainsSubstring("non-distributed") );
            
            if (!env.isGpuAccelerated)
                REQUIRE_THROWS_WITH( createCustomQureg(numQubits,isDensMatr, 0,1,0), ContainsSubstring("non-GPU") );
            
            if (!env.isMultithreaded)
                REQUIRE_THROWS_WITH( createCustomQureg(numQubits,isDensMatr, 0,0,1), ContainsSubstring("non-multithreaded") );
        }

        SECTION( "too few qubits" ) {

            REQUIRE_THROWS_WITH( createCustomQureg(-1, 0,0,0,0), ContainsSubstring("must contain one or more qubits") );
            REQUIRE_THROWS_WITH( createCustomQureg(+0, 0,0,0,0), ContainsSubstring("must contain one or more qubits") );

            int numNodes = getQuESTEnv().numNodes;
            if (numNodes > 2) { // nodes=2 => min=1
                int minNumQubits = getLog2(numNodes);
                int useDistrib = 1;
                REQUIRE_THROWS_WITH( createCustomQureg(minNumQubits-1,0, useDistrib,0,0), ContainsSubstring("each node would contain fewer than one amplitude") );
            }
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createCustomQureg(100, 0, 0,0,0), ContainsSubstring("qindex") );
            REQUIRE_THROWS_WITH( createCustomQureg(50,  1, 0,0,0), ContainsSubstring("qindex") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createCustomQureg(62, 0, 0,0,0), ContainsSubstring("memory would overflow size_t") );
            REQUIRE_THROWS_WITH( createCustomQureg(31, 1, 0,0,0), ContainsSubstring("memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createCustomQureg(50, 0, 0,0,0), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            REQUIRE_THROWS_WITH( createCustomQureg(25, 1, 0,0,0), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createCloneQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        auto cache = GENERATE( getCachedStatevecs(), getCachedDensmatrs() );

        for (auto& [label, qureg]: cache) {

            initRandomPureState(qureg);
            Qureg clone = createCloneQureg(qureg);

            // check identical fields
            REQUIRE( clone.isMultithreaded   == qureg.isMultithreaded );
            REQUIRE( clone.isGpuAccelerated  == qureg.isGpuAccelerated );
            REQUIRE( clone.isDistributed     == qureg.isDistributed );
            REQUIRE( clone.rank              == qureg.rank );
            REQUIRE( clone.numNodes          == qureg.numNodes );
            REQUIRE( clone.logNumNodes       == qureg.logNumNodes );
            REQUIRE( clone.isDensityMatrix   == qureg.isDensityMatrix );
            REQUIRE( clone.numQubits         == qureg.numQubits );
            REQUIRE( clone.numAmps           == qureg.numAmps );
            REQUIRE( clone.logNumAmps        == qureg.logNumAmps );
            REQUIRE( clone.numAmpsPerNode    == qureg.numAmpsPerNode );
            REQUIRE( clone.logNumAmpsPerNode == qureg.logNumAmpsPerNode );
            REQUIRE( clone.logNumColsPerNode == qureg.logNumColsPerNode );

            // check memory is different
            REQUIRE( clone.cpuAmps != qureg.cpuAmps );
            if (clone.isGpuAccelerated)
                REQUIRE( clone.gpuAmps != qureg.gpuAmps );
            if (clone.isDistributed)
                REQUIRE( clone.cpuCommBuffer != qureg.cpuCommBuffer );
            if (clone.isGpuAccelerated && clone.isDistributed)
                REQUIRE( clone.gpuCommBuffer != qureg.gpuCommBuffer );

            // check identical states
            REQUIRE_AGREE( clone, qureg );

            destroyQureg(clone);
        }
    }
}


TEST_CASE( "destroyQureg", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        Qureg qureg = createQureg(5);
        REQUIRE_NOTHROW( destroyQureg(qureg) );
    }

    SECTION( LABEL_VALIDATION ) {

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            Qureg qureg;
            REQUIRE_THROWS_WITH( destroyQureg(qureg), ContainsSubstring("invalid Qureg") );
        }
        #endif
        #endif
    }
}


TEST_CASE( "getQuregAmp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (auto& [label, qureg]: getCachedStatevecs()) {

            qvector ref = getRefStatevec();
            initDebugState(qureg);
            setToDebugState(ref);

            int index = GENERATE_COPY( range(0, (int) qureg.numAmps) );
            CAPTURE( index );

            // not necessarily identical when qureg is GPU-accelerated
            REQUIRE_AGREE( getQuregAmp(qureg, index), ref[index] );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            Qureg qureg;
            REQUIRE_THROWS_WITH( getQuregAmp(qureg,0), ContainsSubstring("invalid Qureg") );
        }
        #endif
        #endif

        SECTION( "invalid index" ) {

            int numQubits = GENERATE( range(1,5) );
            Qureg qureg = createQureg(numQubits);
            int index = GENERATE_COPY( -1, getPow2(numQubits), getPow2(numQubits) + 1 );

            REQUIRE_THROWS_WITH( getQuregAmp(qureg,index), ContainsSubstring("state index") && ContainsSubstring("is invalid") );

            destroyQureg(qureg);
        }

        SECTION( "densitymatrix" ) {

            Qureg qureg = createDensityQureg(1);

            REQUIRE_THROWS_WITH( getQuregAmp(qureg,0), ContainsSubstring("received a density matrix") );

            destroyQureg(qureg);
        }
    }
}


TEST_CASE( "getDensityQuregAmp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (auto& [label, qureg]: getCachedDensmatrs()) {

            qmatrix ref = getRefDensmatr();
            initDebugState(qureg);
            setToDebugState(ref);

            int dim = (int) getPow2(qureg.numQubits);
            int row = getRandomInt(0, dim);
            int col = getRandomInt(0, dim);

            GENERATE( range(0,100) );
            CAPTURE( row, col );

            // not necessarily identical when qureg is GPU-accelerated
            REQUIRE_AGREE( getDensityQuregAmp(qureg, row, col), ref[row][col] );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            Qureg qureg;
            REQUIRE_THROWS_WITH( getDensityQuregAmp(qureg,0,0), ContainsSubstring("invalid Qureg") );
        }
        #endif
        #endif

        SECTION( "invalid indices" ) {

            int numQubits = GENERATE( range(1,5) );
            Qureg qureg = createDensityQureg(numQubits);
            qindex index = GENERATE_COPY( -1, getPow2(numQubits), getPow2(numQubits) + 1 );

            REQUIRE_THROWS_WITH( getDensityQuregAmp(qureg,0,index), ContainsSubstring("A") );
            REQUIRE_THROWS_WITH( getDensityQuregAmp(qureg,index,0), ContainsSubstring("B") );

            destroyQureg(qureg);
        }

        SECTION( "statevector" ) {

            Qureg qureg = createQureg(1);

            REQUIRE_THROWS_WITH( getDensityQuregAmp(qureg,0,0), ContainsSubstring("received a statevector") );

            destroyQureg(qureg);
        }
    }
}


TEST_CASE( "getQuregAmps", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (auto& [label, qureg]: getCachedStatevecs()) {

            qvector ref = getRefStatevec();
            initDebugState(qureg);
            setToDebugState(ref);

            GENERATE( range(0,100) );
            int startInd = getRandomInt(0, qureg.numAmps);
            int numAmps = getRandomInt(0, qureg.numAmps - startInd + 1);

            vector<qcomp> out(numAmps);
            getQuregAmps(out.data(), qureg, startInd, numAmps);

            // not necessarily identical when qureg is GPU-accelerated
            REQUIRE_AGREE( out, getSublist(ref, startInd, numAmps) );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = createQureg(5);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            Qureg bad;
            REQUIRE_THROWS_WITH( getQuregAmps(nullptr,bad,0,0), ContainsSubstring("invalid Qureg") );
        }
        #endif
        #endif

        SECTION( "indices" ) {

            int startInd = GENERATE_COPY( -1, qureg.numAmps, qureg.numAmps + 1 );

            REQUIRE_THROWS_WITH( getQuregAmps(nullptr,qureg,startInd,0), ContainsSubstring("starting basis state index") );
        }

        SECTION( "num amps") {

            int numOut = GENERATE_COPY( -1, qureg.numAmps + 1 );

            REQUIRE_THROWS_WITH( getQuregAmps(nullptr,qureg,0,numOut), ContainsSubstring("number of amplitudes") );
        }

        SECTION( "subrange" ) {

            int numOut = 10;

            REQUIRE_THROWS_WITH( getQuregAmps(nullptr,qureg, qureg.numAmps - numOut, numOut + 1), ContainsSubstring("implies an end index") );
        }
        
        destroyQureg(qureg);
    }
}


TEST_CASE( "getDensityQuregAmps", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (auto& [label, qureg]: getCachedDensmatrs()) {

            qmatrix ref = getRefDensmatr();
            initDebugState(qureg);
            setToDebugState(ref);

            GENERATE( range(0,100) );

            auto startRow = getRandomInt(0, ref.size());
            auto startCol = getRandomInt(0, ref.size());
            auto numRows = getRandomInt(1, ref.size() - startRow);
            auto numCols = getRandomInt(1, ref.size() - startCol);

            qcomp** out = (qcomp**) malloc(numRows * sizeof *out);
            for (int i=0; i<numRows; i++)
                out[i] = (qcomp*) malloc(numCols * sizeof **out);

            getDensityQuregAmps(out, qureg, startRow, startCol, numRows, numCols);

            bool agrees = true;
            for (int r=0; r<numRows && agrees; r++)
                for (int c=0; c<numCols && agrees; c++)
                    agrees = doScalarsAgree(out[r][c], ref[startRow+r][startCol+c]);

            REQUIRE( agrees );

            for (int i=0; i<numRows; i++)
                free(out[i]);
            free(out);
        }
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = createDensityQureg(5);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            Qureg bad;
            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,bad,0,0,0,0), ContainsSubstring("invalid Qureg") );
        }
        #endif
        #endif

        SECTION( "indices" ) {

            int startInd = GENERATE_COPY( -1, qureg.numAmps, qureg.numAmps + 1 );

            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,qureg, startInd,0, 0,0), ContainsSubstring("Either or both of the starting row and column") );
            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,qureg, 0,startInd, 0,0), ContainsSubstring("Either or both of the starting row and column") );
        }

        SECTION( "num amps") {

            int numOut = GENERATE_COPY( -1, qureg.numAmps + 1 );

            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,qureg, 0,0, 0,numOut), ContainsSubstring("Either or both of the number of rows and columns") );
            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,qureg, 0,0, numOut,0), ContainsSubstring("Either or both of the number of rows and columns") );
        }

        SECTION( "subrange" ) {

            int numOut = 10;
            int dim = getPow2(qureg.numQubits);

            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,qureg, dim-numOut,0, numOut + 1,1), ContainsSubstring("combination of starting") && ContainsSubstring("number of rows and columns") );
            REQUIRE_THROWS_WITH( getDensityQuregAmps(nullptr,qureg, 0,dim-numOut, 1,numOut + 1), ContainsSubstring("combination of starting") && ContainsSubstring("number of rows and columns") );
        }
        
        destroyQureg(qureg);
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void reportQuregParams(Qureg qureg);
void reportQureg(Qureg qureg);

void syncQuregToGpu  (Qureg qureg);
void syncQuregFromGpu(Qureg qureg);

void syncSubQuregToGpu  (Qureg qureg, qindex localStartInd, qindex numLocalAmps);
void syncSubQuregFromGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps);
