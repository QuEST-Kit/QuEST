/** @file
 * Unit tests of the matrices module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitmatr Matrices
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
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

#include <cstdlib>
#include <algorithm>

using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[matrices]"



/** 
 * TESTS
 * 
 * @ingroup unitmatr
 * @{
 */


TEST_CASE( "getCompMatr1", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        constexpr int dim = 2;
        qmatrix ref = getRandomMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // 2D compile-time array
            qcomp arr[dim][dim] = {{ref[0][0], ref[0][1]}, {ref[1][0], ref[1][1]}};
            REQUIRE_AGREE( getCompMatr1(arr), ref );

            // nested pointers
            qcomp** ptrs = (qcomp**) malloc(dim * sizeof *ptrs);
            for (int i=0; i<dim; i++) {
                ptrs[i] = (qcomp*) malloc(dim * sizeof **ptrs);
                for (int j=0; j<dim; j++)
                    ptrs[i][j] = ref[i][j];
            }
            REQUIRE_AGREE( getCompMatr1(ptrs), ref );

            // array of pointers
            qcomp* ptrArr[dim];
            for (int i=0; i<dim; i++)
                ptrArr[i] = ptrs[i];
            REQUIRE_AGREE( getCompMatr1(ptrArr), ref );

            // cleanup
            for (int i=0; i<dim; i++)
                free(ptrs[i]);
            free(ptrs);
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // nested vectors
            REQUIRE_AGREE( getCompMatr1(ref), ref );

            // inline vectors
            REQUIRE_AGREE( getCompMatr1({{ref[0][0],ref[0][1]},{ref[1][0],ref[1][1]}}), ref );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "null pointer" ) {

            qcomp** ptr = nullptr;
            REQUIRE_THROWS_WITH( getCompMatr1(ptr), ContainsSubstring("was a null pointer") );

            qcomp* arr[1] = {nullptr};
            REQUIRE_THROWS_WITH( getCompMatr1(arr), ContainsSubstring("contained a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_THROWS_WITH( getCompMatr1({{0,0}}),             ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( getCompMatr1({{0,0},{0,0},{0,0}}), ContainsSubstring("Incompatible number of rows") );
        
            REQUIRE_THROWS_WITH( getCompMatr1({{0,0},{0}}),     ContainsSubstring("incompatible number of elements") );
            REQUIRE_THROWS_WITH( getCompMatr1({{0,0},{0,0,0}}), ContainsSubstring("incompatible number of elements") );
        }
    }
}


TEST_CASE( "getCompMatr2", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        constexpr int dim = 4;
        qmatrix ref = getRandomMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // 2D compile-time array
            qcomp arr[dim][dim];
            for (int i=0; i<dim; i++)
                for (int j=0; j<dim; j++)
                    arr[i][j] = ref[i][j];
            REQUIRE_AGREE( getCompMatr2(arr), ref );

            // nested pointers
            qcomp** ptrs = (qcomp**) malloc(dim * sizeof *ptrs);
            for (int i=0; i<dim; i++) {
                ptrs[i] = (qcomp*) malloc(dim * sizeof **ptrs);
                for (int j=0; j<dim; j++)
                    ptrs[i][j] = ref[i][j];
            }
            REQUIRE_AGREE( getCompMatr2(ptrs), ref );

            // array of pointers
            qcomp* ptrArr[dim];
            for (int i=0; i<dim; i++)
                ptrArr[i] = ptrs[i];
            REQUIRE_AGREE( getCompMatr2(ptrArr), ref );

            // cleanup
            for (int i=0; i<dim; i++)
                free(ptrs[i]);
            free(ptrs);
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // nested vectors
            REQUIRE_AGREE( getCompMatr2(ref), ref );

            // inline vectors
            REQUIRE_AGREE( 
                getCompMatr2({
                    {ref[0][0],ref[0][1],ref[0][2],ref[0][3]},
                    {ref[1][0],ref[1][1],ref[1][2],ref[1][3]},
                    {ref[2][0],ref[2][1],ref[2][2],ref[2][3]},
                    {ref[3][0],ref[3][1],ref[3][2],ref[3][3]}
                }),
                ref
            );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "null pointer" ) {

            qcomp** ptr = nullptr;
            REQUIRE_THROWS_WITH( getCompMatr2(ptr), ContainsSubstring("was a null pointer") );

            qcomp* arr[1] = {nullptr};
            REQUIRE_THROWS_WITH( getCompMatr2(arr), ContainsSubstring("contained a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            qmatrix bad1 = {{0,0,0,0}};
            qmatrix bad2 = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};

            REQUIRE_THROWS_WITH( getCompMatr2(bad1), ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( getCompMatr2(bad2), ContainsSubstring("Incompatible number of rows") );
        
            REQUIRE_THROWS_WITH( getCompMatr2({{0,0,0,0},{0},{0},{0}}),   ContainsSubstring("incompatible number of elements") );
            REQUIRE_THROWS_WITH( getCompMatr2({{0,0,0,0,0},{0},{0},{0}}), ContainsSubstring("incompatible number of elements") );
        }
    }
}


TEST_CASE( "getDiagMatr1", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        constexpr int dim = 2;
        qmatrix ref = getRandomDiagonalMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // compile-time array
            qcomp arr[dim] = {ref[0][0], ref[1][1] };
            REQUIRE_AGREE( getDiagMatr1(arr), ref );

            // pointer
            qvector diags = getDiagonals(ref);
            REQUIRE_AGREE( getDiagMatr1(diags.data()), ref );
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // inline vectors
            REQUIRE_AGREE( getDiagMatr1({ref[0][0],ref[1][1]}), ref );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "null pointer" ) {

            qcomp* ptr = nullptr;
            REQUIRE_THROWS_WITH( getDiagMatr1(ptr), ContainsSubstring("was a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface
            qvector bad1 = {0};
            qvector bad2 = {0,0,0};

            REQUIRE_THROWS_WITH( getDiagMatr1(bad1), ContainsSubstring("Incompatible number of elements") );
            REQUIRE_THROWS_WITH( getDiagMatr1(bad2), ContainsSubstring("Incompatible number of elements") );
        }
    }
}


TEST_CASE( "getDiagMatr2", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        constexpr int dim = 4;
        qmatrix ref = getRandomDiagonalMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // compile-time array
            qcomp arr[dim] = {ref[0][0], ref[1][1], ref[2][2], ref[3][3] };
            REQUIRE_AGREE( getDiagMatr2(arr), ref );

            // pointer
            qvector diags = getDiagonals(ref);
            REQUIRE_AGREE( getDiagMatr2(diags.data()), ref );
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // inline vectors
            REQUIRE_AGREE( getDiagMatr2({ref[0][0], ref[1][1], ref[2][2], ref[3][3] }), ref );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "null pointer" ) {

            qcomp* ptr = nullptr;
            REQUIRE_THROWS_WITH( getDiagMatr2(ptr), ContainsSubstring("was a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface
            qvector bad1 = {0,0,0};
            qvector bad2 = {0,0,0,0,0};

            REQUIRE_THROWS_WITH( getDiagMatr2(bad1), ContainsSubstring("Incompatible number of elements") );
            REQUIRE_THROWS_WITH( getDiagMatr2(bad2), ContainsSubstring("Incompatible number of elements") );
        }
    }
}


TEST_CASE( "getInlineCompMatr1", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix ref = {{1,2},{3_i,4_i}};

        REQUIRE_AGREE( getInlineCompMatr1({{1,2},{3_i,4_i}}), ref );
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_THROWS_WITH( getInlineCompMatr1({{0,0}}),             ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( getInlineCompMatr1({{0,0},{0,0},{0,0}}), ContainsSubstring("Incompatible number of rows") );
        
            REQUIRE_THROWS_WITH( getInlineCompMatr1({{0,0},{0}}),     ContainsSubstring("incompatible number of elements") );
            REQUIRE_THROWS_WITH( getInlineCompMatr1({{0,0},{0,0,0}}), ContainsSubstring("incompatible number of elements") );
        }
    }
}


TEST_CASE( "getInlineCompMatr2", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix ref = {
            {1,2,3,4},
            {5_i,6_i,7_i,8_i},
            {-1,-2,-3,-4},
            {-1_i,-2_i,-3_i,-4_i}};

        REQUIRE_AGREE( 
            getInlineCompMatr2({
                {1,2,3,4},
                {5_i,6_i,7_i,8_i},
                {-1,-2,-3,-4},
                {-1_i,-2_i,-3_i,-4_i}
            }), ref );
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_THROWS_WITH( getInlineCompMatr2({{0,0,0,0}}),                     ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( getInlineCompMatr2({{0,0,0,0},{0,0,0,0},{0,0,0,0}}), ContainsSubstring("Incompatible number of rows") );
        
            REQUIRE_THROWS_WITH( getInlineCompMatr2({{0,0,0,0},{0},{0},{0}}),   ContainsSubstring("incompatible number of elements") );
            REQUIRE_THROWS_WITH( getInlineCompMatr2({{0,0,0,0,0},{0},{0},{0}}), ContainsSubstring("incompatible number of elements") );
        }
    }
}


TEST_CASE( "getInlineDiagMatr1", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix ref = getDiagonalMatrix({1,2_i});

        REQUIRE_AGREE( getInlineDiagMatr1({1,2_i}), ref );
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "invalid dimensions" ) {

            REQUIRE_THROWS_WITH( getInlineDiagMatr1({0,0,0}), ContainsSubstring("Incompatible number of elements") );
        }
    }
}


TEST_CASE( "getInlineDiagMatr2", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix ref = getDiagonalMatrix({1,2_i,-3,-4_i});

        REQUIRE_AGREE( getInlineDiagMatr2({1,2_i,-3,-4_i}), ref );
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "invalid dimensions" ) {

            REQUIRE_THROWS_WITH( getInlineDiagMatr2({0,0,0}),     ContainsSubstring("Incompatible number of elements") );
            REQUIRE_THROWS_WITH( getInlineDiagMatr2({0,0,0,0,0}), ContainsSubstring("Incompatible number of elements") );
        }
    }
}


TEST_CASE( "createCompMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // inessential to link to cached Qureg sizes, but
        // provides a nice ballpark
        int maxNumQubits = getNumCachedQubits();
        int numQubits = GENERATE_COPY( range(1, maxNumQubits+1) );
        CAPTURE( numQubits );

        CompMatr matr = createCompMatr(numQubits);
        qmatrix blank = getZeroMatrix(getPow2(numQubits));

        // default elements
        REQUIRE_AGREE( matr, blank );

        // fields
        REQUIRE( matr.numQubits == numQubits );
        REQUIRE( matr.numRows == getPow2(numQubits) );

        // default properties
        REQUIRE( *(matr.isApproxUnitary)   == -1 ); // unknown
        REQUIRE( *(matr.isApproxHermitian) == -1 ); // unknown
        REQUIRE( *(matr.wasGpuSynced)      ==  0 ); // false

        // pointers
        REQUIRE( matr.cpuElems != nullptr );
        REQUIRE( matr.cpuElemsFlat != nullptr );

        if (getQuESTEnv().isGpuAccelerated)
            REQUIRE( matr.gpuElemsFlat != nullptr );
        else
            REQUIRE( matr.gpuElemsFlat == nullptr );
        
        destroyCompMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createCompMatr(numQubits), ContainsSubstring("must target one or more qubits") );
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createCompMatr(50), ContainsSubstring("maximum which can be addressed by the qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createCompMatr(31), ContainsSubstring("necessary memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createCompMatr(25), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // inessential to link to cached Qureg sizes, but
        // provides a nice ballpark
        int maxNumQubits = getNumCachedQubits();
        int numQubits = GENERATE_COPY( range(1, maxNumQubits+1) );
        CAPTURE( numQubits );

        DiagMatr matr = createDiagMatr(numQubits);
        qmatrix blank = getZeroMatrix(getPow2(numQubits));

        // default elements
        REQUIRE_AGREE( matr, blank );

        // fields
        REQUIRE( matr.numQubits == numQubits );
        REQUIRE( matr.numElems == getPow2(numQubits) );

        // default properties
        REQUIRE( *(matr.isApproxUnitary)       == -1 ); // unknown
        REQUIRE( *(matr.isApproxHermitian)     == -1 ); // unknown
        REQUIRE( *(matr.isApproxNonZero)       == -1 ); // unknown
        REQUIRE( *(matr.isStrictlyNonNegative) == -1 ); // unknown
        REQUIRE( *(matr.wasGpuSynced)          ==  0 ); // false

        // pointers
        REQUIRE( matr.cpuElems != nullptr );

        if (getQuESTEnv().isGpuAccelerated)
            REQUIRE( matr.gpuElems != nullptr );
        else
            REQUIRE( matr.gpuElems == nullptr );
        
        destroyDiagMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createDiagMatr(numQubits), ContainsSubstring("must target one or more qubits") );
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createDiagMatr(100), ContainsSubstring("maximum which can be addressed by the qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createDiagMatr(62), ContainsSubstring("necessary memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createDiagMatr(50), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "createFullStateDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE_COPY( range(1, 20) );
        CAPTURE( numQubits );

        // uses auto deployment
        FullStateDiagMatr matr = createFullStateDiagMatr(numQubits);

        // check state is blank
        bool isBlank = true;
        for (qindex i=0; i<matr.numElemsPerNode && isBlank; i++)
            isBlank = (matr.cpuElems[i] == qcomp(0,0));
        REQUIRE( isBlank );

        // dimensions
        REQUIRE( matr.numQubits == numQubits );
        REQUIRE( matr.numElems == getPow2(numQubits) );
        REQUIRE( matr.numElemsPerNode == matr.numElems / (matr.isDistributed? getQuESTEnv().numNodes : 1) );

        // default properties
        REQUIRE( *(matr.isApproxUnitary)       == -1 ); // unknown
        REQUIRE( *(matr.isApproxHermitian)     == -1 ); // unknown
        REQUIRE( *(matr.isApproxNonZero)       == -1 ); // unknown
        REQUIRE( *(matr.isStrictlyNonNegative) == -1 ); // unknown
        REQUIRE( *(matr.wasGpuSynced)          ==  0 ); // false

        // pointers
        REQUIRE( matr.cpuElems != nullptr );

        int gpu = matr.isGpuAccelerated;
        if (gpu) REQUIRE( matr.gpuElems != nullptr );
        else     REQUIRE( matr.gpuElems == nullptr );
        
        destroyFullStateDiagMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createFullStateDiagMatr(numQubits), ContainsSubstring("must target one or more qubits") );
        }

        SECTION( "too many qubits" ) {

            // overflows qindex in all precisions
            REQUIRE_THROWS_WITH( createFullStateDiagMatr(100), ContainsSubstring("maximum which can be addressed by the qindex type") );

            // overflows size_t in single precision (and ergo also in double and quad)
            REQUIRE_THROWS_WITH( createFullStateDiagMatr(62), ContainsSubstring("memory would overflow size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createFullStateDiagMatr(50), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }

        // this function chooses automatic deployment,
        // so we cannot force illegal distribution
    }
}


TEST_CASE( "createCustomFullStateDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (auto [label, mpi, gpu, omp] : getSupportedDeployments()) {

            DYNAMIC_SECTION( label ) {

                int minNumQubits = std::max({1, getLog2(getQuESTEnv().numNodes)});
                int maxNumQubits = std::min({20, minNumQubits + 1});
                int numQubits = GENERATE_COPY( range(minNumQubits, maxNumQubits+1) );
                CAPTURE( numQubits );

                FullStateDiagMatr matr = createCustomFullStateDiagMatr(numQubits, mpi, gpu, omp);

                // check state is blank
                bool isBlank = true;
                for (qindex i=0; i<matr.numElemsPerNode && isBlank; i++)
                    isBlank = (matr.cpuElems[i] == qcomp(0,0));
                REQUIRE( isBlank );

                // dimensions
                REQUIRE( matr.numQubits == numQubits );
                REQUIRE( matr.numElems == getPow2(numQubits) );
                REQUIRE( matr.numElemsPerNode == matr.numElems / (mpi? getQuESTEnv().numNodes : 1) );

                // accelerations
                REQUIRE( matr.isDistributed    == (int) (mpi && getQuESTEnv().numNodes > 1) );
                REQUIRE( matr.isGpuAccelerated == (int) gpu );
                REQUIRE( matr.isMultithreaded  == (int) omp );

                // default properties
                REQUIRE( *(matr.isApproxUnitary)       == -1 ); // unknown
                REQUIRE( *(matr.isApproxHermitian)     == -1 ); // unknown
                REQUIRE( *(matr.isApproxNonZero)       == -1 ); // unknown
                REQUIRE( *(matr.isStrictlyNonNegative) == -1 ); // unknown
                REQUIRE( *(matr.wasGpuSynced)          ==  0 ); // false

                // pointers
                REQUIRE( matr.cpuElems != nullptr );

                if (gpu) REQUIRE( matr.gpuElems != nullptr );
                else     REQUIRE( matr.gpuElems == nullptr );
                
                destroyFullStateDiagMatr(matr); 
            }
        }

    }

    SECTION( LABEL_VALIDATION ) {
    
        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            for (auto [label, mpi, gpu, omp] : getSupportedDeployments())
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(numQubits, mpi,gpu,omp), ContainsSubstring("must target one or more qubits") );

            if (getQuESTEnv().numNodes > 2) { // 1 node => min=1
                int minNumQubits = getLog2(getQuESTEnv().numNodes);
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(minNumQubits-1, 1,0,0), ContainsSubstring("node would contain fewer than one element") );
            }
        }

        SECTION( "too many qubits" ) {

            for (auto [label, mpi, gpu, omp] : getSupportedDeployments()) {

                // overflows qindex in all precisions
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(100, mpi,gpu,omp), ContainsSubstring("maximum which can be addressed by the qindex type") );

                // overflows size_t in single precision (and ergo also in double and quad)
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(62, mpi,gpu,omp), ContainsSubstring("memory would overflow size_t") );

                // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
                // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
                // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
                // advance or whether it proceeded to malloc() which subsequently failed
                #ifndef SANITIZER_IS_ACTIVE
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(50, mpi,gpu,omp), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
                #endif
            }
        }

        SECTION( "unsupported deployments" ) {

            if (!getQuESTEnv().isDistributed)
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(1, 1,0,0), ContainsSubstring("non-distributed QuEST environment") );

            if (!getQuESTEnv().isGpuAccelerated)
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(1, 0,1,0), ContainsSubstring("non-GPU-accelerated QuEST environment") );

            if (!getQuESTEnv().isMultithreaded)
                REQUIRE_THROWS_WITH( createCustomFullStateDiagMatr(1, 0,0,1), ContainsSubstring("non-multithreaded QuEST environment") );
        }
    }
}


TEST_CASE( "destroyCompMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        CompMatr m = createCompMatr(5);
        REQUIRE_NOTHROW( destroyCompMatr(m) );
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            CompMatr m;
            REQUIRE_THROWS_WITH( destroyCompMatr(m), ContainsSubstring("Invalid CompMatr") && ContainsSubstring("not created") );
        }
        #endif
    }
}


TEST_CASE( "destroyDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        DiagMatr m = createDiagMatr(5);
        REQUIRE_NOTHROW( destroyDiagMatr(m) );
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            DiagMatr m;
            REQUIRE_THROWS_WITH( destroyDiagMatr(m), ContainsSubstring("Invalid DiagMatr") && ContainsSubstring("not created") );
        }
        #endif
    }
}


TEST_CASE( "destroyFullStateDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        FullStateDiagMatr m = createFullStateDiagMatr(5);
        REQUIRE_NOTHROW( destroyFullStateDiagMatr(m) );
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            FullStateDiagMatr m;
            REQUIRE_THROWS_WITH( destroyFullStateDiagMatr(m), ContainsSubstring("Invalid FullStateDiagMatr") );
        }
        #endif
    }
}


TEST_CASE( "syncCompMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        CompMatr matr = createCompMatr(5);

        REQUIRE( *(matr.wasGpuSynced) == 0 );

        SECTION( "overwrites GPU elements" ) {

            // to test that the GPU memory was actually overwritten,
            // we would need a custom accessor of GPU memory, requiring
            // the tests are CUDA-compiled - no thank you mam! It is
            // certain this function works from the other GPU tests.

            SUCCEED( );
        }

        SECTION( "sets was-synced flag" ) {

            *(matr.wasGpuSynced) = 0;

            syncCompMatr(matr);
            REQUIRE( *(matr.wasGpuSynced) == 1 );
        }

        SECTION( "clears numerical flags" ) {

            *(matr).isApproxHermitian = 1;
            *(matr).isApproxUnitary   = 0;

            syncCompMatr(matr);
            REQUIRE( *(matr.isApproxHermitian) == -1 );
            REQUIRE( *(matr.isApproxUnitary)   == -1 );
        }

        destroyCompMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            CompMatr m;
            REQUIRE_THROWS_WITH( syncCompMatr(m), ContainsSubstring("Invalid CompMatr") && ContainsSubstring("not created") );
        }
        #endif
    }
}


TEST_CASE( "syncDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        DiagMatr matr = createDiagMatr(5);
        REQUIRE( *(matr.wasGpuSynced) == 0 );

        SECTION( "overwrites GPU elements" ) {

            // to test that the GPU memory was actually overwritten,
            // we would need a custom accessor of GPU memory, requiring
            // the tests are CUDA-compiled - no thank you mam! It is
            // certain this function works from the other GPU tests.

            SUCCEED( );
        }

        SECTION( "sets was-synced flag" ) {

            *(matr.wasGpuSynced) = 0;

            syncDiagMatr(matr);
            REQUIRE( *(matr.wasGpuSynced) == 1 );
        }

        SECTION( "clears numerical flags" ) {

            *(matr).isApproxHermitian = 1;
            *(matr).isApproxUnitary   = 0;
            *(matr).isApproxNonZero   = 1;
            *(matr).isStrictlyNonNegative = 0;

            syncDiagMatr(matr);
            REQUIRE( *(matr.isApproxHermitian) == -1 );
            REQUIRE( *(matr.isApproxUnitary)   == -1 );
            REQUIRE( *(matr.isApproxNonZero)   == -1 );
            REQUIRE( *(matr.isStrictlyNonNegative) == -1 );
        }

        destroyDiagMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            DiagMatr m;
            REQUIRE_THROWS_WITH( syncDiagMatr(m), ContainsSubstring("Invalid DiagMatr") && ContainsSubstring("not created") );
        }
        #endif
    }
}


TEST_CASE( "syncFullStateDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (auto& [label, matrix]: getCachedFullStateDiagMatrs()) {

            DYNAMIC_SECTION( label ) {

                *(matrix.wasGpuSynced) = 0;
                *(matrix).isApproxHermitian = 1;
                *(matrix).isApproxUnitary   = 0;
                *(matrix).isApproxNonZero   = 1;
                *(matrix).isStrictlyNonNegative = 0;

                syncFullStateDiagMatr(matrix);
                REQUIRE( *(matrix.wasGpuSynced)      == 1 );
                
                REQUIRE( *(matrix.isApproxHermitian) == -1 );
                REQUIRE( *(matrix.isApproxUnitary)   == -1 );
                REQUIRE( *(matrix.isApproxNonZero)   == -1 );
                REQUIRE( *(matrix.isStrictlyNonNegative) == -1 );

                // to test that the GPU memory was actually overwritten,
                // we would need a custom accessor of GPU memory, requiring
                // the tests are CUDA-compiled - no thank you mam! It is
                // certain this function works from the other GPU tests.
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            FullStateDiagMatr m;
            REQUIRE_THROWS_WITH( syncFullStateDiagMatr(m), ContainsSubstring("Invalid FullStateDiagMatr") );
        }
        #endif
    }
}


TEST_CASE( "setCompMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,6) );
        CAPTURE( numQubits );

        CompMatr matr = createCompMatr(numQubits);
        REQUIRE( *(matr.wasGpuSynced) == 0 );

        int dim = getPow2(numQubits);
        qmatrix ref = getRandomMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // nested pointers
            qcomp** ptrs = (qcomp**) malloc(dim * sizeof *ptrs);
            for (int i=0; i<dim; i++) {
                ptrs[i] = (qcomp*) malloc(dim * sizeof **ptrs);
                for (int j=0; j<dim; j++)
                    ptrs[i][j] = ref[i][j];
            }
            setCompMatr(matr, ptrs);
            REQUIRE_AGREE( matr, ref );
            REQUIRE( *(matr.wasGpuSynced) == 1 );

            // cannot test 2D VLAs in this C++ file

            // cleanup
            for (int i=0; i<dim; i++)
                free(ptrs[i]);
            free(ptrs);
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // nested vectors
            setCompMatr( matr, getZeroMatrix(dim) ); // clear
            setCompMatr( matr, ref );
            REQUIRE_AGREE( matr, ref );
            REQUIRE( *(matr.wasGpuSynced) == 1 );
        }

        destroyCompMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        CompMatr matr = createCompMatr(1);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            CompMatr bad;
            qcomp** dummy;
            REQUIRE_THROWS_WITH( setCompMatr(bad, dummy), ContainsSubstring("Invalid CompMatr") && ContainsSubstring("not created") );
        }
        #endif
        #endif

        SECTION( "null pointer" ) {

            qcomp** ptr = nullptr;
            REQUIRE_THROWS_WITH( setCompMatr(matr, ptr), ContainsSubstring("was a null pointer") );

            qcomp* arr[1] = {nullptr};
            REQUIRE_THROWS_WITH( setCompMatr(matr, arr), ContainsSubstring("contained a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_NOTHROW( setCompMatr(matr, {{1,2},{3,4}}) );

            REQUIRE_THROWS_WITH( setCompMatr(matr,{{1,2}}),              ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( setCompMatr(matr, {{1,2},{3,4},{5,6}}), ContainsSubstring("Incompatible number of rows") );
        
            REQUIRE_THROWS_WITH( setCompMatr(matr, {{0,0},{0}}),     ContainsSubstring("incompatible number of elements") );
            REQUIRE_THROWS_WITH( setCompMatr(matr, {{0,0},{0,0,0}}), ContainsSubstring("incompatible number of elements") );
        }

        destroyCompMatr(matr);
    }
}


TEST_CASE( "setDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,6) );
        CAPTURE( numQubits );

        DiagMatr matr = createDiagMatr(numQubits);
        REQUIRE( *(matr.wasGpuSynced) == 0 );

        int dim = getPow2(numQubits);
        qmatrix ref = getRandomDiagonalMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // pointer
            qvector diags = getDiagonals(ref);
            setDiagMatr(matr, diags.data());
            REQUIRE_AGREE( matr, ref );
            REQUIRE( *(matr.wasGpuSynced) == 1 );
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // vector
            setDiagMatr(matr, getZeroVector(dim)); // clear
            setDiagMatr(matr, getDiagonals(ref));
            REQUIRE_AGREE( matr, ref );
            REQUIRE( *(matr.wasGpuSynced) == 1 );
        }

        destroyDiagMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        DiagMatr matr = createDiagMatr(1);

        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            DiagMatr bad;
            qcomp* dummy;
            REQUIRE_THROWS_WITH( setDiagMatr(bad, dummy), ContainsSubstring("Invalid DiagMatr") && ContainsSubstring("not created") );
        }
        #endif

        SECTION( "null pointer" ) {

            qcomp* ptr = nullptr;
            REQUIRE_THROWS_WITH( setDiagMatr(matr, ptr), ContainsSubstring("was a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_NOTHROW( setDiagMatr(matr, {1,2}) );

            REQUIRE_THROWS_WITH( setDiagMatr(matr, {1}),     ContainsSubstring("Incompatible number of elements") );
            REQUIRE_THROWS_WITH( setDiagMatr(matr, {1,2,3}), ContainsSubstring("Incompatible number of elements") );
        }

        destroyDiagMatr(matr);
    }
}


TEST_CASE( "setInlineCompMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        CompMatr matr = createCompMatr(1);
        REQUIRE( *(matr.wasGpuSynced) == 0 );

        setInlineCompMatr( matr, 1, {{1,2},{3,4}} );
        REQUIRE_AGREE( matr, {{1,2},{3,4}} );
        REQUIRE( *(matr.wasGpuSynced) == 1 );

        destroyCompMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        CompMatr matr = createCompMatr(1);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            CompMatr bad;
            REQUIRE_THROWS_WITH( setInlineCompMatr(bad, 1, {{1,2},{3,4}}), ContainsSubstring("Invalid CompMatr") && ContainsSubstring("not created") );
        }
        #endif
        #endif

        SECTION( "mismatching dimension" ) {

            REQUIRE_THROWS_WITH( setInlineCompMatr(matr, 2, {{1,2},{3,4}}), ContainsSubstring("declared number of qubits") && ContainsSubstring("differs") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_NOTHROW( setInlineCompMatr(matr, 1, {{1,2},{3,4}}) );

            REQUIRE_THROWS_WITH( setInlineCompMatr(matr, 1, {{1,2}}),             ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( setInlineCompMatr(matr, 1, {{1,2},{3,4},{5,6}}), ContainsSubstring("Incompatible number of rows") );

            REQUIRE_THROWS_WITH( setInlineCompMatr(matr, 1, {{1},{2}}),         ContainsSubstring("One or more rows contained an incompatible number of elements") );
            REQUIRE_THROWS_WITH( setInlineCompMatr(matr, 1, {{1,2,3},{4,5,6}}), ContainsSubstring("One or more rows contained an incompatible number of elements") );
        }

        destroyCompMatr(matr);
    }
}


TEST_CASE( "setInlineDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        DiagMatr matr = createDiagMatr(1);
        REQUIRE( *(matr.wasGpuSynced) == 0 );

        setInlineDiagMatr( matr, 1, {1,2_i} );
        REQUIRE_AGREE( matr, {{1,0},{0,2_i}} );
        REQUIRE( *(matr.wasGpuSynced) == 1 );

        destroyDiagMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        DiagMatr matr = createDiagMatr(1);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            DiagMatr bad;
            REQUIRE_THROWS_WITH( setInlineDiagMatr(bad, 1, {1,2}), ContainsSubstring("Invalid DiagMatr") && ContainsSubstring("not created") );
        }
        #endif
        #endif

        SECTION( "mismatching dimension" ) {

            REQUIRE_THROWS_WITH( setInlineDiagMatr(matr, 2, {1,2}), ContainsSubstring("declared number of qubits") && ContainsSubstring("differs") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_NOTHROW( setInlineDiagMatr(matr, 1, {1,2}) );

            REQUIRE_THROWS_WITH( setInlineDiagMatr(matr, 1, {1}),     ContainsSubstring("Incompatible number of elements") );
            REQUIRE_THROWS_WITH( setInlineDiagMatr(matr, 1, {1,2,3}), ContainsSubstring("Incompatible number of elements") );
        }

        destroyDiagMatr(matr);
    }
}


TEST_CASE( "createInlineCompMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        CompMatr matr = createInlineCompMatr(1, {{1,2},{3,4}});

        REQUIRE_AGREE( matr, {{1,2},{3,4}} );
        REQUIRE( *(matr.wasGpuSynced) == 1 );

        destroyCompMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createInlineCompMatr(numQubits, {{1}}), ContainsSubstring("must target one or more qubits") );
        }

        SECTION( "mismatching dimension" ) {

            REQUIRE_THROWS_WITH( createInlineCompMatr(2, {{1,2},{3,4}}), ContainsSubstring("Incompatible number of rows") );
        }
    }
}


TEST_CASE( "createInlineDiagMatr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        DiagMatr matr = createInlineDiagMatr(1, {1,2_i});

        REQUIRE_AGREE( matr, {{1,0},{0,2_i}} );
        REQUIRE( *(matr.wasGpuSynced) == 1 );

        destroyDiagMatr(matr);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createInlineDiagMatr(numQubits, {1}), ContainsSubstring("must target one or more qubits") );
        }

        SECTION( "mismatching dimension" ) {

            REQUIRE_THROWS_WITH( createInlineDiagMatr(2, {1,2}), ContainsSubstring("Incompatible number of elements") );
        }
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */


void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems);
void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, std::vector<qcomp> in);

void setInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems, std::vector<qcomp> in);


void setDiagMatrFromMultiVarFunc(DiagMatr out, qcomp (*func)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);

void setDiagMatrFromMultiDimLists(DiagMatr out, void* lists, int* numQubitsPerDim, int numDims);


FullStateDiagMatr createFullStateDiagMatrFromPauliStrSum(PauliStrSum in);

void setFullStateDiagMatrFromPauliStrSum(FullStateDiagMatr out, PauliStrSum in);

void setFullStateDiagMatrFromMultiVarFunc(FullStateDiagMatr out, qcomp (*func)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);

void setFullStateDiagMatrFromMultiDimLists(FullStateDiagMatr out, void* lists, int* numQubitsPerDim, int numDims);


void reportCompMatr1(CompMatr1 matrix);
void reportCompMatr2(CompMatr2 matrix);
void reportCompMatr(CompMatr matrix);
void reportDiagMatr1(DiagMatr1 matrix);
void reportDiagMatr2(DiagMatr2 matrix);
void reportDiagMatr(DiagMatr matrix);
void reportFullStateDiagMatr(FullStateDiagMatr matr);
