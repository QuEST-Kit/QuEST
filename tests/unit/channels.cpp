/** @file
 * Unit tests of the channels module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitchannels Channels
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/random.hpp"

#include <vector>

using Catch::Matchers::ContainsSubstring;
using std::vector;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[channels]"



/** 
 * TESTS
 * 
 * @ingroup unitchannels
 * @{
 */


TEST_CASE( "createKrausMap", TEST_CATEGORY ) {
 
    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,4) ); // 4qb superop = 8qb matrix = 16 statevec
        int numMatrices = GENERATE( 1, 2, 10 );
        CAPTURE( numQubits, numMatrices );

        KrausMap map = createKrausMap(numQubits, numMatrices);

        // verify dimensions
        REQUIRE( map.numQubits == numQubits );
        REQUIRE( map.numMatrices == numMatrices );
        REQUIRE( map.numRows == getPow2(numQubits) );
        
        // verify superoperator dimensions
        REQUIRE( map.superop.numQubits == numQubits );
        REQUIRE( map.superop.numRows == getPow2(2 * numQubits) );

        // verify default fields
        REQUIRE( *(map.isApproxCPTP) == -1 );
        REQUIRE( *(map.superop.wasGpuSynced) == 0 );

        // verify pointers
        REQUIRE( map.matrices != nullptr );
        REQUIRE( map.superop.cpuElems != nullptr );
        REQUIRE( map.superop.cpuElemsFlat != nullptr );
        if (getQuESTEnv().isGpuAccelerated)
            REQUIRE( map.superop.gpuElemsFlat != nullptr );
        else
        REQUIRE( map.superop.gpuElemsFlat == nullptr );

        // verify that all matrices default to zero
        bool isZero = true;
        for (qindex i=0; i<map.numMatrices && isZero; i++)
            for (qindex r=0; r<map.numRows && isZero; r++)
                for (qindex c=0; c<map.numRows && isZero; c++)
                    isZero = (map.matrices[i][r][c] == qcomp(0,0));
        REQUIRE( isZero );

        // verify superoperator defaults to zero
        isZero = true;
        for (qindex r=0; r<map.superop.numRows && isZero; r++)
            for (qindex c=0; c<map.superop.numRows && isZero; c++)
                isZero = (map.superop.cpuElems[r][c] == qcomp(0,0));
        REQUIRE( isZero );

        destroyKrausMap(map);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env is init" ) {

            // no way to test this
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createKrausMap(numQubits,1), ContainsSubstring("one or more qubits") );
        }

        SECTION( "too few operators" ) {

            int numOpers =  GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createKrausMap(1,numOpers), ContainsSubstring("strictly positive number of matrices") );
        }

        SECTION( "too many qubits" ) {

            REQUIRE_THROWS_WITH( createKrausMap(100,1), ContainsSubstring("can be addressed by the qindex type") );

            REQUIRE_THROWS_WITH( createKrausMap(15,1), ContainsSubstring("necessary memory") && ContainsSubstring("would overflow") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createKrausMap(12,1), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }

        SECTION( "too many operators" ) {

            /// @todo
            /// there is currently no validation for that the operator
            /// list was too big and would ergo overflow memory (since
            /// it's just a ridiculous, annoying scenario)

            SUCCEED( );
        }
    }
}


TEST_CASE( "destroyKrausMap", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        KrausMap map = createKrausMap(3, 3);
        REQUIRE_NOTHROW( destroyKrausMap(map) );
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer interferes with un-initialised struct values
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            KrausMap m;
            REQUIRE_THROWS_WITH( destroyKrausMap(m), ContainsSubstring("Invalid KrausMap") && ContainsSubstring("not created") );
        }
        #endif
    }
}


TEST_CASE( "syncKrausMap", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numTargs = GENERATE( range(1,4) ); // 4qb superop = 8qb matrix = 16qb statevec
        int numMatrs = GENERATE( 1, 2, 10 );
        CAPTURE( numTargs, numMatrs );

        KrausMap map = createKrausMap(numTargs, numMatrs);
        REQUIRE( *(map.superop.wasGpuSynced) == 0 );
        
        *(map.isApproxCPTP) = 0;

        // validate the fields are updated
        syncKrausMap(map);
        REQUIRE( *(map.superop.wasGpuSynced) == 1 );
        REQUIRE( *(map.isApproxCPTP) == -1 );

        // validate the superop gets inferred correctly from the kraus map
        auto matrices = getRandomKrausMap(numTargs, numMatrs);
        setKrausMap(map, matrices); // calls sync anyway
        syncKrausMap(map);
        REQUIRE_AGREE( map.superop, getSuperOperator(matrices) );

        // to test that the GPU memory was actually overwritten,
        // we would need a custom accessor of GPU memory, requiring
        // the tests are CUDA-compiled - no thank you mam! It is
        // certain this function works from the other GPU tests.

        destroyKrausMap(map);
    }

    SECTION( LABEL_VALIDATION ) {

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            KrausMap m;
            REQUIRE_THROWS_WITH( syncKrausMap(m), ContainsSubstring("Invalid KrausMap") && ContainsSubstring("not created") );
        }
        #endif
        #endif
    }
}


TEST_CASE( "setKrausMap", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,4) ); // 4qb superop = 8qb matrix = 16qb statevec
        int numOps = GENERATE( 1, 2, 10 );
        CAPTURE( numQubits, numOps );

        KrausMap map = createKrausMap(numQubits, numOps);
        auto matrices = getRandomKrausMap(numQubits, numOps);

        SECTION( LABEL_C_INTERFACE ) {

            qcomp*** ptrs = (qcomp***) malloc(numOps * sizeof *ptrs);
            for (int n=0; n<numOps; n++) {
                ptrs[n] = (qcomp**) malloc(map.numRows * sizeof **ptrs);
                for (int r=0; r<map.numRows; r++)
                    ptrs[n][r] = matrices[n][r].data();
            }

            setKrausMap(map, ptrs);
            REQUIRE( *(map.superop.wasGpuSynced) == 1 );
            REQUIRE_AGREE( map.superop, getSuperOperator(matrices) );

            for (int n=0; n<numOps; n++)
                free(ptrs[n]);
            free(ptrs);
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            *(map.superop.wasGpuSynced) = 0;

            setKrausMap(map, matrices);
            REQUIRE( *(map.superop.wasGpuSynced) == 1 );
            REQUIRE_AGREE( map.superop, getSuperOperator(matrices) );
        }

        // to test that the GPU memory was actually overwritten,
        // we would need a custom accessor of GPU memory, requiring
        // the tests are CUDA-compiled - no thank you mam! It is
        // certain this function works from the other GPU tests.

        destroyKrausMap(map);
    }

    SECTION( LABEL_VALIDATION ) {

        // only C++ interface is validated
        
        int numQubits = 3;
        int numOps = 3;
        KrausMap map = createKrausMap(numQubits, numOps);

        int err = GENERATE( -1, +1 );

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            KrausMap bad;
            REQUIRE_THROWS_WITH( setKrausMap(bad, getRandomKrausMap(numQubits, numOps)), ContainsSubstring("invalid") );
        }
        #endif
        #endif

        SECTION( "inconsistent dimensions" ) {

            REQUIRE_THROWS_WITH( setKrausMap(map, getRandomKrausMap(numQubits+err, numOps)), ContainsSubstring("dimension") );
        }

        SECTION( "inconsistent number of matrices" ) {

            REQUIRE_THROWS_WITH( setKrausMap(map, getRandomKrausMap(numQubits, numOps-err)), ContainsSubstring("number of matrices") );
        }

        destroyKrausMap(map);
    }
}


TEST_CASE( "setInlineKrausMap", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,4) ); // 4qb superop = 8qb matrix = 16qb statevec
        int numOps = GENERATE( 1, 2, 10 );
        CAPTURE( numQubits, numOps );

        KrausMap map = createKrausMap(numQubits, numOps);
        auto matrices = getRandomKrausMap(numQubits, numOps);

        // only the C++ interface can be tested

        SECTION( LABEL_CPP_INTERFACE ) {

            *(map.superop.wasGpuSynced) = 0;

            setInlineKrausMap(map, numQubits, numOps, matrices);
            REQUIRE( *(map.superop.wasGpuSynced) == 1 );
            REQUIRE_AGREE( map.superop, getSuperOperator(matrices) );

            /// @todo test GPU memory is overwritten
        }

        destroyKrausMap(map);
    }

    SECTION( LABEL_VALIDATION ) {

        // only C++ interface can be validated
        
        int numQubits = 3;
        int numOps = 3;
        KrausMap map = createKrausMap(numQubits, numOps);

        int err = GENERATE( -1, +1 );

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            KrausMap bad;
            REQUIRE_THROWS_WITH( setInlineKrausMap(bad, numQubits, numOps, getRandomKrausMap(numQubits, numOps)), ContainsSubstring("invalid") );
        }
        #endif
        #endif

        SECTION( "macro parameters" ) {

            // check macro parameters are consistent
            REQUIRE_THROWS_WITH( setInlineKrausMap(map, numQubits+err, numOps,     getRandomKrausMap(numQubits, numOps)), ContainsSubstring("number of Kraus operators") && ContainsSubstring("qubits") );
            REQUIRE_THROWS_WITH( setInlineKrausMap(map, numQubits,     numOps+err, getRandomKrausMap(numQubits, numOps)), ContainsSubstring("number of Kraus operators") && ContainsSubstring("qubits") );
        }

        SECTION( "dimensions" ) {

            // check lists are correctly sized
            REQUIRE_THROWS_WITH( setInlineKrausMap(map, numQubits, numOps, getRandomKrausMap(numQubits+err, numOps    )), ContainsSubstring("dimension") );
            REQUIRE_THROWS_WITH( setInlineKrausMap(map, numQubits, numOps, getRandomKrausMap(numQubits,     numOps+err)), ContainsSubstring("number of matrices") );
        }

        destroyKrausMap(map);
    }
}


TEST_CASE( "createInlineKrausMap", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,4) ); // 4qb superop = 8qb matrix = 16qb statevec
        int numOps = GENERATE( 1, 2, 10 );
        CAPTURE( numQubits, numOps );

        auto matrices = getRandomKrausMap(numQubits, numOps);

        // only the C++ interface can be tested

        SECTION( LABEL_CPP_INTERFACE ) {

            KrausMap map = createInlineKrausMap(numQubits, numOps, matrices);

            REQUIRE( *(map.superop.wasGpuSynced) == 1 );
            REQUIRE_AGREE( map.superop, getSuperOperator(matrices) );

            /// @todo test GPU memory is overwritten

            destroyKrausMap(map);
        }
    }

    SECTION( LABEL_VALIDATION ) {

        // only C++ interface can be validated

        SECTION( "env not created" ) {

            // no way to test
            SUCCEED( );
        }

        SECTION( "incomaptible number of matrices" ) {

            int err = GENERATE( -1, +1 );

            REQUIRE_THROWS_WITH( createInlineKrausMap(3, 3, getRandomKrausMap(3, 3+err)), ContainsSubstring("matrices") );
        }

        SECTION( "incompatible dimensions" ) {

            int err = GENERATE( -1, +1 );

            REQUIRE_THROWS_WITH( createInlineKrausMap(3, 3, getRandomKrausMap(3+err, 3)), ContainsSubstring("rows") );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createInlineKrausMap(numQubits,1,{{{}}}), ContainsSubstring("one or more qubits") );
        }

        SECTION( "too few operators" ) {

            int numOpers =  GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createInlineKrausMap(1,numOpers,{{{}}}), ContainsSubstring("strictly positive number of matrices") );
        }

        SECTION( "too many qubits" ) {

            REQUIRE_THROWS_WITH( createInlineKrausMap(100,1,{{{}}}), ContainsSubstring("can be addressed by the qindex type") );

            REQUIRE_THROWS_WITH( createInlineKrausMap(15,1,{{{}}}), ContainsSubstring("necessary memory") && ContainsSubstring("would overflow") );

            // cannot check when massive alloc run-time fails since the passed
            // vectors must be of the correct size - and ergo impossibly big!
        }

        SECTION( "too many operators" ) {

            /// @todo
            /// there is currently no validation for that the operator
            /// list was too big and would ergo overflow memory (since
            /// it's just a ridiculous, annoying scenario)

            SUCCEED( );
        }
    }
}



TEST_CASE( "createSuperOp", TEST_CATEGORY ) {
 
    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,4) ); // 4qb superop = 8qb matrix = 16 statevec
        CAPTURE( numQubits );

        SuperOp op = createSuperOp(numQubits);

        // verify dimensions
        REQUIRE( op.numQubits == numQubits );
        REQUIRE( op.numRows == getPow2(2*numQubits) );
        
        // verify default fields
        REQUIRE( *(op.wasGpuSynced) == 0 );

        // verify pointers
        REQUIRE( op.cpuElems != nullptr );
        REQUIRE( op.cpuElemsFlat != nullptr );
        if (getQuESTEnv().isGpuAccelerated)
            REQUIRE( op.gpuElemsFlat != nullptr );
        else
        REQUIRE( op.gpuElemsFlat == nullptr );

        // verify superoperator defaults to zero
        bool isZero = true;
        for (qindex r=0; r<op.numRows && isZero; r++)
            for (qindex c=0; c<op.numRows && isZero; c++)
                isZero = (op.cpuElems[r][c] == qcomp(0,0));
        REQUIRE( isZero );

        destroySuperOp(op);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env is init" ) {

            // no way to test this
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createSuperOp(numQubits), ContainsSubstring("one or more qubits") );
        }

        SECTION( "too many qubits" ) {

            REQUIRE_THROWS_WITH( createSuperOp(100), ContainsSubstring("qindex") );

            REQUIRE_THROWS_WITH( createSuperOp(15), ContainsSubstring("size_t") );

            // no overflows, but definitely exceeds local RAM and fails to allocate; frightens address sanitizer!
            // note the specific error message depends on the what backend the auto-deployer tried to use (e.g.
            // GPU-accel or distributed) and whether memory-probers realised there was insufficient memory in
            // advance or whether it proceeded to malloc() which subsequently failed
            #ifndef SANITIZER_IS_ACTIVE
            REQUIRE_THROWS_WITH( createSuperOp(12), ContainsSubstring("failed") || ContainsSubstring("insufficient available memory") || ContainsSubstring("available GPU memory") );
            #endif
        }
    }
}


TEST_CASE( "syncSuperOp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SuperOp op = createSuperOp(1);
        REQUIRE( *(op.wasGpuSynced) == 0 );

        syncSuperOp(op);
        REQUIRE( *(op.wasGpuSynced) == 1 );

        // to test that the GPU memory was actually overwritten,
        // we would need a custom accessor of GPU memory, requiring
        // the tests are CUDA-compiled - no thank you mam! It is
        // certain this function works from the other GPU tests.

        destroySuperOp(op);
    }

    SECTION( LABEL_VALIDATION ) {

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            SuperOp m;
            REQUIRE_THROWS_WITH( syncSuperOp(m), ContainsSubstring("invalid fields") );
        }
        #endif
        #endif
    }
}


TEST_CASE( "destroySuperOp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SuperOp m = createSuperOp(5);
        REQUIRE_NOTHROW( destroySuperOp(m) );
    }

    SECTION( LABEL_VALIDATION ) {

        // sanitizer interferes with un-initialised struct values
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            SuperOp m;
            REQUIRE_THROWS_WITH( destroySuperOp(m), ContainsSubstring("invalid fields") );
        }
        #endif
    }
}


TEST_CASE( "setSuperOp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = GENERATE( range(1,4) ); // 3qb superop = 6qb matrix = 12qb statevec
        CAPTURE( numQubits );

        SuperOp op = createSuperOp(numQubits);
        REQUIRE( *(op.wasGpuSynced) == 0 );

        int dim = getPow2(2 * numQubits);
        qmatrix ref = getRandomMatrix(dim);

        SECTION( LABEL_C_INTERFACE ) {

            // nested pointers
            qcomp** ptrs = (qcomp**) malloc(dim * sizeof *ptrs);
            for (int i=0; i<dim; i++) {
                ptrs[i] = (qcomp*) malloc(dim * sizeof **ptrs);
                for (int j=0; j<dim; j++)
                    ptrs[i][j] = ref[i][j];
            }
            setSuperOp(op, ptrs);
            REQUIRE_AGREE( op, ref );
            REQUIRE( *(op.wasGpuSynced) == 1 );

            // cannot test 2D VLAs in this C++ file

            // cleanup
            for (int i=0; i<dim; i++)
                free(ptrs[i]);
            free(ptrs);
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // nested vectors
            setSuperOp( op, getZeroMatrix(dim) ); // clear
            setSuperOp( op, ref );
            REQUIRE_AGREE( op, ref );
            REQUIRE( *(op.wasGpuSynced) == 1 );
        }

        destroySuperOp(op);
    }

    SECTION( LABEL_VALIDATION ) {

        SuperOp op = createSuperOp(1);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            SuperOp bad;
            qcomp** dummy;
            REQUIRE_THROWS_WITH( setSuperOp(bad, dummy), ContainsSubstring("invalid fields") );
        }
        #endif
        #endif

        SECTION( "null pointer" ) {

            qcomp** ptr = nullptr;
            REQUIRE_THROWS_WITH( setSuperOp(op, ptr), ContainsSubstring("was a null pointer") );

            qcomp* arr[1] = {nullptr};
            REQUIRE_THROWS_WITH( setSuperOp(op, arr), ContainsSubstring("contained a null pointer") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_NOTHROW( setSuperOp(op, qmatrix(4, qvector(4))) );

            REQUIRE_THROWS_WITH( setSuperOp(op, qmatrix(3, qvector(4))), ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( setSuperOp(op, qmatrix(5, qvector(4))), ContainsSubstring("Incompatible number of rows") );
        
            REQUIRE_THROWS_WITH( setSuperOp(op, qmatrix(4, qvector(3))), ContainsSubstring("Incompatible number of columns") );
            REQUIRE_THROWS_WITH( setSuperOp(op, qmatrix(4, qvector(5))), ContainsSubstring("Incompatible number of columns") );
        }

        destroySuperOp(op);
    }
}


TEST_CASE( "setInlineSuperOp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SuperOp op = createSuperOp(1);
        REQUIRE( *(op.wasGpuSynced) == 0 );

        // can only check C++ interface

        setInlineSuperOp( op, 1, {{1,2,3,4},{5,6,7,8},{-9,-8,-7,-6},{5_i,4_i,3_i,2_i}} );
        REQUIRE_AGREE(    op,    {{1,2,3,4},{5,6,7,8},{-9,-8,-7,-6},{5_i,4_i,3_i,2_i}} );
        REQUIRE( *(op.wasGpuSynced) == 1 );

        destroySuperOp(op);
    }

    SECTION( LABEL_VALIDATION ) {

        SuperOp op = createSuperOp(1);

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            SuperOp bad;
            REQUIRE_THROWS_WITH( setInlineSuperOp(bad, 1, {{}}), ContainsSubstring("invalid fields") );
        }
        #endif
        #endif

        SECTION( "mismatching dimension" ) {

            REQUIRE_THROWS_WITH( setInlineSuperOp(op, 2, {{}}), ContainsSubstring("specified number of qubits") && ContainsSubstring("differs") );
        }

        SECTION( "invalid dimensions" ) {

            // detectable only by the C++ interface

            REQUIRE_NOTHROW( setInlineSuperOp(op, 1, qmatrix(4, qvector(4))) );

            REQUIRE_THROWS_WITH( setInlineSuperOp(op, 1, qmatrix(3, qvector(4))), ContainsSubstring("Incompatible number of rows") );
            REQUIRE_THROWS_WITH( setInlineSuperOp(op, 1, qmatrix(5, qvector(4))), ContainsSubstring("Incompatible number of rows") );

            REQUIRE_THROWS_WITH( setInlineSuperOp(op, 1, qmatrix(4, qvector(3))), ContainsSubstring("Incompatible number of columns") );
            REQUIRE_THROWS_WITH( setInlineSuperOp(op, 1, qmatrix(4, qvector(5))), ContainsSubstring("Incompatible number of columns") );
        }

        destroySuperOp(op);
    }
}


TEST_CASE( "createInlineSuperOp", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SuperOp op = createInlineSuperOp(1, {{1,2,3,4}, {5,6,7,8}, {-9,-8,-7,-6}, {-5_i,-4_i,-3_i,-2_i}});
        REQUIRE_AGREE( op,                  {{1,2,3,4}, {5,6,7,8}, {-9,-8,-7,-6}, {-5_i,-4_i,-3_i,-2_i}});
        REQUIRE( *(op.wasGpuSynced) == 1 );

        destroySuperOp(op);
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // impossible to test
            SUCCEED( );
        }

        SECTION( "too few qubits" ) {

            int numQubits = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createInlineSuperOp(numQubits, {{}}), ContainsSubstring("must act upon one or more qubits") );
        }

        SECTION( "mismatching dimension" ) {

            REQUIRE_THROWS_WITH( createInlineSuperOp(1, {{1,2},{3,4}}), ContainsSubstring("inconsistent with the number of rows") );
        }
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void reportKrausMap(KrausMap map);

void reportSuperOp(SuperOp op);
