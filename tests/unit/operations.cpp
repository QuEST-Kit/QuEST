#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/cache.hpp"
#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"



// TODO:
// some of the below functions may move into utils/

auto contains(std::string str) {
    
    return Catch::Matchers::ContainsSubstring(str, Catch::CaseSensitive::No);
}


/*
 * execute 'function' upon each kind of cached qureg
 * (e.g. distributed, GPU-accelerated, etc) and a
 * reference state (T1 = qvector or qmatrix), when
 * both are initialised in the debug state, and
 * thereafter confirm they approximately agree. Each
 * qureg deployment is featured in a separate test
 * section, so are accounted distinctly.
 */

template <class T1, class T2>
void performTestUponAllQuregDeployments(quregCache& quregs, T1& reference, T2& function) {

    for (auto& [label, qureg]: quregs) {

        DYNAMIC_SECTION( label ) {

            initDebugState(qureg);
            setToDebugState(reference);

            function(qureg, reference);
            REQUIRE_AGREE( qureg, reference );
        }
    }
}


/*
 * we will define templated functions for processing
 * many API signatures which accept zero, one, two
 * or any number (including 0-2) of control and/or
 * target qubits. This enum lists all template values
 */

enum NumQubitsFlag { zero, one, two, any };


/*
 * extract the runtime values of the number of control
 * and target qubtits from their compile-time templates.
 * When their number is permitted to be 'any', we use
 * a generator to successively test all possible numbers.
 * As such, this function is named similar to a Catch2
 * macro so the caller recognises it is a generator.
 */

template <NumQubitsFlag Ctrls>
int GENERATE_NUM_CTRLS(int numFreeQubits) {
    
    DEMAND( Ctrls != two );
    DEMAND( Ctrls != one || numFreeQubits > 0 );

    if constexpr (Ctrls == zero)
        return 0;
    
    if constexpr (Ctrls == one)
        return 1;

    if constexpr (Ctrls == any)
        return GENERATE_COPY( range(0, numFreeQubits) );
}

template <NumQubitsFlag Targs>
int GENERATE_NUM_TARGS(int numFreeQubits) {

    DEMAND( Targs != zero );
    DEMAND( Targs != two || numFreeQubits >= 2 );
    DEMAND( numFreeQubits > 0 );

    if constexpr (Targs == one)
        return 1;

    if constexpr (Targs == two)
        return 2;

    if constexpr(Targs == any)
        return GENERATE_COPY( range(1, numFreeQubits) );
}


/*
 * invoke an API operation (e.g. applyHadamard), passing
 * any args of (ctrls,states,targs) that it accepts, as
 * informed by the template values. We permit explicitly
 * passing numCtrls distinctly from ctrls.size(), so that
 * we input validation tests can pass invalid numCtrls.
 *
 * 'FixedOperation' refers to API operations which effect
 * fixed-size (i.e. 1 or 2 target) non-parameterised
 * matrices upon the state, such as Hadamard or SWAPs.
 */

template <NumQubitsFlag Ctrls, bool HasCtrlStates, NumQubitsFlag Targs, typename F>
void invokeFixedOperation(F operation, Qureg qureg, vector<int> ctrls, vector<int> states, int numCtrls, vector<int> targs) {

    DEMAND( Ctrls != two );
    DEMAND( Targs == one || Targs == two );
    if (HasCtrlStates)
        DEMAND( Ctrls == any );
    if (Ctrls == zero)
        DEMAND( ctrls.size() == 0 );
    if (Ctrls == one)
        DEMAND( ctrls.size() == 1 );

    if constexpr (Ctrls == zero) {
        if constexpr (Targs == one) operation(qureg, targs[0]);
        if constexpr (Targs == two) operation(qureg, targs[0], targs[1]);
    }
    if constexpr (Ctrls == one) {
        if constexpr (Targs == one) operation(qureg, ctrls[0], targs[0]);
        if constexpr (Targs == two) operation(qureg, ctrls[0], targs[0], targs[1]);
    }
    if constexpr (Ctrls == any && ! HasCtrlStates) {
        if constexpr (Targs == one) operation(qureg, ctrls.data(), numCtrls, targs[0]);
        if constexpr (Targs == two) operation(qureg, ctrls.data(), numCtrls, targs[0], targs[1]);
    }
    if constexpr (Ctrls == any && HasCtrlStates) {
        if constexpr (Targs == one) operation(qureg, ctrls.data(), states.data(), numCtrls, targs[0]);
        if constexpr (Targs == two) operation(qureg, ctrls.data(), states.data(), numCtrls, targs[0], targs[1]);
    }
}

template <NumQubitsFlag Ctrls, bool HasCtrlStates, NumQubitsFlag Targs, typename F>
void invokeFixedOperation(F operation, Qureg qureg, vector<int> ctrls, vector<int> states, vector<int> targs) {

    invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, ctrls.size(), targs);
}


/*
 * perform the unit test for an API operation which
 * is "fixed", i.e. effects a fixed-size (1 or 2 target) 
 * non-parameterised matrix upon the state. The control
 * qubits are completely free, so this function can be
 * invoked upon all many-controlled variants of the API.
 */

template <NumQubitsFlag Ctrls, bool HasCtrlStates, NumQubitsFlag Targs, typename F>
void performTestUponFixedOperation(F operation, qmatrix matrix) {
    
    DEMAND( Ctrls != two );
    DEMAND( Targs == one  || Targs == two );
    DEMAND( Ctrls == any  || ! HasCtrlStates );

    auto statevecQuregs = getCachedStatevecs();
    auto densmatrQuregs = getCachedDensmatrs();

    SECTION( "correctness" ) {

        qvector statevecRef = getZeroVector(getPow2(NUM_QUREG_QUBITS));
        qmatrix densmatrRef = getZeroMatrix(getPow2(NUM_QUREG_QUBITS));

        int numCtrls = GENERATE_NUM_CTRLS<Ctrls>(NUM_QUREG_QUBITS);
        int numTargs = GENERATE_NUM_TARGS<Targs>(NUM_QUREG_QUBITS - numCtrls);

        auto listpair = GENERATE_COPY( disjointsublists(range(0,NUM_QUREG_QUBITS), numCtrls, numTargs) );
        auto ctrls = std::get<0>(listpair);
        auto targs = std::get<1>(listpair);
        auto states = (HasCtrlStates)? getRandomInts(0, 2, numCtrls) : vector<int>();

        CAPTURE( ctrls, targs, states );

        auto testFunc = [&]<class T>(Qureg qureg, T& reference) -> void { 
            invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, targs);
            applyReferenceOperator(reference, ctrls, states, targs, matrix);
        };

        SECTION( "statevector"    ) performTestUponAllQuregDeployments(statevecQuregs, statevecRef, testFunc);
        SECTION( "density matrix" ) performTestUponAllQuregDeployments(densmatrQuregs, densmatrRef, testFunc);
    }

    SECTION( "input validation" ) {

        Qureg qureg = statevecQuregs["CPU"];
        vector<int> ctrls = (Ctrls == one)? vector<int>{0} : vector<int>{};
        vector<int> targs = {1, 2};
        vector<int> states = {};

        SECTION( "qureg initialisation" ) {

            Qureg uninit;

            REQUIRE_THROWS_WITH( 
                (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, uninit, ctrls, states, targs)),
                contains("invalid Qureg")
            );
        }

        SECTION( "target index" ) {

            targs[0] = GENERATE( -1, NUM_QUREG_QUBITS );

            REQUIRE_THROWS_WITH( 
                (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, targs)),
                contains("invalid target qubit")
            );
        }

        if constexpr (Ctrls == any) {

            SECTION( "number of controls" ) {

                int numCtrls = GENERATE( -1, NUM_QUREG_QUBITS+1 );

                REQUIRE_THROWS_WITH( 
                    (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, numCtrls, targs)),
                    contains("number of control qubits")
                );
            }

        }

        if constexpr (Ctrls != zero) {

            SECTION( "control indices" ) {

                ctrls = { GENERATE( -1, NUM_QUREG_QUBITS ) };

                REQUIRE_THROWS_WITH( 
                    (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("invalid control qubit")
                );
            }

        }

        if constexpr (Ctrls == any) {

            SECTION( "control repetition" ) {

                ctrls = {0, 0};

                REQUIRE_THROWS_WITH( 
                    (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("control qubits contained duplicates")
                );
            }

        }

        if constexpr (Ctrls != zero) {

            SECTION( "target in controls" ) {

                ctrls = {targs[0]};

                REQUIRE_THROWS_WITH( 
                    (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("qubit appeared among both the control and target qubits")
                );
            }
        }

        if constexpr (HasCtrlStates) {

            SECTION( "control state bits" ) {

                states = { GENERATE( -1, 2 ) };
                ctrls = {0};

                REQUIRE_THROWS_WITH( 
                    (invokeFixedOperation<Ctrls,HasCtrlStates,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("invalid control-state")
                );
            }
        }
    }
}


/*
 * the Z-basis matrices of all API operations for
 * which the effected matrix is fixed-size (1 or 2
 * targets) and non-parameterised, as accepted by
 * function performTestUponFixedOperation()
 */

namespace FixedMatrices {

    qmatrix H = {
        {1/sqrt(2),  1/sqrt(2)},
        {1/sqrt(2), -1/sqrt(2)}};

    qmatrix X = {
        {0, 1},
        {1, 0}};
    
    qmatrix Y = {
        {0, -1_i},
        {1_i, 0}};

    qmatrix Z = {
        {1,  0},
        {0, -1}};

    qmatrix T = {
        {1, 0},
        {0, exp(1_i*M_PI/4)}};

    qmatrix S = {
        {1, 0},
        {0, 1_i}};

    qmatrix SWAP = {
        {1, 0, 0, 0},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1}};

    qmatrix sqrtSWAP = {
        {1, 0, 0, 0},
        {0, (1+1_i)/2, (1-1_i)/2, 0},
        {0, (1-1_i)/2, (1+1_i)/2, 0},
        {0, 0, 0, 1}};
}


/*
 * instantiate unit tests for all ctrl-variants
 * of a "fixed" API operation
 */

#define TEST_FIXED_ANY_CTRL_OPERATION( name, numtargs, matrix ) \
    TEST_CASE( "apply" #name, "[operations]" ) { \
        performTestUponFixedOperation<zero,false,numtargs>( apply ## name, matrix); } \
    TEST_CASE( "applyControlled" #name, "[operations]" ) { \
        performTestUponFixedOperation<one,false,numtargs>( applyControlled ## name, matrix); } \
    TEST_CASE( "applyMultiControlled" #name, "[operations]" ) { \
        performTestUponFixedOperation<any,false,numtargs>( applyMultiControlled ## name, matrix); } \
    TEST_CASE( "applyMultiStateControlled" #name, "[operations]" ) { \
        performTestUponFixedOperation<any,true,numtargs>( applyMultiStateControlled ## name, matrix); } 

TEST_FIXED_ANY_CTRL_OPERATION( Hadamard, one, FixedMatrices::H );
TEST_FIXED_ANY_CTRL_OPERATION( PauliX,   one, FixedMatrices::X );
TEST_FIXED_ANY_CTRL_OPERATION( PauliY,   one, FixedMatrices::Y );
TEST_FIXED_ANY_CTRL_OPERATION( PauliZ,   one, FixedMatrices::Z );
TEST_FIXED_ANY_CTRL_OPERATION( T,        one, FixedMatrices::T );
TEST_FIXED_ANY_CTRL_OPERATION( S,        one, FixedMatrices::S );
TEST_FIXED_ANY_CTRL_OPERATION( Swap,     two, FixedMatrices::SWAP );
TEST_FIXED_ANY_CTRL_OPERATION( SqrtSwap, two, FixedMatrices::sqrtSWAP );


// TODO:
// we will likely change the above design, for
// example modularising the input validation 
// checks, when extending to non-fixed operations

