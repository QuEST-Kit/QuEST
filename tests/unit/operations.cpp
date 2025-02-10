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
 * target qubits. This enum lists all template values.
 * Value 'anystates' is reserved for control qubits,
 * indicating when they must accompany ctrl-states in
 * the API signature.
 * 
 * For example:
 * - applyHadamard:                Ctrls=zero, Targs=one
 * - applyControlledSwap:          Ctrls=one,  Targs=two
 * - applyMultiControlledCompMatr: Ctrls=any,  Targs=any
 * - applyMultiStateControlledT:   Ctrls=anystates, Targs=one
 */

enum NumQubitsFlag { zero, one, two, any, anystates };

void assertNumQubitsFlagsAreValid(NumQubitsFlag ctrlsFlag, NumQubitsFlag targsFlag) {

    DEMAND( 
        ctrlsFlag == zero ||
        ctrlsFlag == one  ||
        ctrlsFlag == any  ||
        ctrlsFlag == anystates );

    DEMAND(
        targsFlag == one ||
        targsFlag == two ||
        targsFlag == any);
}

void assertNumQubitsFlagsValid(NumQubitsFlag ctrlsFlag, NumQubitsFlag targsFlag, vector<int> ctrls, vector<int> states, vector<int> targs) {

    assertNumQubitsFlagsAreValid(ctrlsFlag, targsFlag);

    DEMAND( targs.size() > 0 );

    if (targsFlag == one) 
        DEMAND( targs.size() == 1 );

    if (targsFlag == two)
        DEMAND( targs.size() == 2 );

    if (ctrlsFlag == zero)
        DEMAND( ctrls.size() == 0 );

    if (ctrlsFlag == one)
        DEMAND( ctrls.size() == 1 );

    if (ctrlsFlag == anystates)
        DEMAND( states.size() == ctrls.size() );
    else
        DEMAND( states.size() == 0 );
}


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

    assertNumQubitsFlagsAreValid(Ctrls, one);
    DEMAND( Ctrls != one || numFreeQubits >= 1 );
    DEMAND( numFreeQubits >= 0 );

    if constexpr (Ctrls == zero)
        return 0;
    
    if constexpr (Ctrls == one)
        return 1;

    if constexpr (Ctrls == any || Ctrls == anystates)
        return GENERATE_COPY( range(0, numFreeQubits) );
}

template <NumQubitsFlag Targs>
int GENERATE_NUM_TARGS(int numFreeQubits) {

    assertNumQubitsFlagsAreValid(zero, Targs);
    DEMAND( Targs != one || numFreeQubits >= 1 );
    DEMAND( Targs != two || numFreeQubits >= 2 );
    DEMAND( numFreeQubits > 0 );

    if constexpr (Targs == one)
        return 1;

    if constexpr (Targs == two)
        return 2;

    if constexpr (Targs == any)
        return GENERATE_COPY( range(1, numFreeQubits) );
}


/*
 * invoke an API operation (e.g. applyHadamard), passing
 * any elements of (ctrls,states,targs) that it accepts 
 * (as informed by the template values) along with any
 * additional arguments. We explicitly accept numCtrls
 * and numTargs, rather than inferring them from ctrls
 * and targs, such that invalid numbers (like -1) can be
 * passed by input validation testing.
 * 
 * This big, ugly bespoke function is necessary, rather
 * than a simple variadic template, because the QuEST 
 * API accepts fixed numbers of qubits as individual 
 * arguments, rather than as lists/vectors/pointers.
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, typename F, typename... Args>
void invokeApiOperation(
    F operation, Qureg qureg, 
    vector<int> ctrls, vector<int> states, int numCtrls, 
    vector<int> targs, int numTargs, 
    Args... args
) {
    assertNumQubitsFlagsValid(Ctrls, Targs, ctrls, states, targs);

    if constexpr (Ctrls == zero) {
        if constexpr (Targs == one) operation(qureg, targs[0], args...);
        if constexpr (Targs == two) operation(qureg, targs[0], targs[1], args...);
        if constexpr (Targs == any) operation(qureg, targs.data(), numTargs, args...);
    }
    if constexpr (Ctrls == one) {
        if constexpr (Targs == one) operation(qureg, ctrls[0], targs[0], args...);
        if constexpr (Targs == two) operation(qureg, ctrls[0], targs[0], targs[1], args...);
        if constexpr (Targs == any) operation(qureg, ctrls[0], targs.data(), numTargs, args...);
    }
    if constexpr (Ctrls == any) {
        if constexpr (Targs == one) operation(qureg, ctrls.data(), numCtrls, targs[0], args...);
        if constexpr (Targs == two) operation(qureg, ctrls.data(), numCtrls, targs[0], targs[1], args...);
        if constexpr (Targs == any) operation(qureg, ctrls.data(), numCtrls, targs.data(), numTargs, args...);
    }
    if constexpr (Ctrls == anystates) {
        if constexpr (Targs == one) operation(qureg, ctrls.data(), states.data(), numCtrls, targs[0], args...);
        if constexpr (Targs == two) operation(qureg, ctrls.data(), states.data(), numCtrls, targs[0], targs[1], args...);
        if constexpr (Targs == any) operation(qureg, ctrls.data(), states.data(), numCtrls, targs.data(), numTargs, args...);
    }
}

// overload to avoid passing numCtrls and numTargs

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, typename F, typename... Args>
void invokeApiOperation(F operation, Qureg qureg, vector<int> ctrls, vector<int> states, vector<int> targs, Args... args) {
    invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, ctrls.size(), targs, targs.size(), args...);
}


/*
 * perform the unit test for an API operation which
 * is "fixed", i.e. effects a fixed-size (1 or 2 target) 
 * non-parameterised matrix upon the state. The control
 * qubits are completely free, so this function can be
 * invoked upon all many-controlled variants of the API.
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, typename F>
void performTestUponFixedOperation(F operation, qmatrix matrix) {
    
    assertNumQubitsFlagsAreValid(Ctrls, Targs);

    auto statevecQuregs = getCachedStatevecs();
    auto densmatrQuregs = getCachedDensmatrs();

    SECTION( "correctness" ) {

        qvector statevecRef = getZeroVector(getPow2(NUM_QUREG_QUBITS));
        qmatrix densmatrRef = getZeroMatrix(getPow2(NUM_QUREG_QUBITS));

        int numCtrls = GENERATE_NUM_CTRLS<Ctrls>(NUM_QUREG_QUBITS);
        int numTargs = GENERATE_NUM_TARGS<Targs>(NUM_QUREG_QUBITS - numCtrls); // always 1 or 2

        auto listpair = GENERATE_COPY( disjointsublists(range(0,NUM_QUREG_QUBITS), numCtrls, numTargs) );
        auto ctrls = std::get<0>(listpair);
        auto targs = std::get<1>(listpair);
        auto states = getRandomInts(0, 2, numCtrls); // may be ignored

        CAPTURE( ctrls, targs, states );

        auto testFunc = [&]<class T>(Qureg qureg, T& reference) -> void { 
            invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, targs);
            applyReferenceOperator(reference, ctrls, states, targs, matrix);
        };

        SECTION( "statevector"    ) { performTestUponAllQuregDeployments(statevecQuregs, statevecRef, testFunc); }
        SECTION( "density matrix" ) { performTestUponAllQuregDeployments(densmatrQuregs, densmatrRef, testFunc); }
    }

    SECTION( "input validation" ) {

        Qureg qureg = statevecQuregs["CPU"];
        vector<int> ctrls = (Ctrls == one)? vector<int>{0} : vector<int>{};
        vector<int> targs = (Targs == one)? vector<int>{1} : vector<int>{1, 2};
        vector<int> states(0, ctrls.size());

        SECTION( "qureg initialisation" ) {

            Qureg uninit;

            REQUIRE_THROWS_WITH( 
                (invokeApiOperation<Ctrls,Targs>(operation, uninit, ctrls, states, targs)),
                contains("invalid Qureg")
            );
        }

        SECTION( "target index" ) {

            targs[0] = GENERATE( -1, NUM_QUREG_QUBITS );

            REQUIRE_THROWS_WITH( 
                (invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, targs)),
                contains("invalid target qubit")
            );
        }

        if constexpr (Ctrls == any) {

            SECTION( "number of controls" ) {

                int numCtrls = GENERATE( -1, NUM_QUREG_QUBITS+1 );

                REQUIRE_THROWS_WITH( 
                    (invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, numCtrls, targs, targs.size())),
                    contains("number of control qubits")
                );
            }

        }

        if constexpr (Ctrls != zero) {

            SECTION( "control indices" ) {

                ctrls = { GENERATE( -1, NUM_QUREG_QUBITS ) };

                REQUIRE_THROWS_WITH( 
                    (invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("invalid control qubit")
                );
            }

        }

        if constexpr (Ctrls == any) {

            SECTION( "control repetition" ) {

                ctrls = {0, 0};

                REQUIRE_THROWS_WITH( 
                    (invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("control qubits contained duplicates")
                );
            }

        }

        if constexpr (Ctrls != zero) {

            SECTION( "target in controls" ) {

                ctrls = {targs[0]};

                REQUIRE_THROWS_WITH( 
                    (invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, targs)),
                    contains("qubit appeared among both the control and target qubits")
                );
            }
        }

        if constexpr (Ctrls == anystates) {

            SECTION( "control state bits" ) {

                states = { GENERATE( -1, 2 ) };
                ctrls = {0};

                REQUIRE_THROWS_WITH( 
                    (invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, targs)),
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
        performTestUponFixedOperation<zero,numtargs>( apply ## name, matrix); } \
    TEST_CASE( "applyControlled" #name, "[operations]" ) { \
        performTestUponFixedOperation<one,numtargs>( applyControlled ## name, matrix); } \
    TEST_CASE( "applyMultiControlled" #name, "[operations]" ) { \
        performTestUponFixedOperation<any,numtargs>( applyMultiControlled ## name, matrix); } \
    TEST_CASE( "applyMultiStateControlled" #name, "[operations]" ) { \
        performTestUponFixedOperation<anystates,numtargs>( applyMultiStateControlled ## name, matrix); } 

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

