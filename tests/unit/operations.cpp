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

#include <tuple>

using std::tuple;



// TODO:
// some of the below functions may move into utils/

auto contains(std::string str) {
    
    return Catch::Matchers::ContainsSubstring(str, Catch::CaseSensitive::No);
}



/*
 * reference operator matrices used by testing
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

namespace SinglyParameterisedMatrices {

    auto Rx = [](qreal p) { return getExponentialOfPauliMatrix(p, FixedMatrices::X); };
    auto Ry = [](qreal p) { return getExponentialOfPauliMatrix(p, FixedMatrices::Y); };
    auto Rz = [](qreal p) { return getExponentialOfPauliMatrix(p, FixedMatrices::Z); };
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

void testAllQuregDeployments(quregCache& quregs, auto& reference, auto& function) {

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
 * Template flags for specifying what kind of additional 
 * arguments (in addition to ctrls/states/targs below) are
 * accepted by an API operation, when passing said operation
 * to automated testing facilities. This is NOT consulted
 * when generically invoking the API operation (we instead
 * used variadic templates above for that), but IS used
 * by the testing code to decide how to prepare inputs.
 * 
 * For example:
 * - applyHadamard:           none
 * - applyRotateX:            scalar
 * - applyDiagMatr1:          diagmatr  (Targs=one)
 * - applyDiagMatrPower:      diagpower (Targs=any)
 * - applyControlledCompMatr: compmatr  (Targs=any)
*/

enum ArgsFlag { none, scalar, diagmatr, diagpower, compmatr };


/*
 * Template flags for specifying how many control and
 * target qubits are accepted by an API operation.
 * Value 'anystates' is reserved for control qubits,
 * indicating that ctrls must accompany ctrl-states in
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

void assertNumQubitsFlagsValid(
    NumQubitsFlag ctrlsFlag, NumQubitsFlag targsFlag, 
    vector<int> ctrls, vector<int> states, vector<int> targs
) {
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

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
int GENERATE_NUM_TARGS(int numQuregQubits) {

    assertNumQubitsFlagsAreValid(zero, Targs);
    DEMAND( Targs != one || numQuregQubits >= 1 );
    DEMAND( Targs != two || numQuregQubits >= 2 );
    DEMAND( numQuregQubits > 0 );

    // single choice if #targs is compile-time set
    if constexpr (Targs == one)
        return 1;
    if constexpr (Targs == two)
        return 2;

    if constexpr (Targs == any) {

        // we can target all qubits...
        int maxNumTargs = numQuregQubits;

        // unless we're applying CompMatr in distributed
        // mode. Technically we can support all targets
        // if the Qureg is a density matrix or not
        // distributed, but we safely choose the min here
        if (Args == compmatr)
            maxNumTargs = numQuregQubits - getLog2(getQuESTEnv().numNodes);

        // we must also ensure there is space left for a forced ctrl
        if (Ctrls == one && maxNumTargs == numQuregQubits)
            maxNumTargs -= 1;

        return GENERATE_COPY( range(1, maxNumTargs+1) );
    }
}

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
        return GENERATE_COPY( range(0, numFreeQubits+1) );
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
 * arguments, rather than as lists/vectors/pointers. Note
 * that our use of variadic templates (args) means we do
 * not need to include ArgsFlag as a template parameter.
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs>
void invokeApiOperation(
    auto operation, Qureg qureg, 
    vector<int> ctrls, vector<int> states, int numCtrls, 
    vector<int> targs, int numTargs, 
    auto&... args
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

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs>
void invokeApiOperation(auto operation, Qureg qureg, vector<int> ctrls, vector<int> states, vector<int> targs, auto&... args) {
    invokeApiOperation<Ctrls,Targs>(operation, qureg, ctrls, states, ctrls.size(), targs, targs.size(), args...);
}


/*
 * prepare an API matrix (e.g. CompMatr1), as per
 * the given template parameters, with random
 * elements. This is used for testing API functions
 * which accept matrices.
 */

template <NumQubitsFlag Targs, ArgsFlag Args>
auto getRandomApiMatrix(int numTargs) {

    DEMAND(
        Args == diagmatr ||
        Args == diagpower ||
        Args == compmatr );
    
    qmatrix qm = (Args == compmatr)?
        getRandomUnitary(numTargs) : 
        getRandomDiagonalUnitary(numTargs);

    if constexpr (Args == compmatr && Targs == one) 
        return getCompMatr1(qm);

    if constexpr (Args == compmatr && Targs == two)
        return getCompMatr2(qm);

    if constexpr (Args == compmatr && Targs == any) {
        CompMatr cm = createCompMatr(numTargs); // leaks
        setCompMatr(cm, qm);
        return cm;
    }

    qvector dv = getDiagonals(qm);
    constexpr bool diag = (Args == diagmatr || Args == diagpower);

    if constexpr (diag && Targs == one)
        return getDiagMatr1(dv);

    if constexpr (diag && Targs == two)
        return getDiagMatr2(dv);

    if constexpr (diag && Targs == any) {
        DiagMatr dm = createDiagMatr(numTargs); // leaks
        setDiagMatr(dm, dv);
        return dm;
    }
}


/*
 * chooses random values for the remaining arguments
 * (after controls/states/targets) to API operations,
 * with types informed by the template parameters.
 * For example:
 * - applyRotateX() accepts a scalar
 * - applyCompMatr() accepts a CompMatr
 * - applyDiagMatrPower() accepts a DiagMatr and qcomp
 */

template <NumQubitsFlag Targs, ArgsFlag Args>
auto getRandomRemainingArgs(int numTargs) {

    if constexpr (Args == none)
        return tuple{};

    if constexpr (Args == scalar)
        return tuple{ getRandomReal(-2*M_PI, 2*M_PI) };

    if constexpr (Args == compmatr || Args == diagmatr)
        return tuple{ getRandomApiMatrix<Targs,Args>(numTargs) };

    if constexpr (Args == diagpower)
        return tuple{ getRandomApiMatrix<Targs,Args>(numTargs), getRandomComplex() };
}


/*
 * unpack the given reference operator matrix (a qmatrix) 
 * for an API operation, which is passed to testOperation(),
 * and which will be effected upon the reference state (a 
 * qvector or qmatrix). The type/form of matrixRefGen depends 
 * on the type of API operation, indicated by template parameter.
 */

template <ArgsFlag Args>
auto getReferenceMatrix(auto matrixRefGen, int numTargs, auto additionalArgs) {

    if constexpr (Args == none)
        return matrixRefGen;

    if constexpr (Args == scalar)
        return matrixRefGen(std::get<0>(additionalArgs));

    if constexpr (Args == compmatr || Args == diagmatr)
        return getMatrix(std::get<0>(additionalArgs));

    if constexpr (Args == diagpower) {
        qmatrix diag = getMatrix(std::get<0>(additionalArgs));
        qcomp power = std::get<1>(additionalArgs);
        return getPowerOfDiagonalMatrix(diag, power);
    }
}


/*
 * perform a unit test on an API operation. The
 * template parameters are compile-time clues
 * about what inputs to prepare and pass to the
 * operation, and how its reference matrix (arg
 * matrixRefGen) is formatted.
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperation(auto operation, auto matrixRefGen) {

    assertNumQubitsFlagsAreValid(Ctrls, Targs);

    // use existing cached Quregs
    auto statevecQuregs = getCachedStatevecs();
    auto densmatrQuregs = getCachedDensmatrs();

    SECTION( "correctness" ) {

        qvector statevecRef = getZeroVector(getPow2(NUM_QUREG_QUBITS));
        qmatrix densmatrRef = getZeroMatrix(getPow2(NUM_QUREG_QUBITS));

        // try all possible number of ctrls and targs
        int numTargs = GENERATE_NUM_TARGS<Ctrls,Targs,Args>(NUM_QUREG_QUBITS);
        int numCtrls = GENERATE_NUM_CTRLS<Ctrls>(NUM_QUREG_QUBITS - numTargs);
        
        // try all possible ctrls and targs
        auto listpair = GENERATE_COPY( disjointsublists(range(0,NUM_QUREG_QUBITS), numCtrls, numTargs) );
        auto ctrls = std::get<0>(listpair);
        auto targs = std::get<1>(listpair);

        // randomise control states (if operation accepts them)
        auto states = getRandomInts(0, 2, numCtrls * (Ctrls == anystates));

        // randomise remaining operation parameters
        auto primaryArgs = tuple{ ctrls, states, targs };
        auto furtherArgs = getRandomRemainingArgs<Targs,Args>(numTargs);
        auto matrixRef = getReferenceMatrix<Args>(matrixRefGen, numTargs, furtherArgs);

        // prepare test function which will receive both statevectors and density matrices
        auto testFunc = [&](Qureg qureg, auto& stateRef) -> void { 
            
            // invoke API operation, passing all args (unpacking variadic)
            auto apiFunc = [](auto&&... args) { return invokeApiOperation<Ctrls,Targs>(args...); };
            auto allArgs = std::tuple_cat(tuple{operation, qureg}, primaryArgs, furtherArgs);
            std::apply(apiFunc, allArgs);

            // update reference state
            applyReferenceOperator(stateRef, ctrls, states, targs, matrixRef);
        };

        // report relevant inputs if subsequent tests fail
        CAPTURE( ctrls, targs, states );

        // test API operation on all available deployment combinations (e.g. MPI, MPI+CUDA, etc)
        SECTION( "statevector"    ) { testAllQuregDeployments(statevecQuregs, statevecRef, testFunc); }
        SECTION( "density matrix" ) { testAllQuregDeployments(densmatrQuregs, densmatrRef, testFunc); }
    }
}


/*
 * perform unit tests for the four distinctly-controlled
 * variants of the given API operation (with the specified
 * function name suffix). 'numtargs' indicates the number 
 * of target qubits accepted by the operation, and 'argtype'
 * indicates the types of remaining arguments (if any exist).
 * 'matrixgen' is the matrix representation (of varying 
 * formats) of the operation, against which it will be compared.
 */

#define TEST_ANY_CTRL_OPERATION( namesuffix, numtargs, argtype, matrixgen ) \
    TEST_CASE( "apply" #namesuffix,                     "[operations]" ) { testOperation<zero,     numtargs,argtype>( apply ## namesuffix,                     matrixgen); } \
    TEST_CASE( "applyControlled" #namesuffix,           "[operations]" ) { testOperation<one,      numtargs,argtype>( applyControlled ## namesuffix,           matrixgen); } \
    TEST_CASE( "applyMultiControlled" #namesuffix,      "[operations]" ) { testOperation<any,      numtargs,argtype>( applyMultiControlled ## namesuffix,      matrixgen); } \
    TEST_CASE( "applyMultiStateControlled" #namesuffix, "[operations]" ) { testOperation<anystates,numtargs,argtype>( applyMultiStateControlled ## namesuffix, matrixgen); } 


TEST_ANY_CTRL_OPERATION( CompMatr1, one, compmatr, nullptr );
TEST_ANY_CTRL_OPERATION( CompMatr2, two, compmatr, nullptr );
TEST_ANY_CTRL_OPERATION( CompMatr,  any, compmatr, nullptr );
TEST_ANY_CTRL_OPERATION( DiagMatr1, one, diagmatr, nullptr );
TEST_ANY_CTRL_OPERATION( DiagMatr2, two, diagmatr, nullptr );
TEST_ANY_CTRL_OPERATION( DiagMatr,  any, diagmatr, nullptr );

TEST_ANY_CTRL_OPERATION( DiagMatrPower,  any, diagpower, nullptr );

TEST_ANY_CTRL_OPERATION( Hadamard, one, none, FixedMatrices::H );
TEST_ANY_CTRL_OPERATION( PauliX,   one, none, FixedMatrices::X );
TEST_ANY_CTRL_OPERATION( PauliY,   one, none, FixedMatrices::Y );
TEST_ANY_CTRL_OPERATION( PauliZ,   one, none, FixedMatrices::Z );
TEST_ANY_CTRL_OPERATION( T,        one, none, FixedMatrices::T );
TEST_ANY_CTRL_OPERATION( S,        one, none, FixedMatrices::S );
TEST_ANY_CTRL_OPERATION( Swap,     two, none, FixedMatrices::SWAP );
TEST_ANY_CTRL_OPERATION( SqrtSwap, two, none, FixedMatrices::sqrtSWAP );

TEST_ANY_CTRL_OPERATION( RotateX, one, scalar, SinglyParameterisedMatrices::Rx );
TEST_ANY_CTRL_OPERATION( RotateY, one, scalar, SinglyParameterisedMatrices::Ry );
TEST_ANY_CTRL_OPERATION( RotateZ, one, scalar, SinglyParameterisedMatrices::Rz );
