#define _USE_MATH_DEFINES

#include "quest/include/quest.h"

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


#define TEST_CATEGORY "[unit][operations]"



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

    qmatrix X = getPauliMatrix(1);
    qmatrix Y = getPauliMatrix(2);
    qmatrix Z = getPauliMatrix(3);

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

namespace ParameterisedMatrices {

    auto Rx = [](qreal p) { return getExponentialOfPauliMatrix(p, FixedMatrices::X); };
    auto Ry = [](qreal p) { return getExponentialOfPauliMatrix(p, FixedMatrices::Y); };
    auto Rz = [](qreal p) { return getExponentialOfPauliMatrix(p, FixedMatrices::Z); };

    auto PS  = [](qreal p) { return qmatrix{{1, 0}, {0, getExpI(p)}}; };
    auto PS2 = [](qreal p) { return getControlledMatrix(PS(p), 1); };
}

namespace VariableSizeMatrices {

    auto X  = [](int n) { return getKroneckerProduct(FixedMatrices::X, n); };
    auto PF = [](int n) { return getControlledMatrix(FixedMatrices::Z, n - 1); };
}

namespace VariableSizeParameterisedMatrices {

    auto Z = [](qreal p, int n) { 
        qmatrix m = getKroneckerProduct(FixedMatrices::Z, n);
        return getExponentialOfPauliMatrix(p, m);
    };

    auto PS = [](qreal p, int n) {
        qmatrix m = ParameterisedMatrices::PS(p);
        return getControlledMatrix(m, n - 1);
    };
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

void testQuregIsCorrectOnAllDeployments(quregCache& quregs, auto& reference, auto& function) {

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
 * - applyHadamard:         none
 * - applyRotateX:          scalar
 * - applyRotateAroundAxis: axisrots
 * - applyDiagMatr1:        diagmatr
 * - applyDiagMatrPower:    diagpower
 * - applyCompMatr:         compmatr
 * - applyPauliStr:         paulistr
 * - applyPauliGadgt:       pauligad
*/

enum ArgsFlag { none, scalar, axisrots, diagmatr, diagpower, compmatr, paulistr, pauligad };


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
        targsFlag == zero ||
        targsFlag == one  ||
        targsFlag == two  ||
        targsFlag == any);
}

void assertNumQubitsFlagsValid(
    NumQubitsFlag ctrlsFlag, NumQubitsFlag targsFlag, 
    vector<int> ctrls, vector<int> states, vector<int> targs
) {
    assertNumQubitsFlagsAreValid(ctrlsFlag, targsFlag);

    // we explicitly permit targsFlag=zero while
    // targs.size() != 0, which occurs when targs
    // are supplied to the API through an alternate
    // argument (e.g. a PauliStr)

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
        if constexpr (Targs == zero) operation(qureg, args...);
        if constexpr (Targs == one)  operation(qureg, targs[0], args...);
        if constexpr (Targs == two)  operation(qureg, targs[0], targs[1], args...);
        if constexpr (Targs == any)  operation(qureg, targs.data(), numTargs, args...);
    }
    if constexpr (Ctrls == one) {
        if constexpr (Targs == zero) operation(qureg, ctrls[0], args...);
        if constexpr (Targs == one)  operation(qureg, ctrls[0], targs[0], args...);
        if constexpr (Targs == two)  operation(qureg, ctrls[0], targs[0], targs[1], args...);
        if constexpr (Targs == any)  operation(qureg, ctrls[0], targs.data(), numTargs, args...);
    }
    if constexpr (Ctrls == any) {
        if constexpr (Targs == zero) operation(qureg, ctrls.data(), numCtrls, args...);
        if constexpr (Targs == one)  operation(qureg, ctrls.data(), numCtrls, targs[0], args...);
        if constexpr (Targs == two)  operation(qureg, ctrls.data(), numCtrls, targs[0], targs[1], args...);
        if constexpr (Targs == any)  operation(qureg, ctrls.data(), numCtrls, targs.data(), numTargs, args...);
    }
    if constexpr (Ctrls == anystates) {
        if constexpr (Targs == zero) operation(qureg, ctrls.data(), states.data(), numCtrls, args...);
        if constexpr (Targs == one)  operation(qureg, ctrls.data(), states.data(), numCtrls, targs[0], args...);
        if constexpr (Targs == two)  operation(qureg, ctrls.data(), states.data(), numCtrls, targs[0], targs[1], args...);
        if constexpr (Targs == any)  operation(qureg, ctrls.data(), states.data(), numCtrls, targs.data(), numTargs, args...);
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
        CompMatr cm = createCompMatr(numTargs); // must be freed
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
        DiagMatr dm = createDiagMatr(numTargs); // must be freed
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
auto getRandomRemainingArgs(vector<int> targs) {

    if constexpr (Args == none)
        return tuple{ };

    if constexpr (Args == scalar) {
        qreal angle = getRandomReal(-2*M_PI, 2*M_PI);
        return tuple{ angle };
    }

    if constexpr (Args == axisrots) {
        qreal angle = getRandomReal(-2*M_PI, 2*M_PI);
        qreal x = getRandomReal(-1, 1);
        qreal y = getRandomReal(-1, 1);
        qreal z = getRandomReal(-1, 1);
        return tuple{ angle, x, y, z };
    }

    if constexpr (Args == compmatr || Args == diagmatr) {
        auto matrix = getRandomApiMatrix<Targs,Args>(targs.size()); // allocates heap mem
        return tuple{ matrix };
    }

    if constexpr (Args == diagpower) {
        auto matrix = getRandomApiMatrix<Targs,Args>(targs.size()); // allocates heap mem
        qcomp exponent = getRandomComplex();
        return tuple{ matrix, exponent };
    }

    if constexpr (Args == paulistr) {
        PauliStr str = getRandomPauliStr(targs);
        return tuple{ str };
    }

    if constexpr (Args == pauligad) {
        PauliStr str = getRandomPauliStr(targs);
        qreal angle = getRandomReal(-2*M_PI, 2*M_PI);
        return tuple{ str, angle };
    }
}


template <NumQubitsFlag Targs, ArgsFlag Args>
void freeRemainingArgs(auto args) {

    if constexpr (Targs == any && Args == compmatr)
        destroyCompMatr(std::get<0>(args));

    if constexpr (Targs == any && Args == diagmatr)
        destroyDiagMatr(std::get<0>(args));

    if constexpr (Targs == any && Args == diagpower)
        destroyDiagMatr(std::get<0>(args));
}


/*
 * unpack the given reference operator matrix (a qmatrix) 
 * for an API operation, which is passed to testOperation(),
 * and which will be effected upon the reference state (a 
 * qvector or qmatrix). The type/form of matrixRefGen depends 
 * on the type of API operation, indicated by template parameter.
 */

template <NumQubitsFlag Targs, ArgsFlag Args>
qmatrix getReferenceMatrix(auto matrixRefGen, vector<int> targs, auto additionalArgs) {

    if constexpr (Args == none && Targs != any)
        return matrixRefGen;

    if constexpr (Args == none && Targs == any)
        return matrixRefGen(targs.size());

    if constexpr (Args == scalar && Targs != any) {
        qreal angle = std::get<0>(additionalArgs);
        return matrixRefGen(angle);
    }

    if constexpr (Args == scalar && Targs == any) {
        qreal angle = std::get<0>(additionalArgs);
        return matrixRefGen(angle, targs.size());
    }

    if constexpr (Args == axisrots) {
        qreal angle = std::get<0>(additionalArgs);
        qreal x = std::get<1>(additionalArgs);
        qreal y = std::get<2>(additionalArgs);
        qreal z = std::get<3>(additionalArgs);
        return getExponentialOfNormalisedPauliVector(angle, x, y, z);
    }

    if constexpr (Args == compmatr || Args == diagmatr) {
        auto apiMatrix = std::get<0>(additionalArgs);
        return getMatrix(apiMatrix);
    }

    if constexpr (Args == diagpower) {
        auto apiMatrix = std::get<0>(additionalArgs);
        qmatrix diag = getMatrix(apiMatrix);
        qcomp power = std::get<1>(additionalArgs);
        return getPowerOfDiagonalMatrix(diag, power);
    }

    if constexpr (Args == paulistr) {
        PauliStr str = std::get<0>(additionalArgs);
        return getMatrix(str, targs);
    }

    if constexpr (Args == pauligad) {
        PauliStr str = std::get<0>(additionalArgs);
        qreal angle = std::get<1>(additionalArgs);
        qmatrix matr = getMatrix(str, targs);
        return getExponentialOfPauliMatrix(angle, matr);
    }
}


/*
 * display all/only relevant inputs given to an 
 * API operation when its subsequent test fails.
 * This is like a customisation of CAPTURE, although
 * we must use UNSCOPED_INFO (and ergo re-implement
 * some printers) because our branching makes scopes
 * which end CAPTURE's lifetime.
 */


// TODO: surely this should live somewhere else,
// and/or re-use printer_ functions as much as possible

std::string toString(vector<int> list) {

    std::string out = "{ ";
    for (int& e : list)
        out += std::to_string(e) + " ";
    out += "}";
    return out;
}

extern int paulis_getPauliAt(PauliStr str, int ind);

std::string toString(PauliStr str, vector<int> targs) {

    std::string labels = "IXYZ";

    // ugly but adequate - like me (call me)
    std::string out = "";
    for (int i=targs.size()-1; i>=0; i--)
        out += labels[paulis_getPauliAt(str, targs[i])];
    
    return out;
}

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void CAPTURE_RELEVANT( vector<int> ctrls, vector<int> states, vector<int> targs, auto& args ) {

    // display ctrls
    if (Ctrls == one)
        UNSCOPED_INFO( "control := " << ctrls[0] );
    if (Ctrls == any || Ctrls == anystates )
        UNSCOPED_INFO( "controls := " << toString(ctrls) );

    // display states
    if (Ctrls == anystates)
        UNSCOPED_INFO( "states := " << toString(states) );

    // display targs
    if (Targs == one)
        UNSCOPED_INFO( "target := " << targs[0] );

    if (Targs == two) {
        UNSCOPED_INFO( "target A := " << targs[0] );
        UNSCOPED_INFO( "target B := " << targs[1] );
    }
    if (Targs == any)
        UNSCOPED_INFO( "targets := " << toString(targs) );

    // display rotation angle
    if constexpr (Args == scalar)
        UNSCOPED_INFO( "angle := " << std::get<0>(args) );

    // display rotation angle and axis
    if constexpr (Args == axisrots) {
        UNSCOPED_INFO( "angle := " << std::get<0>(args) );
        UNSCOPED_INFO( "x := " << std::get<1>(args) );
        UNSCOPED_INFO( "y := " << std::get<2>(args) );
        UNSCOPED_INFO( "z := " << std::get<3>(args) );
    }

    // note but don't display API matrices
    if constexpr (Args == compmatr || Args == diagmatr || Args == diagpower)
        UNSCOPED_INFO( "matrix := (omitted)" );

    // display exponent of diagonal matrix
    if constexpr (Args == diagpower) {
        qcomp p = std::get<1>(args);
        UNSCOPED_INFO( "exponent := " << real(p) << " + (" << imag(p) << ")i" );
    }

    // display PauliStr
    if constexpr (Args == paulistr || Args == pauligad)
        UNSCOPED_INFO( "paulis := " << toString(std::get<0>(args), targs) );

    // display PauliStr angle
    if constexpr (Args == pauligad)
        UNSCOPED_INFO( "angle := " << std::get<1>(args) );
}


/*
 * perform a unit test on an API operation. The
 * template parameters are compile-time clues
 * about what inputs to prepare and pass to the
 * operation, and how its reference matrix (arg
 * matrixRefGen) is formatted.
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperation(auto operation, auto matrixRefGen, bool multiplyOnly) {

    assertNumQubitsFlagsAreValid(Ctrls, Targs);

    // use existing cached Quregs
    auto statevecQuregs = getCachedStatevecs();
    auto densmatrQuregs = getCachedDensmatrs();

    SECTION( LABEL_CORRECTNESS ) {

        qvector statevecRef = getZeroVector(getPow2(NUM_UNIT_QUREG_QUBITS));
        qmatrix densmatrRef = getZeroMatrix(getPow2(NUM_UNIT_QUREG_QUBITS));

        // try all possible number of ctrls and targs
        int numTargs = GENERATE_NUM_TARGS<Ctrls,Targs,Args>(NUM_UNIT_QUREG_QUBITS);
        int numCtrls = GENERATE_NUM_CTRLS<Ctrls>(NUM_UNIT_QUREG_QUBITS - numTargs);
        
        // try all possible ctrls and targs
        auto listpair = GENERATE_COPY( disjointsublists(range(0,NUM_UNIT_QUREG_QUBITS), numCtrls, numTargs) );
        vector<int> ctrls = std::get<0>(listpair);
        vector<int> targs = std::get<1>(listpair);

        // randomise control states (if operation accepts them)
        vector<int> states = getRandomInts(0, 1+1, numCtrls * (Ctrls == anystates));

        // randomise remaining operation parameters
        auto primaryArgs = tuple{ ctrls, states, targs };
        auto furtherArgs = getRandomRemainingArgs<Targs,Args>(targs); // may allocate heap memory

        // obtain the reference matrix for this operation 
        qmatrix matrixRef = getReferenceMatrix<Targs,Args>(matrixRefGen, targs, furtherArgs);

        // PauliStr arg replaces target qubit list in API operations
        constexpr NumQubitsFlag RevTargs = (Args==paulistr||Args==pauligad)? zero : Targs;

        // many-target matrices are often non-unitary, so disable validation
        // (this is only really necessary for >= 5 targets on single and quad
        //  precision, but we disable more generally to avoid CI edge-cases)
        (Args == compmatr && numTargs >= 2)?
            setValidationEpsilon(0):
            setValidationEpsilonToDefault();

        // prepare test function which will receive both statevectors and density matrices
        auto testFunc = [&](Qureg qureg, auto& stateRef) -> void { 
            
            // invoke API operation, passing all args (unpacking variadic)
            auto apiFunc = [](auto&&... args) { return invokeApiOperation<Ctrls,RevTargs>(args...); };
            auto allArgs = std::tuple_cat(tuple{operation, qureg}, primaryArgs, furtherArgs);
            std::apply(apiFunc, allArgs);

            // update reference state
            (multiplyOnly)?
                multiplyReferenceOperator(stateRef, ctrls, states, targs, matrixRef):
                applyReferenceOperator(stateRef, ctrls, states, targs, matrixRef);
        };

        // report operation's input parameters if any subsequent test fails
        CAPTURE_RELEVANT<Ctrls,Targs,Args>( ctrls, states, targs, furtherArgs );

        // test API operation on all available deployment combinations (e.g. OMP, MPI, MPI+GPU, etc)
        SECTION( LABEL_STATEVEC ) { testQuregIsCorrectOnAllDeployments(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { testQuregIsCorrectOnAllDeployments(densmatrQuregs, densmatrRef, testFunc); }

        // free any heap-alloated API matrices
        freeRemainingArgs<Targs,Args>(furtherArgs);
    }

    // TODO:
    // input validation!
}

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperation(auto operation, auto matrixRefGen) {

    // default multiplyOnly=false
    testOperation<Ctrls,Targs,Args>(operation, matrixRefGen, false);
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
    TEST_CASE( "apply" #namesuffix,                     TEST_CATEGORY ) { testOperation<zero,     numtargs,argtype>( apply ## namesuffix,                     matrixgen); } \
    TEST_CASE( "applyControlled" #namesuffix,           TEST_CATEGORY ) { testOperation<one,      numtargs,argtype>( applyControlled ## namesuffix,           matrixgen); } \
    TEST_CASE( "applyMultiControlled" #namesuffix,      TEST_CATEGORY ) { testOperation<any,      numtargs,argtype>( applyMultiControlled ## namesuffix,      matrixgen); } \
    TEST_CASE( "applyMultiStateControlled" #namesuffix, TEST_CATEGORY ) { testOperation<anystates,numtargs,argtype>( applyMultiStateControlled ## namesuffix, matrixgen); } 



TEST_ANY_CTRL_OPERATION( PauliStr,    any, paulistr, nullptr );
TEST_ANY_CTRL_OPERATION( PauliGadget, any, pauligad, nullptr );

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

TEST_ANY_CTRL_OPERATION( RotateX, one, scalar, ParameterisedMatrices::Rx );
TEST_ANY_CTRL_OPERATION( RotateY, one, scalar, ParameterisedMatrices::Ry );
TEST_ANY_CTRL_OPERATION( RotateZ, one, scalar, ParameterisedMatrices::Rz );

TEST_ANY_CTRL_OPERATION( RotateAroundAxis, one, axisrots, nullptr );

TEST_ANY_CTRL_OPERATION( MultiQubitNot, any, none, VariableSizeMatrices::X );

TEST_ANY_CTRL_OPERATION( PhaseGadget, any, scalar, VariableSizeParameterisedMatrices::Z );

TEST_CASE( "multiplyPauliStr",     TEST_CATEGORY ) { testOperation<zero,any,paulistr>(multiplyPauliStr,    nullptr, true); }
TEST_CASE( "multiplyPauliGadget",  TEST_CATEGORY ) { testOperation<zero,any,pauligad>(multiplyPauliGadget, nullptr, true); }

TEST_CASE( "multiplyCompMatr1", TEST_CATEGORY ) { testOperation<zero,one,compmatr>(multiplyCompMatr1, nullptr, true); }
TEST_CASE( "multiplyCompMatr2", TEST_CATEGORY ) { testOperation<zero,two,compmatr>(multiplyCompMatr2, nullptr, true); }
TEST_CASE( "multiplyCompMatr",  TEST_CATEGORY ) { testOperation<zero,any,compmatr>(multiplyCompMatr,  nullptr, true); }

TEST_CASE( "multiplyDiagMatr1", TEST_CATEGORY ) { testOperation<zero,one,diagmatr>(multiplyDiagMatr1, nullptr, true); }
TEST_CASE( "multiplyDiagMatr2", TEST_CATEGORY ) { testOperation<zero,two,diagmatr>(multiplyDiagMatr2, nullptr, true); }
TEST_CASE( "multiplyDiagMatr",  TEST_CATEGORY ) { testOperation<zero,any,diagmatr>(multiplyDiagMatr,  nullptr, true); }

TEST_CASE( "multiplyDiagMatrPower",  TEST_CATEGORY ) { testOperation<zero,any,diagpower>(multiplyDiagMatrPower, nullptr, true); }

TEST_CASE( "multiplyMultiQubitNot", TEST_CATEGORY ) { testOperation<zero,any,none>(multiplyMultiQubitNot, VariableSizeMatrices::X, true); }
TEST_CASE( "multiplyPhaseGadget",   TEST_CATEGORY ) { testOperation<zero,any,scalar>(multiplyPhaseGadget, VariableSizeParameterisedMatrices::Z, true); }

TEST_CASE( "applyPhaseFlip",           TEST_CATEGORY ) { testOperation<zero,one,none>(applyPhaseFlip,           VariableSizeMatrices::PF(1)); }
TEST_CASE( "applyTwoQubitPhaseFlip",   TEST_CATEGORY ) { testOperation<zero,two,none>(applyTwoQubitPhaseFlip,   VariableSizeMatrices::PF(2)); }
TEST_CASE( "applyMultiQubitPhaseFlip", TEST_CATEGORY ) { testOperation<zero,any,none>(applyMultiQubitPhaseFlip, VariableSizeMatrices::PF); }

TEST_CASE( "applyPhaseShift",           TEST_CATEGORY ) { testOperation<zero,one,scalar>(applyPhaseShift,           ParameterisedMatrices::PS); }
TEST_CASE( "applyTwoQubitPhaseShift",   TEST_CATEGORY ) { testOperation<zero,two,scalar>(applyTwoQubitPhaseShift,   ParameterisedMatrices::PS2); }
TEST_CASE( "applyMultiQubitPhaseShift", TEST_CATEGORY ) { testOperation<zero,any,scalar>(applyMultiQubitPhaseShift, VariableSizeParameterisedMatrices::PS); }



/*
 * TODO:
 * UNTESTED FUNCTIONS BELOW
 */

void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);
void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);
void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);
void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);

void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);

void applySuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop);


int applyQubitMeasurement(Qureg qureg, int target);

int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability);

qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome);

void applyQubitProjector(Qureg qureg, int target, int outcome);

qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);

qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);

qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);

void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);


void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets);

void applyFullQuantumFourierTransform(Qureg qureg);