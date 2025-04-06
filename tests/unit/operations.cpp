/** @file
 * Unit tests of the operations module. Beware that because the 
 * operation functions have so much interface and test-semantic 
 * overlap (e.g. the logic of control qubits, of control-states,
 * of density matrix variants), this file has opted to make 
 * extensive use of generics and metaprogramming, to avoid the
 * code duplication of defining each function an independent test.
 * This may have been a mistake, and the file now full of spaghetti
 * comprised of advanced, misused C++ facilities. View at your own
 * peril!
 *
 * @author Tyson Jones
 * 
 * @defgroup unitops Operations
 * @ingroup unittests
 */

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
#include "tests/utils/measure.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"

#include <tuple>

using std::tuple;
using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[operations]"


/*
 * reference operator matrices used by testing
 */

namespace FixedMatrices {

    qmatrix H = {
        {1/std::sqrt(2),  1/std::sqrt(2)},
        {1/std::sqrt(2), -1/std::sqrt(2)}};

    qmatrix X = getPauliMatrix(1);
    qmatrix Y = getPauliMatrix(2);
    qmatrix Z = getPauliMatrix(3);

    qreal PI = 3.14159265358979323846;
    qmatrix T = {
        {1, 0},
        {0, std::exp(1_i * PI/4)}};

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

    auto PS  = [](qreal p) { return qmatrix{{1, 0}, {0, std::exp(p*1_i)}}; };
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

void TEST_ON_CACHED_QUREGS(quregCache quregs, auto& reference, auto& function) {

    for (auto& [label, qureg]: quregs) {

        DYNAMIC_SECTION( label ) {

            // no need to validate whether qureg successfully
            // enters the debug state here, because the below
            // serial setToDebugState() is guaranteed to succeed
            initDebugState(qureg);
            setToDebugState(reference);

            function(qureg, reference);
            REQUIRE_AGREE( qureg, reference );
        }
    }
}

void TEST_ON_CACHED_QUREG_AND_MATRIX(quregCache quregs, matrixCache matrices, auto apiFunc, auto refState, auto refMatr, auto refFunc) {

    for (auto& [labelA, qureg]: quregs) {
        for (auto& [labelB, matrix]: matrices) {

            // skip illegal (distributed matrix, local qureg) combo
            if (matrix.isDistributed && ! qureg.isDistributed)
                continue;

            DYNAMIC_SECTION( labelA + LABEL_DELIMITER + labelB ) {

                // set qureg and reference to debug
                initDebugState(qureg);
                setToDebugState(refState);

                // set API matrix to pre-initialised ref matrix
                setFullStateDiagMatr(matrix, 0, getDiagonals(refMatr));

                // API and reference functions should produce agreeing states
                apiFunc(qureg, matrix);
                refFunc(refState, refMatr);
                REQUIRE_AGREE( qureg, refState );
            }
        }
    }
}


/*
 * simply avoids boilerplate
 */

#define PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef ) \
    int numQubits = getNumCachedQubits(); \
    auto statevecQuregs = getCachedStatevecs(); \
    auto densmatrQuregs = getCachedDensmatrs(); \
    qvector statevecRef = getZeroVector(getPow2(numQubits)); \
    qmatrix densmatrRef = getZeroMatrix(getPow2(numQubits));


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
 * the given template parameters. Depending on the
 * elemsFlag, the elements can be initialised to
 * zero (0), the identity matrix (1), or randomly (2).
 * This is used for testing API functions which accept 
 * matrices, in a function-agnostic way.
 */

template <NumQubitsFlag Targs, ArgsFlag Args>
auto getRandomOrIdentityApiMatrix(int numTargs, int elemsFlag) {

    DEMAND(
        Args == diagmatr ||
        Args == diagpower ||
        Args == compmatr );
    DEMAND( 
        elemsFlag == 0 ||
        elemsFlag == 1 ||
        elemsFlag == 2 );

    qmatrix qm;
    if (elemsFlag == 0)
        qm = getZeroMatrix(getPow2(numTargs));
    if (elemsFlag == 1)
        qm = getIdentityMatrix(getPow2(numTargs));
    if (elemsFlag == 2)
        qm = (Args == compmatr)?
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

template <NumQubitsFlag Targs, ArgsFlag Args> auto getZeroApiMatrix    (int numTargs) { return getRandomOrIdentityApiMatrix<Targs,Args>(numTargs, 0); }
template <NumQubitsFlag Targs, ArgsFlag Args> auto getIdentityApiMatrix(int numTargs) { return getRandomOrIdentityApiMatrix<Targs,Args>(numTargs, 1); }
template <NumQubitsFlag Targs, ArgsFlag Args> auto getRandomApiMatrix  (int numTargs) { return getRandomOrIdentityApiMatrix<Targs,Args>(numTargs, 2); }


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
        qreal angle = getRandomPhase();
        return tuple{ angle };
    }

    if constexpr (Args == axisrots) {
        qreal angle = getRandomPhase();
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
        DiagMatr matrix = getRandomApiMatrix<Targs,Args>(targs.size()); // allocates heap mem
        qcomp exponent = qcomp(getRandomReal(-3, 3), 0); // real for unitarity
        return tuple{ matrix, exponent };
    }

    if constexpr (Args == paulistr) {
        PauliStr str = getRandomPauliStr(targs);
        return tuple{ str };
    }

    if constexpr (Args == pauligad) {
        PauliStr str = getRandomPauliStr(targs);
        qreal angle = getRandomPhase();
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
        qreal angle  = std::get<1>(additionalArgs);
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


// @todo surely this should live somewhere else,
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
        UNSCOPED_INFO( "exponent := " << std::real(p) << " + (" << std::imag(p) << ")i" );
    }

    // display PauliStr
    if constexpr (Args == paulistr || Args == pauligad)
        UNSCOPED_INFO( "paulis := " << toString(std::get<0>(args), targs) );

    // display PauliStr angle
    if constexpr (Args == pauligad)
        UNSCOPED_INFO( "angle := " << std::get<1>(args) );
}


/*
 * test the correctness of an API operation. The
 * template parameters are compile-time clues
 * about what inputs to prepare and pass to the
 * operation, and how its reference matrix (arg
 * matrixRefGen) is formatted.
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperationCorrectness(auto operation, auto matrixRefGen, bool multiplyOnly) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    // try all possible number of ctrls and targs
    int numTargs = GENERATE_NUM_TARGS<Ctrls,Targs,Args>(numQubits);
    int numCtrls = GENERATE_NUM_CTRLS<Ctrls>(numQubits - numTargs);
    
    // either try all possible ctrls and targs, or randomise them
    auto ctrlsAndTargs = GENERATE_CTRLS_AND_TARGS( numQubits, numCtrls, numTargs );
    vector<int> ctrls = std::get<0>(ctrlsAndTargs);
    vector<int> targs = std::get<1>(ctrlsAndTargs);

    // randomise control states (if operation accepts them)
    vector<int> states = getRandomOutcomes(numCtrls * (Ctrls == anystates));

    // randomise remaining operation parameters
    auto primaryArgs = tuple{ ctrls, states, targs };
    auto furtherArgs = getRandomRemainingArgs<Targs,Args>(targs); // may allocate heap memory

    // obtain the reference matrix for this operation 
    qmatrix matrixRef = getReferenceMatrix<Targs,Args>(matrixRefGen, targs, furtherArgs);

    // PauliStr arg replaces target qubit list in API operations
    constexpr NumQubitsFlag RevTargs = (Args==paulistr||Args==pauligad)? zero : Targs;

    // disabling unitary-check validation for compmatr, since it's hard to
    // generate numerically-precise random unitaries upon many qubits, or
    // upon few qubits are single-precision. So we disable completely until
    // we re-implement 'input validation' checks which force us to fix thresholds
    (Args == compmatr)?
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
    SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
    SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }

    // free any heap-alloated API matrices and restore epsilon
    freeRemainingArgs<Targs,Args>(furtherArgs);
    setValidationEpsilonToDefault();
}


/*
 * test the input valdiation of an API operation. 
 * The template parameters are compile-time clues
 * about what inputs are accepted by the operation
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
auto getFixedCtrlsStatesTargs(int numQubits) {

    // default each to empty
    vector<int> targs, ctrls, states;

    // assign when non-empty
    if constexpr (Targs == one) targs = {0};
    if constexpr (Targs == two) targs = {0,1};
    if constexpr (Targs == any) targs = {0,1,2};
    if constexpr (Ctrls == one)       ctrls = {3};
    if constexpr (Ctrls == any)       ctrls = {3,4};
    if constexpr (Ctrls == anystates) ctrls = {3,4};
    if constexpr (Ctrls == anystates) states = {0,0};

    DEMAND( numQubits >= targs.size() + ctrls.size() );

    return tuple{ ctrls, states, targs };
}

template <NumQubitsFlag Targs, ArgsFlag Args>
auto getFixedRemainingArgs(vector<int> targs) {

    // getPauliStr uses gives length-3 hardcoded string
    if constexpr (Args == paulistr || Args == pauligad)
        DEMAND( targs.size() == 3 );

    if constexpr (Args == none)      return tuple{ };
    if constexpr (Args == scalar)    return tuple{ 0 }; // angle
    if constexpr (Args == axisrots)  return tuple{ 0, 1,1,1 }; // (angle,x,y,z)
    if constexpr (Args == compmatr)  return tuple{ getIdentityApiMatrix<Targs,Args>(targs.size()) }; // id
    if constexpr (Args == diagmatr)  return tuple{ getIdentityApiMatrix<Targs,Args>(targs.size()) }; // id
    if constexpr (Args == diagpower) return tuple{ getIdentityApiMatrix<Targs,Args>(targs.size()), qcomp(1,0) }; // (id, exponent)
    if constexpr (Args == paulistr)  return tuple{ getPauliStr("XXX", targs) };    // XXX
    if constexpr (Args == pauligad)  return tuple{ getPauliStr("XXX", targs), 0 }; // (XXX, angle)
}

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperationValidation(auto operation, bool multiplyOnly) {

    // use any cached Qureg
    Qureg qureg = getCachedStatevecs().begin()->second;

    // in lieu of preparing random inputs like testOperationCorrectness()
    // above, we instead obtain simple, fixed, compatible inputs
    auto [ctrls,states,targs] = getFixedCtrlsStatesTargs<Ctrls,Targs,Args>(qureg.numQubits);
    auto furtherArgs = getFixedRemainingArgs<Targs,Args>(targs);

    // calling apiFunc() will pass the above args with their call-time values
    auto apiFunc = [&]() {
        constexpr NumQubitsFlag RevTargs = (Args==paulistr||Args==pauligad)? zero : Targs;
        auto func = [](auto&&... allArgs) { return invokeApiOperation<Ctrls,RevTargs>(allArgs...); };
        std::apply(func, std::tuple_cat(tuple{operation, qureg, ctrls, states, targs}, furtherArgs));
    };

    // convenience vars
    int numQubits = qureg.numQubits;
    int numTargs = (int) targs.size();
    int numCtrls = (int) ctrls.size();

    /// @todo
    /// below, we return from intendedly skipped SECTIONS which
    /// appears to work (does not corrupt test statistics, and
    /// does not attribute skipped tests to having passed the
    /// section) but is an undocumented Catch2 trick. Check safe!

    SECTION( "qureg uninitialised" ) {

        // spoof uninitialised value
        qureg.numQubits = -123;
        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("invalid Qureg") );
    }

    SECTION( "invalid target" ) {

        // not applicable (PauliStr already made from targs)
        if (Args == paulistr || Args == pauligad)
            return;

        // sabotage a target
        int ind = GENERATE_COPY( range(0,numTargs) );
        int val = GENERATE_COPY( -1, numQubits );
        targs[ind] = val;

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("Invalid target qubit") );
    }

    SECTION( "invalid non-Identity Pauli index" ) {

        if (Args != paulistr && Args != pauligad)
            return;

        PauliStr badStr = getPauliStr("X", {numQubits + 1});
        if constexpr (Args == paulistr) furtherArgs = tuple{ badStr };
        if constexpr (Args == pauligad) furtherArgs = tuple{ badStr, 1 };

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("highest-index non-identity Pauli operator") && ContainsSubstring("exceeds the maximum target") );
    }

    SECTION( "invalid control" ) {

        if (numCtrls == 0)
            return;

        // sabotage a ctrl
        int ind = GENERATE_COPY( range(0,numCtrls) );
        int val = GENERATE_COPY( -1, numQubits );
        ctrls[ind] = val;

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("Invalid control qubit") );
    }

    SECTION( "control and target collision" ) {

        if (numCtrls == 0)
            return;

        // sabotage a ctrl
        int targInd = GENERATE_COPY( range(0,numTargs) );
        int ctrlInd = GENERATE_COPY( range(0,numCtrls) );
        ctrls[ctrlInd] = targs[targInd];

        if (Args==paulistr||Args==pauligad)
            REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("control qubit overlaps a non-identity Pauli operator") );
        else
            REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("qubit appeared among both the control and target qubits") );
    }

    SECTION( "control states" ) {

        if (states.empty())
            return;

        int ind = GENERATE_COPY( range(0,numCtrls) );
        int val = GENERATE( -1, 2 );
        states[ind] = val;

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("invalid control-state") );
    }

    SECTION( "repetition in controls" ) {

        if (numCtrls < 2)
            return;

        int ind = GENERATE_COPY( range(1,numCtrls) );
        ctrls[ind] = ctrls[0];

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("control qubits contained duplicates") );
    }

    SECTION( "repetition in targets" ) {

        // not applicable to Pauli functions (getPauliStr would throw)
        if (Args==paulistr||Args==pauligad)
            return;

        if (numTargs < 2)
            return;

        int ind = GENERATE_COPY( range(1,numTargs) );
        targs[ind] = targs[0];

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("target qubits contained duplicates") );
    }

    SECTION( "number of targets" ) {

        // not applicable to Pauli functions (getPauliStr would run fine,
        // and the runtime error would be about the non-identity index)
        if (Args == paulistr || Args == pauligad)
            return;

        if (Targs != any)
            return;
        
        // too few (cannot test less than 0)
        targs = {};
        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("Must specify one or more targets") );

        // too many; exceeds Qureg
        targs = getRange(qureg.numQubits+1);
        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("number of target qubits") && ContainsSubstring("exceeds the number of qubits in the Qureg") );

        // note numTargs + numCtrls > numQubits is caught by
        // invalid index or overlapping (ctrls,targs) validation
    }

    SECTION( "mismatching matrix size" ) {

        // only relevant to variable-sized matrices
        if (Targs != any)
            return;
        if (!(Args == compmatr || Args == diagmatr || Args == diagpower))
            return;
        
        DEMAND( numTargs > 1 );
        DEMAND( numCtrls + numTargs < numQubits ); 

        targs.push_back(numQubits - 1);
        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("matrix has an inconsistent size") );

        targs.pop_back();
        targs.pop_back();
        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("matrix has an inconsistent size") );
    }

    SECTION( "matrix unitarity" ) {

        // only relevant to matrix functions...
        if (Args != compmatr && Args != diagmatr && Args != diagpower)
            return;

        // which enforce unitarity
        if (multiplyOnly)
            return;

        if constexpr (Args == compmatr || Args == diagmatr)
            furtherArgs = tuple{ getZeroApiMatrix<Targs,Args>(targs.size()) };
        if constexpr (Args == diagpower)
            furtherArgs = tuple{ getZeroApiMatrix<Targs,Args>(targs.size()), 1 };
    
        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("unitary") );
    }

    SECTION( "matrix uninitialised" ) {

        // sabotage matrix struct field
        if constexpr (Args == compmatr || Args == diagmatr || Args == diagpower)
            std::get<0>(furtherArgs).numQubits = -1;
        else
            return;

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("Invalid") );

        // must correct field so that subsequent destructor doesn't whine!
        if constexpr (Args == compmatr || Args == diagmatr || Args == diagpower)
            std::get<0>(furtherArgs).numQubits = numTargs;
    }

    SECTION( "matrix unsynced" ) {

        if (!getQuESTEnv().isGpuAccelerated)
            return;

        // only relevant to variable-size matrix functions
        if constexpr (Targs == any && (Args == compmatr || Args == diagmatr || Args == diagpower))
            *(std::get<0>(furtherArgs).wasGpuSynced) = 0;
        else
            return; // avoid empty test

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("sync") );
    }

    SECTION( "targeted amps fit in node" ) {

        // can only be validated when environment AND qureg
        // are distributed (over more than 1 node, of course)
        if (qureg.numNodes < 2)
            return;

        // can only be validated if forced ctrl qubits permit 
        // enough remaining targets
        int minNumCtrls = (Ctrls == one)? 1 : 0;
        int minNumTargs = numQubits - qureg.logNumNodes + 1;
        int maxNumTargs = numQubits - minNumCtrls;
        if (minNumTargs > maxNumTargs)
            return;

        // only relevant to >=2-targ dense matrices, and further
        // only testable with any-targ variable-size matrices, since
        // 4x4 matrices might always be permissable by Qureg distribution
        if constexpr (Args == compmatr && Targs == any) {

            // free existing matrix to avoid leak
            destroyCompMatr(std::get<0>(furtherArgs));

            // try all illegally sized matrices
            int numNewTargs = GENERATE_COPY( range(minNumTargs, maxNumTargs+1) );
            targs = getRange(numNewTargs);

            // ensure no overlap with ctrls; just get rid of them, EXCEPT when the API
            // function expects explicitly one ctrl which we must always supply
            ctrls = vector<int>(minNumCtrls, numQubits - 1); // {} or {last}
            states = {};

            // create the new illegaly-sized matrix, which will be destroyed at test-case end
            CompMatr matr = getIdentityApiMatrix<Targs,Args>(numNewTargs);
            furtherArgs = tuple{ matr };

            CAPTURE( minNumCtrls, numNewTargs, numQubits - minNumCtrls, ctrls, targs );
            REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("cannot simultaneously store") && ContainsSubstring("remote amplitudes") );

        } else {
            return; // avoid empty test
        }
    }

    SECTION( "non-unitary exponent" ) {

        // not relevant for functions which do not assert unitarity
        if (multiplyOnly)
            return;

        if constexpr (Args == diagpower)
            furtherArgs = tuple{ std::get<0>(furtherArgs), qcomp(1,1) };
        else
            return; // avoid empty test

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("exponent was not approximately real") );
    }

    SECTION( "diverging exponent" ) {

        // when being applied as a unitary, abs(elem)=1 so there's no
        // possibility of divergence (we'd merely trigger isUnitary)
        if (!multiplyOnly)
            return;

        if constexpr (Args == diagpower)
            furtherArgs = tuple{ getZeroApiMatrix<Targs,Args>(numTargs), qcomp(-1,0) };
        else
            return; // avoid empty test

        REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("divergences") );
    }

    SECTION( "zero axis rotation" ) {

        if constexpr (Args == axisrots) {
            furtherArgs = tuple{ 0, 0, 0, 0 }; // (angle,x,y,z)
        } else
            return; // avoid empty test

            REQUIRE_THROWS_WITH( apiFunc(), ContainsSubstring("zero vector") );
    }

    freeRemainingArgs<Targs,Args>(furtherArgs);
}


/*
 * fully test an API operation, on compatible
 * inputs as indicated by the template flags
 */

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperation(auto operation, auto matrixRefGen, bool multiplyOnly) {

    assertNumQubitsFlagsAreValid(Ctrls, Targs);

    SECTION( LABEL_CORRECTNESS ) { 
        testOperationCorrectness<Ctrls,Targs,Args>(operation, matrixRefGen, multiplyOnly); 
    }

    SECTION( LABEL_VALIDATION ) { 
        testOperationValidation<Ctrls,Targs,Args>(operation, multiplyOnly); 
    }
}

template <NumQubitsFlag Ctrls, NumQubitsFlag Targs, ArgsFlag Args>
void testOperation(auto operation, auto matrixRefGen) {

    bool multiplyOnly = false;
    testOperation<Ctrls,Targs,Args>(operation, matrixRefGen, multiplyOnly);
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

// when every API operation had NO overloads, it was sufficient to
// call the below simple macro. Alas, since C++ overloads have been
// added, the below passing of a function (e.g. applyControlledHadamard)
// is ambigious, not specifying whether it is the int* or vector<int>
// version. We now have to explicit cast the function to its specific
// C-compatible version. Alas, to your imminent horror, this is.. erm...

// #define TEST_ALL_CTRL_OPERATIONS( namesuffix, numtargs, argtype, matrixgen ) \
//     TEST_CASE( "apply" #namesuffix,                     TEST_CATEGORY ) { testOperation<zero,     numtargs,argtype>( apply ## namesuffix,                     matrixgen); } \
//     TEST_CASE( "applyControlled" #namesuffix,           TEST_CATEGORY ) { testOperation<one,      numtargs,argtype>( applyControlled ## namesuffix,           matrixgen); } \
//     TEST_CASE( "applyMultiControlled" #namesuffix,      TEST_CATEGORY ) { testOperation<any,      numtargs,argtype>( applyMultiControlled ## namesuffix,      matrixgen); } \
//     TEST_CASE( "applyMultiStateControlled" #namesuffix, TEST_CATEGORY ) { testOperation<anystates,numtargs,argtype>( applyMultiStateControlled ## namesuffix, matrixgen); } 


/*
 * define macros to pre-processor-time generate the API function
 * signatures to give to static_cast to disambiguate the function
 * from its C++ overloads, as per the above comment. The absolute
 * macro nightmare below results from not being able to propagate
 * templates to static_cast<...> which is not actually templated!
 * May God forgive me for the misdeeds commited here below.
 */

// produces the control qubit arguments of a function signature,
// with a trailing comma (since ctrls always preceed more arguments)
#define GET_FUNC_CTRL_SUB_SIG( numctrls ) GET_FUNC_CTRL_SUB_SIG_##numctrls
#define GET_FUNC_CTRL_SUB_SIG_zero
#define GET_FUNC_CTRL_SUB_SIG_one       int,
#define GET_FUNC_CTRL_SUB_SIG_any       int*,int,
#define GET_FUNC_CTRL_SUB_SIG_anystates int*,int*,int,

// produces the target qubit arguments of a function signature, complicated
// by paulistr and pauligad functions not ever passign explicit targ lists.
// trailing comma attached only when targs exist, and final args exist.
// beware, macros must never have spaces between the name and open-paranthesis!
#define GET_FUNC_TARG_SUB_SIG( numtargs, argtype )  GET_FUNC_TARG_SUB_SIG_##argtype( numtargs )
#define GET_FUNC_TARG_SUB_SIG_none(      numtargs ) GET_FUNC_TARG_SUB_SIG_##numtargs
#define GET_FUNC_TARG_SUB_SIG_scalar(    numtargs ) GET_FUNC_TARG_SUB_SIG_##numtargs ,
#define GET_FUNC_TARG_SUB_SIG_axisrots(  numtargs ) GET_FUNC_TARG_SUB_SIG_##numtargs ,
#define GET_FUNC_TARG_SUB_SIG_compmatr(  numtargs ) GET_FUNC_TARG_SUB_SIG_##numtargs ,
#define GET_FUNC_TARG_SUB_SIG_diagmatr(  numtargs ) GET_FUNC_TARG_SUB_SIG_##numtargs ,
#define GET_FUNC_TARG_SUB_SIG_diagpower( numtargs ) GET_FUNC_TARG_SUB_SIG_##numtargs ,
#define GET_FUNC_TARG_SUB_SIG_paulistr(  numtargs )
#define GET_FUNC_TARG_SUB_SIG_pauligad(  numtargs )
#define GET_FUNC_TARG_SUB_SIG_one  int
#define GET_FUNC_TARG_SUB_SIG_two  int,int
#define GET_FUNC_TARG_SUB_SIG_any  int*,int

// produces the final arguments of a function signature (no trailing comma).
#define GET_FUNC_ARGS_SUB_SIG( numtargs, argtype ) GET_FUNC_ARGS_SUB_SIG_##numtargs##_##argtype
#define GET_FUNC_ARGS_SUB_SIG_one_none
#define GET_FUNC_ARGS_SUB_SIG_two_none
#define GET_FUNC_ARGS_SUB_SIG_any_none
#define GET_FUNC_ARGS_SUB_SIG_one_scalar qreal
#define GET_FUNC_ARGS_SUB_SIG_two_scalar qreal
#define GET_FUNC_ARGS_SUB_SIG_any_scalar qreal
#define GET_FUNC_ARGS_SUB_SIG_one_compmatr CompMatr1
#define GET_FUNC_ARGS_SUB_SIG_two_compmatr CompMatr2
#define GET_FUNC_ARGS_SUB_SIG_any_compmatr CompMatr
#define GET_FUNC_ARGS_SUB_SIG_one_diagmatr DiagMatr1
#define GET_FUNC_ARGS_SUB_SIG_two_diagmatr DiagMatr2
#define GET_FUNC_ARGS_SUB_SIG_any_diagmatr DiagMatr
#define GET_FUNC_ARGS_SUB_SIG_one_diagpower DiagMatr1,qcomp
#define GET_FUNC_ARGS_SUB_SIG_two_diagpower DiagMatr2,qcomp
#define GET_FUNC_ARGS_SUB_SIG_any_diagpower DiagMatr, qcomp
#define GET_FUNC_ARGS_SUB_SIG_one_axisrots qreal,qreal,qreal,qreal
#define GET_FUNC_ARGS_SUB_SIG_any_paulistr PauliStr
#define GET_FUNC_ARGS_SUB_SIG_any_pauligad PauliStr,qreal

// produces the control-qubit-related prefix of a function name
#define GET_FUNC_NAME_PREFIX( numctrls ) GET_FUNC_NAME_PREFIX_##numctrls
#define GET_FUNC_NAME_PREFIX_zero      apply
#define GET_FUNC_NAME_PREFIX_one       applyControlled
#define GET_FUNC_NAME_PREFIX_any       applyMultiControlled
#define GET_FUNC_NAME_PREFIX_anystates applyMultiStateControlled

// produces a function name from the control qubits and the suffix, e.g. (any,T) -> applyMultiControlledT
#define GET_FUNC_NAME(numctrls, suffix) GET_FUNC_NAME_INNER(GET_FUNC_NAME_PREFIX(numctrls), suffix)
#define GET_FUNC_NAME_INNER(A, B)       GET_FUNC_NAME_INNER_INNER(A, B)
#define GET_FUNC_NAME_INNER_INNER(A, B) A##B

// converts the output of GET_FUNC_NAME() to a string, e.g. (any,T) -> "applyMultiControlledT"
#define GET_FUNC_NAME_STR(numctrls, suffix) GET_FUNC_NAME_STR_INNER( GET_FUNC_NAME(numctrls,suffix) )
#define GET_FUNC_NAME_STR_INNER(expr)       GET_FUNC_NAME_STR_INNER_INNER(expr)
#define GET_FUNC_NAME_STR_INNER_INNER(symbol) #symbol

// produces the signature of a function, e.g. (any,one,diagmatr) -> void(*)(Qureg, int*,int, int,DiagMatr1),
// which is the signature of applyMultiControlled(Qureg qureg, int* ctrls, int numCtrls, int targ, DiagMatr1 m);
// NOTE:
// THIS CURRENT EXCLUDES PAULISTR AND PAULIGAD argtype
#define GET_FUNC_SIG( numctrls, numtargs, argtype ) \
    void(*) (                                       \
        Qureg,                                      \
        GET_FUNC_CTRL_SUB_SIG( numctrls )           \
        GET_FUNC_TARG_SUB_SIG( numtargs, argtype )  \
        GET_FUNC_ARGS_SUB_SIG( numtargs, argtype )  \
    )

// produces a function name, casted to its explicit C-argument form, disambiguated from its C++ overload.
// e.g. (DiagMatrPower, zero, any, diagpower) -> static_cast<void(*)(Qureg,int*,int,DiagMatr,qcomp)>(applyDiagMatrPower)>
#define GET_CASTED_FUNC( namesuffix, numctrls, numtargs, argtype ) \
    static_cast< GET_FUNC_SIG(numctrls, numtargs, argtype) > (     \
        GET_FUNC_NAME(numctrls, namesuffix) )

// defines a Catch2 test-case for the implied function
#define TEST_CASE_OPERATION( namesuffix, numctrls, numtargs, argtype, matrixgen ) \
    TEST_CASE( GET_FUNC_NAME_STR(numctrls, namesuffix), TEST_CATEGORY ) {         \
        testOperation<numctrls, numtargs, argtype>(                               \
            GET_CASTED_FUNC(namesuffix, numctrls, numtargs, argtype),             \
            matrixgen);                                                           \
    }
 
// automate the testing of a function for all its controlled variants
#define TEST_ALL_CTRL_OPERATIONS( namesuffix, numtargs, argtype, matrixgen ) \
    TEST_CASE_OPERATION( namesuffix, zero,      numtargs, argtype, matrixgen ) \
    TEST_CASE_OPERATION( namesuffix, one,       numtargs, argtype, matrixgen ) \
    TEST_CASE_OPERATION( namesuffix, any,       numtargs, argtype, matrixgen ) \
    TEST_CASE_OPERATION( namesuffix, anystates, numtargs, argtype, matrixgen )



/** 
 * TESTS
 * 
 * @ingroup unitops
 * @{
 */


/*
 * controlled operations
 */

TEST_ALL_CTRL_OPERATIONS( PauliStr,    any, paulistr, nullptr );
TEST_ALL_CTRL_OPERATIONS( PauliGadget, any, pauligad, nullptr );
TEST_ALL_CTRL_OPERATIONS( CompMatr1, one, compmatr, nullptr );
TEST_ALL_CTRL_OPERATIONS( CompMatr2, two, compmatr, nullptr );
TEST_ALL_CTRL_OPERATIONS( CompMatr,  any, compmatr, nullptr );
TEST_ALL_CTRL_OPERATIONS( DiagMatr1, one, diagmatr, nullptr );
TEST_ALL_CTRL_OPERATIONS( DiagMatr2, two, diagmatr, nullptr );
TEST_ALL_CTRL_OPERATIONS( DiagMatr,  any, diagmatr, nullptr );
TEST_ALL_CTRL_OPERATIONS( DiagMatrPower, any, diagpower, nullptr );
TEST_ALL_CTRL_OPERATIONS( Hadamard, one, none, FixedMatrices::H );
TEST_ALL_CTRL_OPERATIONS( PauliX,   one, none, FixedMatrices::X );
TEST_ALL_CTRL_OPERATIONS( PauliY,   one, none, FixedMatrices::Y );
TEST_ALL_CTRL_OPERATIONS( PauliZ,   one, none, FixedMatrices::Z );
TEST_ALL_CTRL_OPERATIONS( T,        one, none, FixedMatrices::T );
TEST_ALL_CTRL_OPERATIONS( S,        one, none, FixedMatrices::S );
TEST_ALL_CTRL_OPERATIONS( Swap,     two, none, FixedMatrices::SWAP );
TEST_ALL_CTRL_OPERATIONS( SqrtSwap, two, none, FixedMatrices::sqrtSWAP );
TEST_ALL_CTRL_OPERATIONS( RotateX, one, scalar, ParameterisedMatrices::Rx );
TEST_ALL_CTRL_OPERATIONS( RotateY, one, scalar, ParameterisedMatrices::Ry );
TEST_ALL_CTRL_OPERATIONS( RotateZ, one, scalar, ParameterisedMatrices::Rz );
TEST_ALL_CTRL_OPERATIONS( RotateAroundAxis, one, axisrots, nullptr );
TEST_ALL_CTRL_OPERATIONS( MultiQubitNot, any, none, VariableSizeMatrices::X );
TEST_ALL_CTRL_OPERATIONS( PhaseGadget, any, scalar, VariableSizeParameterisedMatrices::Z );


/*
 * non-controlled operations with no C++ overloads
 */

TEST_CASE( "multiplyPauliStr",        TEST_CATEGORY ) { testOperation<zero,any,paulistr>(multiplyPauliStr,    nullptr, true); }
TEST_CASE( "multiplyPauliGadget",     TEST_CATEGORY ) { testOperation<zero,any,pauligad>(multiplyPauliGadget, nullptr, true); }
TEST_CASE( "multiplyCompMatr1",       TEST_CATEGORY ) { testOperation<zero,one,compmatr>(multiplyCompMatr1,   nullptr, true); }
TEST_CASE( "multiplyCompMatr2",       TEST_CATEGORY ) { testOperation<zero,two,compmatr>(multiplyCompMatr2,   nullptr, true); }
TEST_CASE( "multiplyDiagMatr1",       TEST_CATEGORY ) { testOperation<zero,one,diagmatr>(multiplyDiagMatr1,   nullptr, true); }
TEST_CASE( "multiplyDiagMatr2",       TEST_CATEGORY ) { testOperation<zero,two,diagmatr>(multiplyDiagMatr2,   nullptr, true); }
TEST_CASE( "applyPhaseFlip",          TEST_CATEGORY ) { testOperation<zero,one,none>  (applyPhaseFlip,          VariableSizeMatrices::PF(1)); }
TEST_CASE( "applyTwoQubitPhaseFlip",  TEST_CATEGORY ) { testOperation<zero,two,none>  (applyTwoQubitPhaseFlip,  VariableSizeMatrices::PF(2)); }
TEST_CASE( "applyPhaseShift",         TEST_CATEGORY ) { testOperation<zero,one,scalar>(applyPhaseShift,         ParameterisedMatrices::PS); }
TEST_CASE( "applyTwoQubitPhaseShift", TEST_CATEGORY ) { testOperation<zero,two,scalar>(applyTwoQubitPhaseShift, ParameterisedMatrices::PS2); }


/*
 * non-controlled operations which have a C++ overload
 * (because they accept qubit lists which become vector),
 * and so which require explicit casting to resolve the
 * compiler ambiguity
 */

TEST_CASE( "multiplyCompMatr",  TEST_CATEGORY ) { 
    auto func = static_cast<void(*)(Qureg, int*, int, CompMatr)>(multiplyCompMatr);
    testOperation<zero,any,compmatr>(func, nullptr, true); 
}

TEST_CASE( "multiplyDiagMatr",  TEST_CATEGORY ) {
    auto func = static_cast<void(*)(Qureg, int*, int, DiagMatr)>(multiplyDiagMatr);
    testOperation<zero,any,diagmatr>(func, nullptr, true);
}

TEST_CASE( "multiplyDiagMatrPower",  TEST_CATEGORY ) {
    auto func = static_cast<void(*)(Qureg, int*, int, DiagMatr, qcomp)>(multiplyDiagMatrPower);
    testOperation<zero,any,diagpower>(func, nullptr, true);
}

TEST_CASE( "multiplyMultiQubitNot",  TEST_CATEGORY ) {
    auto func = static_cast<void(*)(Qureg, int*, int)>(multiplyMultiQubitNot);
    testOperation<zero,any,none>(func, VariableSizeMatrices::X, true);
}

TEST_CASE( "multiplyPhaseGadget",  TEST_CATEGORY ) {
    auto func = static_cast<void(*)(Qureg, int*, int, qreal)>(multiplyPhaseGadget);
    testOperation<zero,any,scalar>(func, VariableSizeParameterisedMatrices::Z, true);
}

TEST_CASE( "applyMultiQubitPhaseFlip",  TEST_CATEGORY ) {
    auto func = static_cast<void(*)(Qureg, int*, int)>(applyMultiQubitPhaseFlip);
    testOperation<zero,any,none>(func, VariableSizeMatrices::PF);
}

TEST_CASE( "applyMultiQubitPhaseShift",  TEST_CATEGORY ) {
    auto func = static_cast<void(*)(Qureg, int*, int, qreal)>(applyMultiQubitPhaseShift);
    testOperation<zero,any,scalar>(func, VariableSizeParameterisedMatrices::PS);
}


/*
 * operations which need custom logic
 */

TEST_CASE( "applyQuantumFourierTransform", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targs = GENERATE_TARGS( numQubits, numTargs );

        CAPTURE( targs );

        SECTION( LABEL_STATEVEC ) { 

            auto testFunc = [&](Qureg qureg, qvector& ref) {
                applyQuantumFourierTransform(qureg, targs.data(), targs.size());
                ref = getDisceteFourierTransform(ref, targs);
            };

            TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc);
        }

        SECTION( LABEL_DENSMATR ) { 

            // prepare a random mixture
            auto states = getRandomOrthonormalStateVectors(numQubits, getRandomInt(1,10));
            auto probs = getRandomProbabilities(states.size());

            auto testFunc = [&](Qureg qureg, qmatrix& ref) {

                // overwrite the Qureg debug state set by caller to above mixture
                setQuregToReference(qureg, getMixture(states, probs));
                applyQuantumFourierTransform(qureg, targs.data(), targs.size());
                
                ref = getZeroMatrix(ref.size());
                for (size_t i=0; i<states.size(); i++) {
                    qvector vec = getDisceteFourierTransform(states[i], targs);
                    ref += probs[i] * getOuterProduct(vec, vec);
                }
            };

            TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc);
        }
    }

    /// @todo input validation
}


TEST_CASE( "applyFullQuantumFourierTransform", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );

        SECTION( LABEL_STATEVEC ) { 

            auto testFunc = [&](Qureg qureg, qvector& ref) {
                applyFullQuantumFourierTransform(qureg);
                ref = getDisceteFourierTransform(ref);
            };

            TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc);
        }

        SECTION( LABEL_DENSMATR ) { 

            // prepare a random mixture
            auto states = getRandomOrthonormalStateVectors(numQubits, getRandomInt(1,10));
            auto probs = getRandomProbabilities(states.size());

            auto testFunc = [&](Qureg qureg, qmatrix& ref) {

                // overwrite the Qureg debug state set by caller to above mixture
                setQuregToReference(qureg, getMixture(states, probs));
                applyFullQuantumFourierTransform(qureg);
                
                ref = getZeroMatrix(ref.size());
                for (size_t i=0; i<states.size(); i++) {
                    qvector vec = getDisceteFourierTransform(states[i]);
                    ref += probs[i] * getOuterProduct(vec, vec);
                }
            };

            TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc);
        }
    }

    /// @todo input validation
}


TEST_CASE( "applyQubitProjector", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        int target = GENERATE_COPY( range(0,numQubits) );
        int outcome = GENERATE( 0, 1 );

        qmatrix projector = getProjector(outcome);

        auto testFunc = [&](Qureg qureg, auto& ref) {
            applyQubitProjector(qureg, target, outcome);
            applyReferenceOperator(ref, {target}, projector);
        };

        CAPTURE( target, outcome );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "applyMultiQubitProjector", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_TARGS( numQubits, numTargs );
        auto outcomes = getRandomOutcomes(numTargs);

        qmatrix projector = getProjector(targets, outcomes, numQubits);

        auto testFunc = [&](Qureg qureg, auto& ref) {
            applyMultiQubitProjector(qureg, targets.data(), outcomes.data(), numTargs);
            applyReferenceOperator(ref, projector);
        };

        CAPTURE( targets, outcomes );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "applyForcedQubitMeasurement", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        int target = GENERATE_COPY( range(0,numQubits) );
        int outcome = GENERATE( 0, 1 );

        qmatrix projector = getProjector(outcome);

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // overwrite caller's setting of initDebugState, since
            // that precludes outcomes=|0><0| due to zero-probability
            setToRandomState(ref);
            setQuregToReference(qureg, ref);

            // compare the probabilities...
            qreal apiProb = applyForcedQubitMeasurement(qureg, target, outcome);
            qreal refProb = getReferenceProbability(ref, {target}, {outcome});
            REQUIRE_AGREE( apiProb, refProb );

            // and the post-projection states (caller calls subsequent REQUIRE_AGREE)
            applyReferenceOperator(ref, {target}, projector);
            ref /= (qureg.isDensityMatrix)?
                refProb : std::sqrt(refProb);
        };

        CAPTURE( target, outcome );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "applyForcedMultiQubitMeasurement", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_TARGS( numQubits, numTargs );
        auto outcomes = getRandomOutcomes(numTargs);

        qmatrix projector = getProjector(targets, outcomes, numQubits);

        // this test may randomly request a measurement outcome which
        // is illegally unlikely, triggering validation; we merely
        // disable such validation and hope divergences don't break the test!
        setValidationEpsilon(0);

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // overwrite caller's setting of initDebugState, since
            // that precludes outcomes=|0><0| due to zero-probability
            setToRandomState(ref);
            setQuregToReference(qureg, ref);

            // compare the probabilities...
            qreal apiProb = applyForcedMultiQubitMeasurement(qureg, targets.data(), outcomes.data(), numTargs);
            qreal refProb = getReferenceProbability(ref, targets, outcomes);
            REQUIRE_AGREE( apiProb, refProb );

            // and the post-measurement states (caller calls subsequent REQUIRE_AGREE)
            applyReferenceOperator(ref, projector);
            ref /= (qureg.isDensityMatrix)?
                refProb : std::sqrt(refProb);
        };

        CAPTURE( targets, outcomes );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }

        setValidationEpsilonToDefault();
    }

    /// @todo input validation
}


TEST_CASE( "applyMultiQubitMeasurement", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_TARGS( numQubits, numTargs );

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // overwrite caller's setting of initDebugState, since
            // sampling requires the outcome probs are normalised
            setToRandomState(ref);
            setQuregToReference(qureg, ref);
            
            // the output API state...
            qindex apiOut = applyMultiQubitMeasurement(qureg, targets.data(), numTargs);

            // informs the projector which determines the post-measurement reference
            auto apiOutBits = getBits(apiOut, numTargs);
            qmatrix projector = getProjector(targets, apiOutBits, numQubits);
            applyReferenceOperator(ref, projector);
            qreal refProb = getReferenceProbability(ref, targets, apiOutBits);
            ref /= (qureg.isDensityMatrix)?
                refProb : std::sqrt(refProb);
        };

        CAPTURE( targets );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "applyMultiQubitMeasurementAndGetProb", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        int numTargs = GENERATE_COPY( range(1,numQubits+1) );
        auto targets = GENERATE_TARGS( numQubits, numTargs );

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // overwrite caller's setting of initDebugState, since
            // sampling requires the outcome probs are normalised
            setToRandomState(ref);
            setQuregToReference(qureg, ref);

            // compare the measurement probability...
            qreal apiProb = -1;
            qindex apiOut = applyMultiQubitMeasurementAndGetProb(qureg, targets.data(), numTargs, &apiProb);
            auto apiOutBits = getBits(apiOut, numTargs);
            qreal refProb = getReferenceProbability(ref, targets, apiOutBits);
            REQUIRE_AGREE( apiProb, refProb );
            
            // and the post-measurement states (caller calls subsequent REQUIRE_AGREE)
            qmatrix projector = getProjector(targets, apiOutBits, numQubits);
            applyReferenceOperator(ref, projector);
            ref /= (qureg.isDensityMatrix)?
                refProb : std::sqrt(refProb);
        };

        CAPTURE( targets );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "applyQubitMeasurement", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        int target = GENERATE_COPY( range(0,numQubits) );

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // overwrite caller's setting of initDebugState, since
            // sampling requires the outcome probs are normalised
            setToRandomState(ref);
            setQuregToReference(qureg, ref);

            // the output API state...
            int apiOut = applyQubitMeasurement(qureg, target);

            // informs the projector which determines the post-measurement reference
            qmatrix projector = getProjector(apiOut);
            applyReferenceOperator(ref, {target}, projector);
            qreal refProb = getReferenceProbability(ref, {target}, {apiOut});
            ref /= (qureg.isDensityMatrix)?
                refProb : std::sqrt(refProb);
        };

        CAPTURE( target );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "applyQubitMeasurementAndGetProb", TEST_CATEGORY ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );
        int target = GENERATE_COPY( range(0,numQubits) );

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // overwrite caller's setting of initDebugState, since
            // sampling requires the outcome probs are normalised
            setToRandomState(ref);
            setQuregToReference(qureg, ref);

            // compare the measurement probability...
            qreal apiProb = -1;
            int apiOut = applyQubitMeasurementAndGetProb(qureg, target, &apiProb);
            qreal refProb = getReferenceProbability(ref, {target}, {apiOut});
            REQUIRE_AGREE( apiProb, refProb );
            
            // and the post-measurement states (caller calls subsequent REQUIRE_AGREE)
            qmatrix projector = getProjector(apiOut);
            applyReferenceOperator(ref, {target}, projector);
            ref /= (qureg.isDensityMatrix)?
                refProb : std::sqrt(refProb);            
        };

        CAPTURE( target );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


TEST_CASE( "multiplyFullStateDiagMatr", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    PREPARE_TEST( numQubits, cachedSV, cachedDM, refSV, refDM );

    auto cachedMatrs = getCachedFullStateDiagMatrs();

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalMatrix(getPow2(numQubits));
        auto apiFunc = multiplyFullStateDiagMatr;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector& state, qmatrix matr) { multiplyReferenceOperator(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedSV, cachedMatrs, apiFunc, refSV, refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix& state, qmatrix matr) { multiplyReferenceOperator(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedDM, cachedMatrs, apiFunc, refDM, refMatr, refFunc);
        }
    }

    /// @todo input validation
}


TEST_CASE( "multiplyFullStateDiagMatrPower", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    PREPARE_TEST( numQubits, cachedSV, cachedDM, refSV, refDM );

    auto cachedMatrs = getCachedFullStateDiagMatrs();

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalMatrix(getPow2(numQubits));
        qcomp exponent = getRandomComplex();

        auto apiFunc = [&](Qureg qureg, FullStateDiagMatr matr) { 
            return multiplyFullStateDiagMatrPower(qureg, matr, exponent);
        };

        CAPTURE( exponent );
        
        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector& state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, exponent);
                multiplyReferenceOperator(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedSV, cachedMatrs, apiFunc, refSV, refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix& state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, exponent);
                multiplyReferenceOperator(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedDM, cachedMatrs, apiFunc, refDM, refMatr, refFunc);
        }
    }

    /// @todo input validation
}


TEST_CASE( "applyFullStateDiagMatr", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    PREPARE_TEST( numQubits, cachedSV, cachedDM, refSV, refDM );

    auto cachedMatrs = getCachedFullStateDiagMatrs();

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalUnitary(numQubits);
        auto apiFunc = applyFullStateDiagMatr;

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector& state, qmatrix matr) { applyReferenceOperator(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedSV, cachedMatrs, apiFunc, refSV, refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix& state, qmatrix matr) { applyReferenceOperator(state, matr); };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedDM, cachedMatrs, apiFunc, refDM, refMatr, refFunc);
        }
    }

    /// @todo input validation
}


TEST_CASE( "applyFullStateDiagMatrPower", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    PREPARE_TEST( numQubits, cachedSV, cachedDM, refSV, refDM );

    auto cachedMatrs = getCachedFullStateDiagMatrs();

    SECTION( LABEL_CORRECTNESS ) {

        qmatrix refMatr = getRandomDiagonalUnitary(numQubits);

        // supplying a complex exponent requires disabling
        // numerical validation to relax unitarity
        bool testRealExp = GENERATE( true, false );
        qcomp exponent = (testRealExp)?
            qcomp(getRandomReal(-2, 2), 0):
            getRandomComplex();

        auto apiFunc = [&](Qureg qureg, FullStateDiagMatr matr) { 
            return applyFullStateDiagMatrPower(qureg, matr, exponent);
        };

        CAPTURE( exponent );

        GENERATE( range(0, TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS) );

        if (!testRealExp)
            setValidationEpsilon(0);

        SECTION( LABEL_STATEVEC ) {

            auto refFunc = [&] (qvector& state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, exponent);
                applyReferenceOperator(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedSV, cachedMatrs, apiFunc, refSV, refMatr, refFunc);
        }

        SECTION( LABEL_DENSMATR ) {

            auto refFunc = [&] (qmatrix& state, qmatrix matr) { 
                matr = getPowerOfDiagonalMatrix(matr, exponent);
                applyReferenceOperator(state, matr);
            };

            TEST_ON_CACHED_QUREG_AND_MATRIX( cachedDM, cachedMatrs, apiFunc, refDM, refMatr, refFunc);
        }

        setValidationEpsilonToDefault();
    }

    /// @todo input validation
}


TEST_CASE( "multiplyPauliStrSum", TEST_CATEGORY LABEL_MIXED_DEPLOY_TAG ) {

    PREPARE_TEST( numQubits, statevecQuregs, densmatrQuregs, statevecRef, densmatrRef );

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTerms = GENERATE_COPY( 1, 2, 10 );

        PauliStrSum sum = createRandomPauliStrSum(numQubits, numTerms);

        auto testFunc = [&](Qureg qureg, auto& ref) {

            // must use (and ergo make) an identically-deployed workspace
            Qureg workspace = createCloneQureg(qureg);
            multiplyPauliStrSum(qureg, sum, workspace);
            destroyQureg(workspace);

            ref = getMatrix(sum, numQubits) * ref;
        };

        CAPTURE( numTerms );
        SECTION( LABEL_STATEVEC ) { TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, testFunc); }
        SECTION( LABEL_DENSMATR ) { TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, testFunc); }
    }

    /// @todo input validation
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);
