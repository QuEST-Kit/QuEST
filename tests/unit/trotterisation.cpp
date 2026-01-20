/** @file
 * Unit tests of the trotterisation module.
 *
 * @author Tyson Jones
 * @author Vasco Ferreira (initial Pauli permutation tests)
 * @author Maurice Jamieson (real and imaginary time evolution tests)
 * @author Oliver Thomson Brown (real and imaginary time evolution tests)
 * 
 * @defgroup unittrotter Trotterisation
 * @ingroup unittests
 */

#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/random.hpp"

#include <vector>
#include <string>

using std::vector;
using std::string;
using namespace Catch::Matchers;

/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[trotterisation]"

void TEST_ON_CACHED_QUREGS(quregCache quregs, auto& refFunc, auto& regularFunc, auto& randFunc) {
    for (auto& [label, refQureg]: quregs) {

        DYNAMIC_SECTION( label ) {
            initDebugState(refQureg);

            Qureg regularQureg = createCloneQureg(refQureg);
            Qureg randQureg = createCloneQureg(refQureg);

            refFunc(refQureg);
            regularFunc(regularQureg);
            randFunc(randQureg);

            double regularDistance = calcDistance(regularQureg, refQureg);
            double randDistance = calcDistance(randQureg, refQureg);

            REQUIRE( randDistance < regularDistance );

            destroyQureg(regularQureg);
            destroyQureg(randQureg);
        }
    }
}

/*
 * Prepare a Hamiltonian H under which dynamical
 * evolution will be simulated via Trotterisation
 * of unitary-time evolution operator e^(-itH).
 * If the Hamiltonian was fixed/known in advance,
 * we could instead use createInlinePauliStrSum()
 *
 * (Adapted from dynamics.cpp, @author Tyson Jones)
 */
PauliStrSum createHeisenbergHamiltonian(int numQubits) {

    // we prepare a Heisenberg XYZ spin-ring Hamiltonian,
    // i.e. H = -1/2 sum( Jx XX + Jy YY + Jz ZZ + h Z )
    // upon all nearest neighbour qubits, with periodicity.
    // The coefficients must be real for H to be Hermitian
    // and ergo its time-evolution operator to be unitary,
    // although they must be represented with a qcomp type.
    vector<string> operators = {"XX", "YY", "ZZ", "Z"};
    vector<qcomp> coefficients = {.1, .2, .3, .4}; // Jx,Jy,Jz,h

    // we will populate the below vectors with 4*numQubits
    // elements which we could pre-allocate with .reserve,
    // but we might incur Donald Knuth's justified wrath.
    vector<PauliStr> allStrings;
    vector<qcomp> allCoeffs;
    
    // prepare all XX + YY + ZZ
    for (int p=0; p<3; p++) {
        for (int i=0; i<numQubits; i++) {

            // A_i, A_i+1
            vector<int> targs = {i, (i+1)%numQubits};
            PauliStr str = getPauliStr(operators[p], targs);

            allStrings.push_back(str);
            allCoeffs.push_back(coefficients[p]);
        }
    }

    // prepare Z
    for (int i=0; i<numQubits; i++) {
        allStrings.push_back(getPauliStr(operators[3], {i}));
        allCoeffs.push_back(coefficients[3]);
    }

    // must be freed by caller
    return createPauliStrSum(allStrings, allCoeffs);
}

/*
 * Prepare the observable operator O under which the
 * evolved state (under H above) will be measured.
 * If this were one term (a single tensor product of
 * Pauli operators), we could return instead a PauliStr
 * but we here return an arbitrary weighted sum thereof.
 *
 * (Adapted from dynamics.cpp, @author Tyson Jones)
 */
PauliStrSum createAlternatingPauliObservable(int numQubits) {

    // we prepare a weighted sum of alternating Paulis
    // upon each qubit, i.e. 1 X0 + 2 Y1 + 3 Z2 + 1 X3 + ...
    // where the coefficients are real such that the
    // output observable is Hermitian.

    vector<PauliStr> strings(numQubits);
    vector<qcomp> coeffs(numQubits);

    for (int i=0; i<numQubits; i++) {
        strings[i] = getPauliStr({"XYZ"[i%3]}, {i});
        coeffs[i] = getQcomp(i%4 + 1, 0);
    }

    // must be freed by caller
    return createPauliStrSum(strings, coeffs);
}

/*
 * Constructs a PauliStrSum representing a 1D Hamiltonian of the form
 * H = - \mu \sum^{N}_{j} Z_{j} - J \sum_{<ij>}^{N} Z_{i}Z_{j}
 * where,
 * \mu = magField,
 * J = interactionStrength,
 * <ij> indicates nearest-neighbour interactions only,
 * and boundary conditions are periodic such that site N-1 interacts with site 0.
 *
 * The asymmetricBias term can be used to break the symmetry of the system
 * in order to 'choose' a preferred antiferromagnetic state, and ensure repeatable
 * predictable outcomes.
 * It adds a term of the form:
 * -BZ_{0}
 */
PauliStrSum createIsingHamiltonian(int numQubits, qreal magField, 
                                   qreal interactionStrength, qreal asymmetricBias) {
    const int NTERMS = 2 * numQubits + 1;
    
    vector<qcomp> coeffs;
    vector<PauliStr> pauli_terms;
    coeffs.reserve(NTERMS);
    pauli_terms.reserve(NTERMS);
    
    for (int i = 0; i < numQubits; ++i) {
        pauli_terms.push_back(getPauliStr("Z", {i}));
        coeffs.push_back(getQcomp(-magField, 0));
        
        int next = (i + 1) % numQubits;
        pauli_terms.push_back(getPauliStr("ZZ", {i, next}));
        coeffs.push_back(getQcomp(-interactionStrength, 0));
    }
    
    pauli_terms.push_back(getPauliStr("Z", {0}));
    coeffs.push_back(getQcomp(-asymmetricBias, 0));
    
    return createPauliStrSum(pauli_terms, coeffs);
}


/* 
 * TESTS
 */


/**
 * @todo
 * Basic validation for randomisation, should be expanded and merged
 * once the Trotterisation function tests have been implemented.
 */
TEST_CASE( "randomisedTrotter", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = getNumCachedQubits();
        int numTerms = 25;
        int reps = 50;
        double time = 1.0;

        int refOrder = 4;
        int order = GENERATE_COPY(1, 2);

        GENERATE( range(0, 10) );
        PauliStrSum sum = createRandomPauliStrSum(numQubits, numTerms);

        auto refFunc = [&](Qureg qureg) { applyTrotterizedUnitaryTimeEvolution(qureg, sum, time, refOrder, reps, false); };
        auto regularFunc = [&](Qureg qureg) { applyTrotterizedUnitaryTimeEvolution(qureg, sum, time, order, reps, false); };
        auto randFunc = [&](Qureg qureg) { applyTrotterizedUnitaryTimeEvolution(qureg, sum, time, order, reps, true); };
        
        TEST_ON_CACHED_QUREGS(getCachedDensmatrs(), refFunc, regularFunc, randFunc);
        TEST_ON_CACHED_QUREGS(getCachedStatevecs(), refFunc, regularFunc, randFunc);

        destroyPauliStrSum(sum);

    }
}

/*
* Time evolution tests
* @todo Add Pauli permutation variants
*/ 
TEST_CASE( "applyTrotterizedUnitaryTimeEvolution", TEST_CATEGORY ) { 

    // BEWARE: this test creates a new Qureg below which will have
    // deployments chosen by the auto-deployer; it is ergo unpredictable
    // whether it will be multithreaded, GPU-accelerated or distributed.
    // This test is ergo checking only a single, unspecified deployment,
    // unlike other tests which check all deployments. This is tolerable
    // since (non-randomised) Trotterisation is merely invoking routines
    // (Pauli gadgets) already independently tested across deployments

    SECTION( LABEL_CORRECTNESS ) {

        int numQubits = 20;
        Qureg qureg = createQureg(numQubits);
        initPlusState(qureg);
        bool permutePaulis = false;
        
        PauliStrSum hamil = createHeisenbergHamiltonian(numQubits);
        PauliStrSum observ = createAlternatingPauliObservable(numQubits);
        
        qreal dt = 0.1;
        int order = 4;
        int reps = 5;
        int steps = 10;
        
        // nudge the epsilon used by internal validation functions up a bit
        // as the time evolution operation plays badly with single precision
        // Defaults for validation epsilon are:
        //  - 1E-5 at single precision
        //  - 1E-12 at double precision
        //  - 1E-13 at quad precision
        qreal initialValidationEps = getValidationEpsilon();
        setValidationEpsilon(2 * initialValidationEps);

        /*
        * Tolerance for floating-point comparisons
        * Note that the underlying numerics are sensitive to the float
        * precision AND to the number of threads. As such we set quite 
        * large epsilon values to account for the worst-case scenario which 
        * is single precision, single thread. The baseline for these results
        * is double precision, multiple threads.
        *
        * Values (assuming default initialValidationEps) are:
        * Single precision:
        *   obsEps = 0.03
        *   normEps = 0.001
        *
        * Double precision:
        *   obsEps = 3E-9
        *   normEps = 1E-10
        *
        * Quad precision:
        *   obsEps = 3E-10
        *   normEps = 1E-11
        */
        qreal obsEps = 3E3 * initialValidationEps;
        qreal normEps = 100 * initialValidationEps;
       
        vector<qreal> refObservables = {
            19.26827777028073,
            20.34277275871839,
            21.21120737889526,
            21.86585902741717,
            22.30371711358924,
            22.52644660547882,
            22.54015748825067,
            22.35499202583118,
            21.9845541501027,
            21.44521638719462
        };
        
        for (int i = 0; i < steps; i++) {
            applyTrotterizedUnitaryTimeEvolution(qureg, hamil, dt, order, reps, permutePaulis);
            qreal expec = calcExpecPauliStrSum(qureg, observ);
            
            REQUIRE_THAT( expec, WithinAbs(refObservables[i], obsEps) );
        }
        
        // Verify state remains normalized
        REQUIRE_THAT( calcTotalProb(qureg), WithinAbs(1.0, normEps) );

        // Restore validation epsilon
        setValidationEpsilon(initialValidationEps);

        destroyQureg(qureg);
        destroyPauliStrSum(hamil);
        destroyPauliStrSum(observ);
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedStatevec();
        PauliStrSum hamil = createHeisenbergHamiltonian(qureg.numQubits);
        bool permutePaulis = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH( 
                applyTrotterizedUnitaryTimeEvolution(badQureg, hamil, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badHamil = hamil;
            badHamil.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, badHamil, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "hamiltonian not hermitian" ) {

            vector<PauliStr> strings;
            vector<qcomp> coeffs;
            strings.push_back(getPauliStr("X", {0}));
            coeffs.push_back(getQcomp(1.0, 1.0));  
            PauliStrSum nonHermitian = createPauliStrSum(strings, coeffs);

            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, nonHermitian, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeHamil = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, largeHamil, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeHamil);
        }

        SECTION( "invalid trotter order (zero)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, 0, 5, permutePaulis),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter order (negative)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, -2, 5, permutePaulis),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter order (odd, not 1)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, 3, 5, permutePaulis),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps (zero)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, 4, 0, permutePaulis),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "invalid trotter reps (negative)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, 4, -3, permutePaulis),
                ContainsSubstring("repetitions")
            );
        }

        destroyPauliStrSum(hamil);
    }
}


TEST_CASE( "applyTrotterizedImaginaryTimeEvolution", TEST_CATEGORY ) {

    // BEWARE: this test creates a new Qureg below which will have
    // deployments chosen by the auto-deployer; it is ergo unpredictable
    // whether it will be multithreaded, GPU-accelerated or distributed.
    // This test is ergo checking only a single, unspecified deployment,
    // unlike other tests which check all deployments. This is tolerable
    // since (non-randomised) Trotterisation is merely invoking routines
    // (Pauli gadgets) already independently tested across deployments

    SECTION( LABEL_CORRECTNESS ) {
           
        int numQubits = 16;
        qreal tau = 0.1;
        int order = 6;
        int reps = 5;
        int steps = 10;
        bool permutePaulis = false;
        
        // Tolerance for ground state amplitude
        qreal eps = 1E-2;
        
        // Ground state: all qubits align down (driven by strong magnetic field)
        {
            Qureg qureg = createQureg(numQubits);
            initPlusState(qureg);
            
            PauliStrSum ising = createIsingHamiltonian(numQubits, 10.0, 0.0, 0.0);
            
            for (int i = 0; i < steps; ++i) {
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, tau, order, reps, permutePaulis);
                setQuregToRenormalized(qureg);
            }
            
            qcomp amp = getQuregAmp(qureg, 0);
            qreal amp_mag = amp.real() * amp.real() + amp.imag() * amp.imag();
            
            REQUIRE_THAT( amp_mag, WithinAbs(1.0, eps) );
            
            for (long long i = 1; i < (1LL << numQubits); i++) {
                qcomp other_amp = getQuregAmp(qureg, i);
                qreal other_mag = other_amp.real() * other_amp.real() + 
                                 other_amp.imag() * other_amp.imag();
                REQUIRE( other_mag < eps );
            }
            
            destroyQureg(qureg);
            destroyPauliStrSum(ising);
        }
        
        // Ground state: all qubits align up (driven by opposite magnetic field)
        {
            Qureg qureg = createQureg(numQubits);
            initPlusState(qureg);
            
            PauliStrSum ising = createIsingHamiltonian(numQubits, -10.0, 0.0, 0.0);
            
            for (int i = 0; i < steps; ++i) {
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, tau, order, reps, permutePaulis);
                setQuregToRenormalized(qureg);
            }
            
            long long last_state = (1LL << numQubits) - 1;
            qcomp amp = getQuregAmp(qureg, last_state);
            qreal amp_mag = amp.real() * amp.real() + amp.imag() * amp.imag();
            
            REQUIRE_THAT( amp_mag, WithinAbs(1.0, eps) );
            
            for (long long i = 0; i < (1LL << numQubits); i++) {
                if (i == last_state) continue;
                qcomp other_amp = getQuregAmp(qureg, i);
                qreal other_mag = other_amp.real() * other_amp.real() + 
                                 other_amp.imag() * other_amp.imag();
                REQUIRE( other_mag < eps );
            }
            
            destroyQureg(qureg);
            destroyPauliStrSum(ising);
        }
        
        // Ground state: all qubits align down (driven by ferromagnetic interactions and bias)
        {
            Qureg qureg = createQureg(numQubits);
            initPlusState(qureg);
            
            PauliStrSum ising = createIsingHamiltonian(numQubits, 0.0, 10.0, 10.0);
            
            for (int i = 0; i < steps; ++i) {
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, tau, order, reps, permutePaulis);
                setQuregToRenormalized(qureg);
            }
            
            qcomp amp = getQuregAmp(qureg, 0);
            qreal amp_mag = amp.real() * amp.real() + amp.imag() * amp.imag();
            
            REQUIRE_THAT( amp_mag, WithinAbs(1.0, eps) );
            
            for (long long i = 1; i < (1LL << numQubits); i++) {
                qcomp other_amp = getQuregAmp(qureg, i);
                qreal other_mag = other_amp.real() * other_amp.real() + 
                                 other_amp.imag() * other_amp.imag();
                REQUIRE( other_mag < eps );
            }
            
            destroyQureg(qureg);
            destroyPauliStrSum(ising);
        }
        
        // Ground state: alternating pattern (driven by antiferromagnetic interactions)
        {
            Qureg qureg = createQureg(numQubits);
            initPlusState(qureg);
            
            PauliStrSum ising = createIsingHamiltonian(numQubits, 0.0, -10.0, 10.0);
            
            for (int i = 0; i < steps; ++i) {
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, tau, order, reps, permutePaulis);
                setQuregToRenormalized(qureg);
            }
            
            unsigned long long idx = 0;
            for (int i = 0; i < numQubits / 2; ++i) {
                idx += (1ULL << (2*i + 1));
            }
            
            qcomp amp = getQuregAmp(qureg, idx);
            qreal amp_mag = amp.real() * amp.real() + amp.imag() * amp.imag();
            
            REQUIRE_THAT( amp_mag, WithinAbs(1.0, eps) );
            
            for (long long i = 0; i < (1LL << numQubits); i++) {
                if (i == idx) continue;
                qcomp other_amp = getQuregAmp(qureg, i);
                qreal other_mag = other_amp.real() * other_amp.real() + 
                                 other_amp.imag() * other_amp.imag();
                REQUIRE( other_mag < eps );
            }
            
            destroyQureg(qureg);
            destroyPauliStrSum(ising);
        }
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedStatevec();
        PauliStrSum ising = createIsingHamiltonian(qureg.numQubits, 1.0, 1.0, 0.0);
        bool permutePaulis = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(badQureg, ising, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badIsing = ising;
            badIsing.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, badIsing, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeIsing = createIsingHamiltonian(qureg.numQubits+1, 1.0, 1.0, 0.0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, largeIsing, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeIsing);
        }

        SECTION( "hamiltonian not hermitian" ) {

            vector<PauliStr> strings;
            vector<qcomp> coeffs;
            strings.push_back(getPauliStr("X", {0}));
            coeffs.push_back(getQcomp(1.0, 1.0));  
            PauliStrSum nonHermitian = createPauliStrSum(strings, coeffs);

            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, nonHermitian, 0.1, 4, 5, permutePaulis),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order (zero)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, 0.1, 0, 5, permutePaulis),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter order (negative)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, 0.1, -2, 5, permutePaulis),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter order (odd, not 1)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, 0.1, 3, 5, permutePaulis),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps (zero)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, 0.1, 4, 0, permutePaulis),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "invalid trotter reps (negative)" ) {

            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, ising, 0.1, 4, -3, permutePaulis),
                ContainsSubstring("repetitions")
            );
        }

        destroyPauliStrSum(ising);
    }
}


/**
 * @todo
 * UNTESTED FUNCTIONS (NOT YET VALIDATED BY REFERENCE TESTS)
 */

void applyTrotterizedNonUnitaryPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps, bool permutePaulis);

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedControlledPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps, bool permutePaulis);

void applyTrotterizedNoisyTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal* damps, PauliStr* jumps, int numJumps, qreal time, int order, int reps, bool permutePaulis);
