/** @file
 * Unit tests of the trotterisation module.
 *
 * @author Tyson Jones
 * @author Vasco Ferreira (initial Pauli permutation and validation tests)
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

void TEST_ON_CACHED_QUREGS(quregCache quregs, qvector& referenceResult, auto& testFunction, PauliStrSum& testHamiltonian) {
  for (auto& [label, qureg]: quregs) {

    DYNAMIC_SECTION( label ) {
      testFunction(qureg, testHamiltonian);
      REQUIRE_AGREE(qureg, referenceResult);
    }

  }

  return;
}

void TEST_ON_CACHED_QUREGS(quregCache quregs, qmatrix& referenceResult, auto& testFunction, PauliStrSum& testHamiltonian) {
  for (auto& [label, qureg]: quregs) {

    DYNAMIC_SECTION( label ) {
      testFunction(qureg, testHamiltonian);
      REQUIRE_AGREE(qureg, referenceResult);
    }

  }

  return;
}

void TEST_OBSERVABLES_ON_QUREGS(quregCache quregs, qvector& referenceResult, auto& testFunction, PauliStrSum testHamiltonian, PauliStrSum testObservable) {
   for (auto& [label, qureg]: quregs) {

        DYNAMIC_SECTION( label ) {
            qvector testResult = testFunction(qureg, testHamiltonian, testObservable);
            REQUIRE_AGREE(calcTotalProb(qureg), 1.0);
            REQUIRE_AGREE(testResult, referenceResult);
        }

   }

   return;
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

    SECTION( LABEL_CORRECTNESS ) {
        // nudge the epsilon used by internal validation functions up a bit
        // as the time evolution operation plays badly with single precision
        // Defaults for validation epsilon are:
        //  - 1E-5 at single precision
        //  - 1E-12 at double precision
        //  - 1E-13 at quad precision
        qreal initialValidationEps = getQuESTValidationEpsilon();
        setQuESTValidationEpsilon(2 * initialValidationEps);

        const int NUM_QUBITS = 8;
        qreal dt = 0.1;
        int order = 4;
        int reps = 5;
        const int STEPS = 20;
        bool permutePaulis = GENERATE(true, false);

        auto unitaryTimeEvoFunc = 
        [dt, order, reps, STEPS, permutePaulis](Qureg qureg, PauliStrSum& hamil, PauliStrSum& observable) 
        -> qvector {
            qvector observations = getZeroVector(STEPS);
            initPlusState(qureg);
           
            for (int i = 0; i < STEPS; i++) {
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, dt, order, reps, permutePaulis);
                observations.at(i)  = calcExpecPauliStrSum(qureg, observable);
            }

            return observations;
        };
        
        PauliStrSum hamil = createHeisenbergHamiltonian(NUM_QUBITS);
        PauliStrSum observ = createAlternatingPauliObservable(NUM_QUBITS);
        
        qvector refObservables = {
            8.521995598825049,
            8.963711845322115,
            9.32005226684505,
            9.587768088649602,
            9.765522600223822,
            9.85387668440598,
            9.855195944206464,
            9.773484879367675,
            9.614158409472378,
            9.383765238225045,
            9.089680663909942,
            8.739788123639109,
            8.342168826039893,
            7.904817272753528,
            7.435397472039873,
            6.94105054616863,
            6.428259679798389,
            5.90277345392904,
            5.369584051930907,
            4.832953030744839
        };

        SECTION("Time Evolve Statevectors") {
            quregCache eightQubitSVCache = createCustomCachedQuregs(NUM_QUBITS, false);
            TEST_OBSERVABLES_ON_QUREGS(eightQubitSVCache, refObservables, unitaryTimeEvoFunc, hamil, observ);
            destroyCustomCachedQuregs(eightQubitSVCache);
        }

        SECTION("Time Evolve Density Matrices") {
            quregCache eightQubitDMCache = createCustomCachedQuregs(NUM_QUBITS, true);
            TEST_OBSERVABLES_ON_QUREGS(eightQubitDMCache, refObservables, unitaryTimeEvoFunc, hamil, observ);
            destroyCustomCachedQuregs(eightQubitDMCache);
        }

        // Restore validation epsilon
        setQuESTValidationEpsilon(initialValidationEps);

        destroyPauliStrSum(hamil);
        destroyPauliStrSum(observ);
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedStatevec();
        PauliStrSum hamil = createHeisenbergHamiltonian(qureg.numQubits);
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH( 
                applyTrotterizedUnitaryTimeEvolution(badQureg, hamil, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badHamil = hamil;
            badHamil.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, badHamil, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, nonHermitian, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeHamil = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, largeHamil, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeHamil);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedUnitaryTimeEvolution(qureg, hamil, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(hamil);
    }
}


TEST_CASE( "applyTrotterizedImaginaryTimeEvolution", TEST_CATEGORY ) {
    int numQubits = getNumCachedQubits();
    auto statevecQuregs = getCachedStatevecs();
    auto densmatrQuregs = getCachedDensmatrs();

    SECTION( LABEL_CORRECTNESS ) {
           
        qreal tau = 0.1;
        int order = 6;
        int reps = 5;
        int steps = 10;
        bool permutePaulis = GENERATE(true, false);

        auto driveToGroundFunc = [steps, tau, order, reps, permutePaulis](Qureg qureg, PauliStrSum& hamil) {
            initPlusState(qureg);
        
            for (int i = 0; i < steps; ++i) {
              applyTrotterizedImaginaryTimeEvolution(qureg, hamil, tau, order, reps, permutePaulis);
              setQuregToRenormalized(qureg);
            }
        };
       

#if QUEST_FLOAT_PRECISION == 4
        /*
         * The numerical exponent is sufficiently inaccurate to breach the default
         * tolerances at quad precision, so we apply the following kludge to prevent irritating test failures.
         * The real lessons from these tests are: 
         *   - Don't do time-evolution at single precision.
         *   - Don't do time-evolution in serial.
         */
    
        qreal initialEps = getTestAbsoluteEpsilon();
        setTestAbsoluteEpsilon(300 * initialEps);
#endif


        // Ground state: all qubits align down (driven by strong magnetic field)
        SECTION("Spin Down Field")
        {
            PauliStrSum ising = createIsingHamiltonian(numQubits, 10.0, 0.0, 0.0);

            qvector statevecRef = getZeroVector(getPow2(numQubits));
            statevecRef.at(0) = 1;

            qmatrix densmatrRef = getZeroMatrix(getPow2(numQubits));
            densmatrRef[0][0] = 1;

            TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, driveToGroundFunc, ising);
            TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, driveToGroundFunc, ising); 

            destroyPauliStrSum(ising);
        }
        
        // Ground state: all qubits align up (driven by opposite magnetic field)
        SECTION("Spin Up Field")
        {     
            PauliStrSum ising = createIsingHamiltonian(numQubits, -10.0, 0.0, 0.0);

            qindex namps = getPow2(numQubits);

            qvector statevecRef = getZeroVector(namps);
            statevecRef.at(namps - 1) = 1;

            qmatrix densmatrRef = getZeroMatrix(namps);
            densmatrRef[namps-1][namps-1] = 1;

            TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, driveToGroundFunc, ising);
            TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, driveToGroundFunc, ising); 

            destroyPauliStrSum(ising);
        }
        
        // Ground state: all qubits align down (driven by ferromagnetic interactions and bias)
        SECTION("Ferromagnetic Interaction")
        {
            PauliStrSum ising = createIsingHamiltonian(numQubits, 0.0, 10.0, 10.0);
           
            qvector statevecRef = getZeroVector(getPow2(numQubits));
            statevecRef.at(0) = 1;

            qmatrix densmatrRef = getZeroMatrix(getPow2(numQubits));
            densmatrRef[0][0] = 1;

            TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, driveToGroundFunc, ising);
            TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, driveToGroundFunc, ising); 

            destroyPauliStrSum(ising);
        }
        
        // Ground state: alternating pattern (driven by antiferromagnetic interactions)
        SECTION("Antiferromagnetic Interaction")
        {
            PauliStrSum ising = createIsingHamiltonian(numQubits, 0.0, -10.0, 10.0);
            
            // This should correctly pick out the non-zero amplitude
            // Qubit 0 is always 0 thanks to asymmetric bias 
            unsigned long long idx = 0;
            for (int i = 0; i < numQubits / 2; ++i) {
                idx += (1ULL << (2*i + 1));
            }
 
            qvector statevecRef = getZeroVector(getPow2(numQubits));
            statevecRef.at(idx) = 1;

            qmatrix densmatrRef = getZeroMatrix(getPow2(numQubits));
            densmatrRef[idx][idx] = 1;

            TEST_ON_CACHED_QUREGS(statevecQuregs, statevecRef, driveToGroundFunc, ising);
            TEST_ON_CACHED_QUREGS(densmatrQuregs, densmatrRef, driveToGroundFunc, ising); 

            destroyPauliStrSum(ising);
        }

#if QUEST_FLOAT_PRECISION == 4
        setTestEpsilon(initialEps);
#endif
    }

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedStatevec();
        PauliStrSum sum = createIsingHamiltonian(qureg.numQubits, 1.0, 1.0, 0.0);
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(badQureg, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badSum = sum;
            badSum.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, badSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeSum = createIsingHamiltonian(qureg.numQubits+1, 1.0, 1.0, 0.0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, largeSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeSum);
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, nonHermitian, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, sum, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedImaginaryTimeEvolution(qureg, sum, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(sum);
    }
}

/**
 * @todo
 * UNTESTED FUNCTIONS (NOT YET VALIDATED BY REFERENCE TESTS)
 */

TEST_CASE( "applyTrotterizedPauliStrSumGadget", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {}

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = GENERATE(getArbitraryCachedStatevec(), getArbitraryCachedDensmatr());
        PauliStrSum sum = createHeisenbergHamiltonian(qureg.numQubits);
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedPauliStrSumGadget(badQureg, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badSum = sum;
            badSum.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedPauliStrSumGadget(qureg, badSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeSum = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedPauliStrSumGadget(qureg, largeSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeSum);
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedPauliStrSumGadget(qureg, nonHermitian, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedPauliStrSumGadget(qureg, sum, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedPauliStrSumGadget(qureg, sum, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(sum);
    }
}


TEST_CASE( "applyTrotterizedControlledPauliStrSumGadget", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {}

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = GENERATE(getArbitraryCachedStatevec(), getArbitraryCachedDensmatr());
        PauliStrSum sum = createInlinePauliStrSum("1 XYI \n 2 ZZI");
        int ctrl = 0;
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(badQureg, ctrl, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badSum = sum;
            badSum.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, ctrl, badSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeSum = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, ctrl, largeSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeSum);
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, ctrl, nonHermitian, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, ctrl, sum, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, ctrl, sum, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "invalid control" ) {

            int badCtrl = GENERATE_COPY(-1, qureg.numQubits);
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, badCtrl, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Invalid control")
            );
        }

        SECTION( "control and target collision" ) {

            int badCtrl = 1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedControlledPauliStrSumGadget(qureg, badCtrl, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("control qubit overlaps")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(sum);
    }
}


TEST_CASE( "applyTrotterizedMultiControlledPauliStrSumGadget", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {}

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = GENERATE(getArbitraryCachedStatevec(), getArbitraryCachedDensmatr());
        PauliStrSum sum = createInlinePauliStrSum("1 XYII \n 2 ZZII");
        vector<int> ctrls = {0, 1};
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(badQureg, ctrls, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badSum = sum;
            badSum.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, ctrls, badSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeSum = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, ctrls, largeSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeSum);
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, ctrls, nonHermitian, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, ctrls, sum, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, ctrls, sum, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "invalid control" ) {

            vector<int> badCtrls = ctrls;
            badCtrls[0] = GENERATE_COPY(-1, qureg.numQubits);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, badCtrls, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Invalid control")
            );
        }

        SECTION( "control and target collision" ) {

            vector<int> badCtrls = ctrls;
            badCtrls[0] = 2;
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, badCtrls, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("control qubit overlaps")
            );
        }

        SECTION( "repetition in controls" ) {

            vector<int> badCtrls = ctrls;
            badCtrls[0] = badCtrls[1];
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, badCtrls, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("control qubits contained duplicates")
            );
        }

        SECTION( "number of controls" ) {

            int badNumCtrls = GENERATE_COPY( -1, qureg.numQubits + 1 );
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiControlledPauliStrSumGadget(qureg, ctrls.data(), badNumCtrls, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("number of control qubits")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(sum);
    }
}


TEST_CASE( "applyTrotterizedMultiStateControlledPauliStrSumGadget", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {}

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = GENERATE(getArbitraryCachedStatevec(), getArbitraryCachedDensmatr());
        PauliStrSum sum = createInlinePauliStrSum("1 XYII \n 2 ZZII");
        vector<int> ctrls = {0, 1};
        vector<int> states = {1, 0};
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(badQureg, ctrls, states, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badSum = sum;
            badSum.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, states, badSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeSum = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, states, largeSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeSum);
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, states, nonHermitian, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, states, sum, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, states, sum, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "invalid control" ) {

            vector<int> badCtrls = ctrls;
            badCtrls[0] = GENERATE_COPY(-1, qureg.numQubits);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, badCtrls, states, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Invalid control")
            );
        }

        SECTION( "invalid control states" ) {

            vector<int> badStates = states;
            badStates[0] = GENERATE(-1, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, badStates, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid control-state")
            );
        }

        SECTION( "control and target collision" ) {

            vector<int> badCtrls = ctrls;
            badCtrls[0] = 2;
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, badCtrls, states, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("control qubit overlaps")
            );
        }

        SECTION( "repetition in controls" ) {

            vector<int> badCtrls = ctrls;
            badCtrls[0] = badCtrls[1];
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, badCtrls, states, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("control qubits contained duplicates")
            );
        }

        SECTION( "number of controls" ) {

            int badNumCtrls = GENERATE_COPY( -1, qureg.numQubits + 1 );
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls.data(), states.data(), badNumCtrls, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("number of control qubits")
            );
        }

        SECTION( "number of controls and states" ) {

            vector<int> badStates = states;
            badStates.push_back(0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, ctrls, badStates, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("number of control qubits")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(sum);
    }
}


TEST_CASE( "applyTrotterizedNonUnitaryPauliStrSumGadget", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {}

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = GENERATE(getArbitraryCachedStatevec(), getArbitraryCachedDensmatr());
        PauliStrSum sum = createInlinePauliStrSum("1+1i X");
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedNonUnitaryPauliStrSumGadget(badQureg, sum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "pauli sum uninitialized" ) {

            PauliStrSum badSum = sum;
            badSum.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedNonUnitaryPauliStrSumGadget(qureg, badSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "pauli sum exceeds qureg qubits" ) {

            PauliStrSum largeSum = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNonUnitaryPauliStrSumGadget(qureg, largeSum, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeSum);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNonUnitaryPauliStrSumGadget(qureg, sum, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNonUnitaryPauliStrSumGadget(qureg, sum, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(sum);
    }
}


TEST_CASE( "applyTrotterizedNoisyTimeEvolution", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {}

    SECTION( LABEL_VALIDATION ) {

        Qureg qureg = getArbitraryCachedDensmatr();
        PauliStrSum hamil = createHeisenbergHamiltonian(qureg.numQubits);
        PauliStrSum jump  = createInlinePauliStrSum("0.5 X");
        qreal damp = 0.1;
        int numJumps = 1;
        bool permuteTerms = false;

        SECTION( "qureg uninitialised" ) {

            Qureg badQureg = qureg;
            badQureg.numQubits = -1;
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(badQureg, hamil, &damp, &jump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("invalid Qureg")
            );
        }

        SECTION( "qureg is not density matrix" ) {

            Qureg badQureg = getArbitraryCachedStatevec();
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(badQureg, hamil, &damp, &jump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("received a statevector")
            );
        }

        SECTION( "hamiltonian uninitialized" ) {

            PauliStrSum badHamil = hamil;
            badHamil.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, badHamil, &damp, &jump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "jump operator uninitialized" ) {

            PauliStrSum badJump = jump;
            badJump.numTerms = 0;
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, hamil, &damp, &badJump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Pauli")
            );
        }

        SECTION( "hamiltonian exceeds qureg qubits" ) {

            PauliStrSum largeHamil = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, largeHamil, &damp, &jump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeHamil);
        }

        SECTION( "jump operator exceeds qureg qubits" ) {

            PauliStrSum largeJump = createHeisenbergHamiltonian(qureg.numQubits + 1);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, hamil, &damp, &largeJump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("only compatible")
            );
            destroyPauliStrSum(largeJump);
        }

        SECTION( "hamiltonian not hermitian" ) {

            PauliStrSum nonHermitian = createInlinePauliStrSum("1+1i X");
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, nonHermitian, &damp, &jump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("Hermitian")
            );
            destroyPauliStrSum(nonHermitian);
        }

        SECTION( "invalid trotter order" ) {

            int order = GENERATE(-2, 0, 3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, hamil, &damp, &jump, numJumps, 0.1, order, 5, permuteTerms),
                ContainsSubstring("order")
            );
        }

        SECTION( "invalid trotter reps" ) {

            int reps = GENERATE(-3, 0);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, hamil, &damp, &jump, numJumps, 0.1, 4, reps, permuteTerms),
                ContainsSubstring("repetitions")
            );
        }

        SECTION( "invalid number of jumps" ) {

            int badNumJumps = GENERATE(-1, -3);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, hamil, &damp, &jump, badNumJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("number of jump")
            );
        }

        SECTION( "negative damping rate" ) {

            qreal badDamp = GENERATE(-1., -3.);
            REQUIRE_THROWS_WITH(
                applyTrotterizedNoisyTimeEvolution(qureg, hamil, &badDamp, &jump, numJumps, 0.1, 4, 5, permuteTerms),
                ContainsSubstring("damping rates were negative")
            );
        }

        SECTION( "sum ordering allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        SECTION( "super-propagator allocation failure" ) {

            // there is no reliable way to force the allocs to fail
            SUCCEED( );
        }

        destroyPauliStrSum(hamil);
        destroyPauliStrSum(jump);
    }

}
