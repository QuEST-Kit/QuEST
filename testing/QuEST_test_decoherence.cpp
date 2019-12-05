

#include "catch.hpp"
#include "QuEST.h"
#include "QuEST_test_utils.hpp"
#include <random>

/** The default number of qubits in the density matrices created for unit testing.
 * Creation of non-NUM_QUBITS sized Quregs should be justified in a comment. 
 * Note that the smaller this number is, the fewer nodes can be employed in 
 * distribution testing, since each node must contain at least one amplitude.
 * Furthermore, the larger this number is, the greater the deviation of correct 
 * results from their expected value, due to numerical error; this is especially 
 * apparent for density matrices.
 */
#define NUM_QUBITS 5

/** Prepares the needed data structures for unit testing. This creates 
 * the QuEST environment, a density matrix of the size 'numQb',
 * and corresponding QQMatrix instance for analytic comparison.
 * numQb should be NUM_QUBITS unless motivated otherwise.
 */
#define PREPARE_TEST(env, qureg, ref, numQb) \
    QuESTEnv env = createQuESTEnv(); \
    Qureg qureg = createDensityQureg(numQb, env); \
    initDebugState(qureg); \
    QMatrix ref = toQMatrix(qureg);

/** Destroys the data structures made by PREPARE_TEST */
#define CLEANUP_TEST(env, qureg) \
    destroyQureg(qureg, env); \
    destroyQuESTEnv(env);

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;



TEST_CASE( "mixDamping", "[decoherence]" ) {
    
    PREPARE_TEST(env, qureg, ref, NUM_QUBITS);
    qreal prob = getRandomReal(0, 1);

    SECTION( "correctness " ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        mixDamping(qureg, target, prob);
        
        // ref -> kraus0 ref kraus0^dagger + kraus1 ref kraus1^dagger
        QMatrix kraus0{{1,0},{0,sqrt(1-prob)}};
        QMatrix rho0 = ref;
        applyReferenceOp(rho0, target, kraus0);
        QMatrix kraus1{{0,sqrt(prob)},{0,0}};
        QMatrix rho1 = ref;
        applyReferenceOp(rho1, target, kraus1);
        ref = rho0 + rho1;
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "validation ") {
        
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixDamping(qureg, target, prob), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDamping(qureg, 0, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixDamping(qureg, 0, 1.1), Contains("Probabilities") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( mixDamping(vec, 0, 0), Contains("density matrices") );
            destroyQureg(vec, env);
        }
    }
    CLEANUP_TEST(env, qureg);
}



TEST_CASE( "mixDephasing", "[decoherence]" ) {
    
    PREPARE_TEST(env, qureg, ref, NUM_QUBITS);
    qreal prob = getRandomReal(0, 1/2.);

    SECTION( "correctness " ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        mixDephasing(qureg, target, prob);
        
        // ref -> (1 - prob) ref + prob Z ref Z
        QMatrix phaseRef = ref;
        applyReferenceOp(phaseRef, target, QMatrix{{1,0},{0,-1}}); // Z ref Z
        ref = ((1 - prob) * ref) + (prob * phaseRef);
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "validation ") {
        
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixDephasing(qureg, target, prob), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDephasing(qureg, 0, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixDephasing(qureg, 0, .6), Contains("probability") && Contains("cannot exceed 1/2") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( mixDephasing(vec, 0, 0), Contains("density matrices") );
            destroyQureg(vec, env);
        }
    }
    CLEANUP_TEST(env, qureg);
}



TEST_CASE( "mixDepolarising", "[decoherence]" ) {
    
    PREPARE_TEST(env, qureg, ref, NUM_QUBITS);
    qreal prob = getRandomReal(0, 3/4.);

    SECTION( "correctness " ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        mixDepolarising(qureg, target, prob);
        
        QMatrix xRef = ref;
        applyReferenceOp(xRef, target, QMatrix{{0,1},{1,0}}); // X ref X
        QMatrix yRef = ref;
        applyReferenceOp(yRef, target, QMatrix{{0,-1i},{1i,0}}); // Y ref Y
        QMatrix zRef = ref;
        applyReferenceOp(zRef, target, QMatrix{{1,0},{0,-1}}); // Z ref Z
        ref = ((1 - prob) * ref) + ((prob/3.) * ( xRef + yRef + zRef));
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "validation ") {
        
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, target, prob), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, 0, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, 0, .76), Contains("probability") && Contains("cannot exceed 3/4") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( mixDepolarising(vec, 0, 0), Contains("density matrices") );
            destroyQureg(vec, env);
        }
    }
    CLEANUP_TEST(env, qureg);
}



TEST_CASE( "mixPauli", "[decoherence]") {
    
    PREPARE_TEST(env, qureg, ref, NUM_QUBITS);

    // randomly generate valid pauli-error probabilities
    qreal probs[3];
    qreal max0 = 1/2.;                 // satisfies p1 < 1 - py
    probs[0] = getRandomReal(0, max0);
    qreal max1 = (max0 - probs[0])/2.; // p2 can use half of p1's "unused space"
    probs[1] = getRandomReal(0, max1);
    qreal max2 = (max1 - probs[1])/2.; // p3 can use half of p2's "unused space"
    probs[2] = getRandomReal(0, max2);
        
    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        
        // uniformly randomly assign probs (bound to target)
        int inds[3] = {0,1,2};
        std::shuffle(inds,inds+3, std::default_random_engine(1E5 * target));
        qreal probX = probs[inds[0]];           // seed:target shows no variation
        qreal probY = probs[inds[1]];
        qreal probZ = probs[inds[2]];

        mixPauli(qureg, target, probX, probY, probZ);
        
        QMatrix xRef = ref;
        applyReferenceOp(xRef, target, QMatrix{{0,1},{1,0}}); // X ref X
        QMatrix yRef = ref;
        applyReferenceOp(yRef, target, QMatrix{{0,-1i},{1i,0}}); // Y ref Y
        QMatrix zRef = ref;
        applyReferenceOp(zRef, target, QMatrix{{1,0},{0,-1}}); // Z ref Z
        ref = ((1 - probX - probY - probZ) * ref) +
              (probX * xRef) + (probY * yRef) + (probZ * zRef);
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "validation" ) {
        
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target, 0, 0, 0), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
                
            int target = 0;
            
            // probs clearly must be in [0, 1]
            REQUIRE_THROWS_WITH( mixPauli(qureg, target, -.1,   0,   0), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,   0, -.1,   0), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,   0,   0, -.1), Contains("Probabilities") );
            
            // max single-non-zero-prob is 0.5
            REQUIRE_THROWS_WITH( mixPauli(qureg, target, .6,  0,  0), Contains("cannot exceed the probability") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,  0, .6,  0), Contains("cannot exceed the probability") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,  0,  0, .6), Contains("cannot exceed the probability") );
            
            // must satisfy px, py, pz < 1 - px - py - pz
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,  .3,  .3, .3), Contains("cannot exceed the probability") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( mixPauli(vec, 0, 0, 0, 0), Contains("density matrices") );
            destroyQureg(vec, env);
        }
    }
    CLEANUP_TEST(env, qureg);
}



TEST_CASE( "mixKrausMap" ) {
    
    PREPARE_TEST(env, qureg, ref, NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        int numOps = GENERATE( range(1,5) ); // max 4 inclusive
        std::vector<QMatrix> matrs = getRandomKrausMap(1, numOps);
        
        ComplexMatrix2 ops[numOps];
        for (int i=0; i<numOps; i++)
            ops[i] = toComplexMatrix2(matrs[i]);
        mixKrausMap(qureg, target, ops, numOps);
        
        // set ref -> K_i ref K_i^dagger
        QMatrix matrRefs[numOps];
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], target, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref) );
        
    }
    SECTION( "validation" ) {
        
        SECTION( "number of operators" ) {
            
            int numOps = GENERATE( 0, 5 );
            REQUIRE_THROWS_WITH( mixKrausMap(qureg, 0, NULL, numOps), Contains("operators") );
        }
        SECTION( "trace preserving" ) {
            
            // valid Kraus map
            int numOps = GENERATE( range(1,5) ); // max 4 inclusive
            std::vector<QMatrix> matrs = getRandomKrausMap(1, numOps);
            ComplexMatrix2 ops[numOps];
            for (int i=0; i<numOps; i++)
                ops[i] = toComplexMatrix2(matrs[i]);
                
            // make invalid
            ops[GENERATE_REF( range(0,numOps) )].real[0][0] = 0;
            REQUIRE_THROWS_WITH( mixKrausMap(qureg, 0, ops, numOps), Contains("trace preserving") );
            
        }
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixKrausMap(qureg, target, NULL, 1), Contains("Invalid target qubit") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( mixKrausMap(vec, 0, NULL, 1), Contains("density matrices") );
            destroyQureg(vec, env);
        }
        SECTION( "operators fit in node" ) {
            
            qureg.numAmpsPerChunk = 3; // min 4
            REQUIRE_THROWS_WITH( mixKrausMap(qureg, 0, NULL, 1), Contains("targets too many qubits") );
        }        
    }
    CLEANUP_TEST(env, qureg);
}
