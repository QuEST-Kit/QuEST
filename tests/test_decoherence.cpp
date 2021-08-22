
#include "catch.hpp"
#include "QuEST.h"
#include "utilities.hpp"
#include <random>

/** Prepares a density matrix in the debug state, and the reference QMatrix 
 */
#define PREPARE_TEST(qureg, ref) \
    Qureg qureg = createDensityQureg(NUM_QUBITS, QUEST_ENV); \
    initDebugState(qureg); \
    QMatrix ref = toQMatrix(qureg);

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;



/** @sa mixDamping
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDamping", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);

    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        qreal prob = getRandomReal(0, 1);
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
            REQUIRE_THROWS_WITH( mixDamping(qureg, target, 0), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDamping(qureg, 0, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixDamping(qureg, 0, 1.1), Contains("Probabilities") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixDamping(vec, 0, 0), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixDensityMatrix
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDensityMatrix", "[decoherence]" ) {
    
    Qureg qureg1 = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    Qureg qureg2 = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    initDebugState(qureg1);
    initDebugState(qureg2);
    QMatrix ref1 = toQMatrix(qureg1);
    QMatrix ref2 = toQMatrix(qureg2);
    
    SECTION( "correctness" ) {
        
        // test p in {0, 1} and 10 random values in (0,1)
        qreal prob = GENERATE( 0., 1., take(10, random(0.,1.)) );
        mixDensityMatrix(qureg1, prob, qureg2);
        
        // ensure target qureg modified correctly
        ref1 = (1-prob)*ref1 + (prob)*ref2;
        REQUIRE( areEqual(qureg1, ref1) );
        
        // enure other qureg was not modified
        REQUIRE( areEqual(qureg2, ref2) );
    }
    SECTION( "input validation" ) {
        
        SECTION( "probabilities") {
            
            qreal prob = GENERATE( -0.1, 1.1 );
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg1, prob, qureg2), Contains("Probabilities") );
        }
        SECTION( "density matrices" ) {
            
            // one is statevec 
            Qureg state1 = createQureg(qureg1.numQubitsRepresented, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg1, 0, state1), Contains("density matrices") );
            REQUIRE_THROWS_WITH( mixDensityMatrix(state1, 0, qureg1), Contains("density matrices") );
            
            // both are statevec
            Qureg state2 = createQureg(qureg1.numQubitsRepresented, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixDensityMatrix(state1, 0, state2), Contains("density matrices") );
            
            destroyQureg(state1, QUEST_ENV);
            destroyQureg(state2, QUEST_ENV);
        }
        SECTION( "matching dimensions" ) {
            
            Qureg qureg3 = createDensityQureg(1 + qureg1.numQubitsRepresented, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg1, 0, qureg3), Contains("Dimensions") );
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg3, 0, qureg1), Contains("Dimensions") );
            destroyQureg(qureg3, QUEST_ENV);
        }
    }
    destroyQureg(qureg1, QUEST_ENV);
    destroyQureg(qureg2, QUEST_ENV);
}



/** @sa mixDephasing
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDephasing", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);

    SECTION( "correctness " ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        qreal prob = getRandomReal(0, 1/2.);
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
            REQUIRE_THROWS_WITH( mixDephasing(qureg, target, 0), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDephasing(qureg, 0, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixDephasing(qureg, 0, .6), Contains("probability") && Contains("cannot exceed 1/2") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixDephasing(vec, 0, 0), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixDepolarising
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDepolarising", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);

    SECTION( "correctness " ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        qreal prob = getRandomReal(0, 3/4.);
        mixDepolarising(qureg, target, prob);
        
        QMatrix xRef = ref;
        applyReferenceOp(xRef, target, QMatrix{{0,1},{1,0}}); // X ref X
        QMatrix yRef = ref;
        applyReferenceOp(yRef, target, QMatrix{{0,-qcomp(0,1)},{qcomp(0,1),0}}); // Y ref Y
        QMatrix zRef = ref;
        applyReferenceOp(zRef, target, QMatrix{{1,0},{0,-1}}); // Z ref Z
        ref = ((1 - prob) * ref) + ((prob/3.) * ( xRef + yRef + zRef));
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "validation ") {
        
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, target, 0), Contains("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, 0, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, 0, .76), Contains("probability") && Contains("cannot exceed 3/4") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixDepolarising(vec, 0, 0), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixMultiQubitKrausMap
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixMultiQubitKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    // figure out max-num (inclusive) targs allowed by hardware backend
    // (each node must contain as 2^(2*numTargs) amps)
    int maxNumTargs = calcLog2(qureg.numAmpsPerChunk) / 2;
    
    SECTION( "correctness" ) {
        
        /* note that this function incurs a stack overhead when numTargs < 4,
         * and a heap overhead when numTargs >= 4
         */
         
        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) ); // inclusive upper bound
        
        // note this is very expensive to try every arrangement (2 min runtime for numTargs=5 alone)
        int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        
        // try the min and max number of operators, and 2 random numbers 
        // (there are way too many to try all!)
        int maxNumOps = (2*numTargs)*(2*numTargs);
        int numOps = GENERATE_COPY( 1, maxNumOps, take(2,random(1,maxNumOps)) );
        
        // use a new random map
        std::vector<QMatrix> matrs = getRandomKrausMap(numTargs, numOps);
                
        // create map in QuEST datatypes
        ComplexMatrixN ops[numOps];
        for (int i=0; i<numOps; i++) {
            ops[i] = createComplexMatrixN(numTargs);
            toComplexMatrixN(matrs[i], ops[i]);
        }
                
        mixMultiQubitKrausMap(qureg, targs, numTargs, ops, numOps);
                
        // set ref -> K_i ref K_i^dagger
        QMatrix matrRefs[numOps];
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], targs, numTargs, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref, 1E2*REAL_EPS) );
        
        // cleanup QuEST datatypes
        for (int i=0; i<numOps; i++)
            destroyComplexMatrixN(ops[i]);
    }
    SECTION( "input validation" ) {
        
        SECTION( "repetition of target" ) {
            
            // make valid targets
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
                
            // duplicate one
            int badInd = GENERATE( range(0,NUM_QUBITS) );
            int copyInd = GENERATE_COPY( filter([=](int i){ return i!=badInd; }, range(0,NUM_QUBITS)) );
            targs[badInd] = targs[copyInd];
            
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, NUM_QUBITS, NULL, 1), Contains("target qubits") && Contains("unique") );
        }
        SECTION( "qubit indices" ) { 
            
            // make valid targets
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
            
            // make one invalid 
            targs[GENERATE( range(0,NUM_QUBITS) )] = GENERATE( -1, NUM_QUBITS );
            
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, NUM_QUBITS, NULL, 1), Contains("Invalid target qubit") );
        }
        SECTION( "number of operators" ) {
            
            int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            int maxNumOps = (2*numTargs)*(2*numTargs);
            int numOps = GENERATE_REF( -1, 0, maxNumOps + 1 );
                        
            // make valid targets to avoid triggering target validation
            int targs[numTargs];
            for (int i=0; i<numTargs; i++)
                targs[i] = i;
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, numTargs, NULL, numOps), Contains("operators may be specified") );
        }
        SECTION( "initialisation of operators" ) {
            
            /* compilers don't auto-initialise to NULL; the below circumstance 
             * only really occurs when 'malloc' returns NULL in createComplexMatrixN, 
             * which actually triggers its own validation. Hence this test is useless 
             * currently.
             */
             
            int numTargs = NUM_QUBITS;
            int numOps = (2*numTargs)*(2*numTargs);
            
            // no need to initialise ops, but set their attribs correct to avoid triggering other validation
            ComplexMatrixN ops[numOps];
            for (int i=0; i<numOps; i++)
                ops[i].numQubits = numTargs;
            
            // make one of the max-ops explicitly NULL
            ops[GENERATE_COPY( range(0,numTargs) )].real = NULL;
             
            // make valid targets to avoid triggering target validation
            int targs[numTargs];
            for (int i=0; i<numTargs; i++)
                targs[i] = i;
               
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, numTargs, ops, numOps), Contains("ComplexMatrixN") && Contains("created") );
        }
        SECTION( "dimension of operators" ) {
        
            // make valid (dimension-wise) max-qubits Kraus map
            int numTargs = NUM_QUBITS;
            int numOps = (2*numTargs)*(2*numTargs);
            ComplexMatrixN ops[numOps];
            for (int i=0; i<numOps; i++)
                ops[i] = createComplexMatrixN(numTargs);
            
            // make one have wrong-dimensions 
            int badInd = GENERATE_COPY( range(0,numTargs) );
            destroyComplexMatrixN(ops[badInd]);
            ops[badInd] = createComplexMatrixN(numTargs - 1);
            
            // make valid targets to avoid triggering target validation
            int targs[numTargs];
            for (int i=0; i<numTargs; i++)
                targs[i] = i;
                
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, numTargs, ops, numOps), Contains("same number of qubits") );
            
            for (int i=0; i<numOps; i++)
                destroyComplexMatrixN(ops[i]);
        }
        SECTION( "trace preserving" ) {
            
            int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            int maxNumOps = (2*numTargs) * (2*numTargs);
            int numOps = GENERATE_COPY( 1, 2, maxNumOps );
            
            // generate a valid map
            std::vector<QMatrix> matrs = getRandomKrausMap(numTargs, numOps);
            ComplexMatrixN ops[numOps];
            for (int i=0; i<numOps; i++) {
                ops[i] = createComplexMatrixN(numTargs);
                toComplexMatrixN(matrs[i], ops[i]);
            }
            
            // make only one invalid
            ops[GENERATE_REF( range(0,numOps) )].real[0][0] = 0;
            
            // make valid targets to avoid triggering target validation
            int targs[numTargs];
            for (int i=0; i<numTargs; i++)
                targs[i] = i;
            
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, numTargs, ops, numOps), Contains("trace preserving") );

            for (int i=0; i<numOps; i++)
                destroyComplexMatrixN(ops[i]);
        }
        SECTION( "density-matrix" ) {
            
            Qureg statevec = createQureg(NUM_QUBITS, QUEST_ENV);
            
            // make valid targets to avoid triggering target validation
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
                
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(statevec, targs, NUM_QUBITS, NULL, 1), Contains("valid only for density matrices") );
            destroyQureg(statevec, QUEST_ENV);
            
        }
        SECTION( "operator fits in node" ) {
            
            // each node requires (2 numTargs)^2 amplitudes
            int minAmps = (2*NUM_QUBITS) * (2*NUM_QUBITS);
            
            // make valid targets to avoid triggering target validation
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
            
            // make a simple Identity map    
            ComplexMatrixN ops[] = {createComplexMatrixN(NUM_QUBITS)};
            for (int i=0; i<(1<<NUM_QUBITS); i++)
                ops[0].real[i][i] = 1;
            
            // fake a smaller qureg 
            qureg.numAmpsPerChunk = minAmps - 1;
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, NUM_QUBITS, ops, 1), Contains("targets too many qubits") && Contains("cannot all fit") );
            
            destroyComplexMatrixN(ops[0]);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixPauli
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixPauli", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
        
    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        
        // randomly generate valid pauli-error probabilities
        qreal probs[3];
        qreal max0 = 1/2.;                 // satisfies p1 < 1 - py
        probs[0] = getRandomReal(0, max0);
        qreal max1 = (max0 - probs[0])/2.; // p2 can use half of p1's "unused space"
        probs[1] = getRandomReal(0, max1);
        qreal max2 = (max1 - probs[1])/2.; // p3 can use half of p2's "unused space"
        probs[2] = getRandomReal(0, max2);
        
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
        applyReferenceOp(yRef, target, QMatrix{{0,-qcomp(0,1)},{qcomp(0,1),0}}); // Y ref Y
        QMatrix zRef = ref;
        applyReferenceOp(zRef, target, QMatrix{{1,0},{0,-1}}); // Z ref Z
        ref = ((1 - probX - probY - probZ) * ref) +
              (probX * xRef) + (probY * yRef) + (probZ * zRef);
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "input validation" ) {
        
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
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixPauli(vec, 0, 0, 0, 0), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixKrausMap
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
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
        
        REQUIRE( areEqual(qureg, ref, 10*REAL_EPS) );
    }
    SECTION( "input validation" ) {
        
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
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixKrausMap(vec, 0, NULL, 1), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
        SECTION( "operators fit in node" ) {
            
            qureg.numAmpsPerChunk = 3; // min 4
            REQUIRE_THROWS_WITH( mixKrausMap(qureg, 0, NULL, 1), Contains("targets too many qubits") );
        }        
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixTwoQubitDephasing
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixTwoQubitDephasing", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int targ1 = GENERATE( range(0,NUM_QUBITS) );
        int targ2 = GENERATE_COPY( filter([=](int t){ return t!=targ1; }, range(0,NUM_QUBITS)) );
        qreal prob = getRandomReal(0, 3/4.);
        
        mixTwoQubitDephasing(qureg, targ1, targ2, prob);
        
        // ref -> (1 - prob) ref + prob/3 (Z1 ref Z1 + Z2 ref Z2 + Z1 Z2 ref Z1 Z2)
        QMatrix zMatr{{1,0},{0,-1}};
        QMatrix z1Ref = ref;
        applyReferenceOp(z1Ref, targ1, zMatr); // Z1 ref Z1
        QMatrix z2Ref = ref;
        applyReferenceOp(z2Ref, targ2, zMatr); // Z2 ref Z2
        QMatrix z1z2Ref = ref;
        applyReferenceOp(z1z2Ref, targ1, zMatr);
        applyReferenceOp(z1z2Ref, targ2, zMatr); // Z1 Z2 ref Z1 Z2
        ref = ((1 - prob) * ref) + (prob/3.) * (z1Ref + z2Ref + z1z2Ref);
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int targ = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, 0, targ, 0), Contains("Invalid target") );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, targ, 0, 0), Contains("Invalid target") );
        }
        SECTION( "target collision" ) {
            
            int targ = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, targ, targ, 0), Contains("target") && Contains("unique") );
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, 0, 1, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, 0, 1, 3/4. + .01), Contains("probability") && Contains("cannot exceed 3/4") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(vec, 0, 1, 0), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixTwoQubitDepolarising
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixTwoQubitDepolarising", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int targ1 = GENERATE( range(0,NUM_QUBITS) );
        int targ2 = GENERATE_COPY( filter([=](int t){ return t!=targ1; }, range(0,NUM_QUBITS)) );
        qreal prob = getRandomReal(0, 15/16.);
        
        mixTwoQubitDepolarising(qureg, targ1, targ2, prob);
        
        QMatrix paulis[4] = {
            QMatrix{{1,0},{0,1}},       // I 
            QMatrix{{0,1},{1,0}},       // X
            QMatrix{{0,-qcomp(0,1)},{qcomp(0,1),0}},    // Y
            QMatrix{{1,0},{0,-1}}       // Z
        };
        
        int targs[2] = {targ1, targ2};
        QMatrix refInit = ref;
        ref = (1 - (16/15.)*prob) * ref;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                QMatrix term = refInit;
                applyReferenceOp(term, targs, 2,
                    getKroneckerProduct(paulis[i], paulis[j]));
                ref += (prob/15.) * term;
            }
        }
        
        REQUIRE( areEqual(qureg, ref, 1E4*REAL_EPS) );
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int targ = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, 0, targ, 0), Contains("Invalid target") );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, targ, 0, 0), Contains("Invalid target") );
        }
        SECTION( "target collision" ) {
            
            int targ = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, targ, targ, 0), Contains("target") && Contains("unique") );
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, 0, 1, -.1), Contains("Probabilities") );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, 0, 1, 15/16. + .01), Contains("probability") && Contains("cannot exceed 15/16") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(vec, 0, 1, 0), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
    }
    destroyQureg(qureg, QUEST_ENV);
}



/** @sa mixTwoQubitKrausMap
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "mixTwoQubitKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int targ1 = GENERATE( range(0,NUM_QUBITS) );
        int targ2 = GENERATE_COPY( filter([=](int t){ return t!=targ1; }, range(0,NUM_QUBITS)) );
        int numOps = GENERATE( range(1,17) ); // max 16 inclusive
        std::vector<QMatrix> matrs = getRandomKrausMap(2, numOps);
        
        ComplexMatrix4 ops[numOps];
        for (int i=0; i<numOps; i++)
            ops[i] = toComplexMatrix4(matrs[i]);
        mixTwoQubitKrausMap(qureg, targ1, targ2, ops, numOps);
        
        // set ref -> K_i ref K_i^dagger
        int targs[2] = {targ1, targ2};
        QMatrix matrRefs[numOps];
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], targs, 2, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref, 10*REAL_EPS) );
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of operators" ) {
            
            int numOps = GENERATE( 0, 17 );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,1, NULL, numOps), Contains("operators") );
        }
        SECTION( "trace preserving" ) {
            
            // valid Kraus map
            int numOps = GENERATE( range(1,16) );
            std::vector<QMatrix> matrs = getRandomKrausMap(2, numOps);
            ComplexMatrix4 ops[numOps];
            for (int i=0; i<numOps; i++)
                ops[i] = toComplexMatrix4(matrs[i]);
                
            // make only one of the ops at a time invalid
            ops[GENERATE_REF( range(0,numOps) )].real[0][0] = 0;
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,1, ops, numOps), Contains("trace preserving") );
        }
        SECTION( "target collision" ) {
            
            int target = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, target, target, NULL, 1), Contains("target qubits") && Contains("unique") );
        }
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,target, NULL, 1), Contains("Invalid target qubit") );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, target,0, NULL, 1), Contains("Invalid target qubit") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(vec, 0,1, NULL, 1), Contains("density matrices") );
            destroyQureg(vec, QUEST_ENV);
        }
        SECTION( "operators fit in node" ) {
            
            qureg.numAmpsPerChunk = 15; // min 16
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,1, NULL, 1), Contains("targets too many qubits") );
        }        
    }
    destroyQureg(qureg, QUEST_ENV);
}

