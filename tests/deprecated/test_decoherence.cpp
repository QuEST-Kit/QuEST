/** @file
 * Ported tests of the deprecated QuEST v3 interface,
 * unit testing the "decoherence" module.
 * 
 * This file should be excluded from doxygen parsing so 
 * as not to conflict with the doc of the v4 unit tests.
 * 
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_random.hpp>

// must define preprocessors to enable quest's
// deprecated v3 API, and disable the numerous
// warnings issued by its compilation
#define INCLUDE_DEPRECATED_FUNCTIONS 1
#define DISABLE_DEPRECATION_WARNINGS 1

#include "test_utilities.hpp"

#include <random>
#include <vector>
#include <algorithm>

using std::vector;

/** Prepares a density matrix in the debug state, and the reference QMatrix 
 */
#define PREPARE_TEST(qureg, ref) \
    Qureg qureg = createForcedDensityQureg(NUM_QUBITS); \
    initDebugState(qureg); \
    QMatrix ref = toQMatrix(qureg); \
    assertQuregAndRefInDebugState(qureg, ref); \
    setValidationEpsilon(REAL_EPS);

/* allows concise use of ContainsSubstring in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::ContainsSubstring;



/** @sa mixDamping
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDamping", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);

    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0, NUM_QUBITS) );
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
            REQUIRE_THROWS_WITH( mixDamping(qureg, target, 0), ContainsSubstring("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDamping(qureg, 0, -.1), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixDamping(qureg, 0, 1.1), ContainsSubstring("probability is invalid") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_NOTHROW( mixDamping(vec, 0, 0) ); // zero-prob allowed in v4
            REQUIRE_THROWS_WITH( mixDamping(vec, 0, .1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixDensityMatrix
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDensityMatrix", "[decoherence]" ) {
    
    Qureg qureg1 = createForcedDensityQureg(NUM_QUBITS);
    Qureg qureg2 = createForcedDensityQureg(NUM_QUBITS);
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
        
        SECTION( "probabilities" ) {
            
            qreal prob = GENERATE( -0.1, 1.1 );
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg1, prob, qureg2), ContainsSubstring("probability is invalid") );
        }
        SECTION( "density matrices" ) {
            
            // one is statevec 
            Qureg state1 = createQureg(qureg1.numQubits, getQuESTEnv());
            REQUIRE_THROWS_WITH( mixDensityMatrix(state1, 0, qureg1), ContainsSubstring("The first Qureg") && ContainsSubstring("must be a density matrix") );

            // in v4, the 2nd arg can be a statevector

                // REQUIRE_THROWS_WITH( mixDensityMatrix(qureg1, 0, state1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            
            // both are statevec
            Qureg state2 = createQureg(qureg1.numQubits, getQuESTEnv());
            REQUIRE_THROWS_WITH( mixDensityMatrix(state1, 0, state2), ContainsSubstring("The first Qureg") && ContainsSubstring("must be a density matrix") );
            
            destroyQureg(state1, getQuESTEnv());
            destroyQureg(state2, getQuESTEnv());
        }
        SECTION( "matching dimensions" ) {
            
            Qureg qureg3 = createDensityQureg(1 + qureg1.numQubits, getQuESTEnv());
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg1, 0, qureg3), ContainsSubstring("inconsistent number of qubits") );
            REQUIRE_THROWS_WITH( mixDensityMatrix(qureg3, 0, qureg1), ContainsSubstring("inconsistent number of qubits") );
            destroyQureg(qureg3, getQuESTEnv());
        }
    }
    destroyQureg(qureg1, getQuESTEnv());
    destroyQureg(qureg2, getQuESTEnv());
}



/** @sa mixDephasing
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixDephasing", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);

    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        qreal prob = getRandomReal(0, 1/2.);
        mixDephasing(qureg, target, prob);
        
        // ref -> (1 - prob) ref + prob Z ref Z
        QMatrix phaseRef = ref;
        applyReferenceOp(phaseRef, target, QMatrix{{1,0},{0,-1}}); // Z ref Z
        ref = ((1 - prob) * ref) + (prob * phaseRef);
        
        REQUIRE( areEqual(qureg, ref) );
    }
    SECTION( "validation" ) {
        
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixDephasing(qureg, target, 0), ContainsSubstring("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDephasing(qureg, 0, -.1), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixDephasing(qureg, 0, .6), ContainsSubstring("probability") && ContainsSubstring("1/2") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_NOTHROW( mixDephasing(vec, 0, 0) ); // zero prob ok in v4
            REQUIRE_THROWS_WITH( mixDephasing(vec, 0, 0.1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixDepolarising
 * @ingroup deprecatedtests 
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
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, target, 0), ContainsSubstring("Invalid target") );
            
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, 0, -.1), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixDepolarising(qureg, 0, .76), ContainsSubstring("probability") && ContainsSubstring("3/4") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_NOTHROW( mixDepolarising(vec, 0, 0) ); // zero-prob ok in v4
            REQUIRE_THROWS_WITH( mixDepolarising(vec, 0, 0.1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixMultiQubitKrausMap
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixMultiQubitKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    // figure out max-num (inclusive) targs allowed by hardware backend
    // (each node must contain as 2^(2*numTargs) amps)
    int maxNumTargs = qureg.logNumAmpsPerNode / 2;
    
    SECTION( "correctness" ) {
        
        /* note that this function incurs a stack overhead when numTargs < 4,
         * and a heap overhead when numTargs >= 4
         */
         
        // try every size (qubit wise) map
        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) ); // inclusive upper bound
        
        // previously, we tried every unique set of targets, via:
        //      int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        // alas, this is too slow for CI, so we instead try a fixed number of random sets:
        GENERATE( range(0, 10) );
        vector<int> targs(numTargs);
        setRandomTargets(targs, NUM_QUBITS);
        
        // try the min and max number of operators, and 2 random numbers 
        // (there are way too many to try all!)
        int maxNumOps = (2*numTargs)*(2*numTargs);
        int numOps = GENERATE_COPY( 1, maxNumOps, take(2,random(1,maxNumOps)) );
        
        // use a new random map
        vector<QMatrix> matrs = getRandomKrausMap(numTargs, numOps);
                
        // create map in QuEST datatypes
        vector<ComplexMatrixN> ops(numOps);
        for (int i=0; i<numOps; i++) {
            ops[i] = createCompMatr(numTargs);
            setCompMatr(ops[i], matrs[i]);
        }
                
        mixMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps);
        
        // set ref -> K_i ref K_i^dagger
        vector<QMatrix> matrRefs(numOps);
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], targs.data(), numTargs, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref, 1E2*REAL_EPS) );
        
        // cleanup QuEST datatypes
        for (int i=0; i<numOps; i++)
            destroyCompMatr(ops[i]);
    }
    SECTION( "input validation" ) {

        // spoof a list of properly-initialised CompMatr to avoid seg-fault in deprecation API.
        // note some functions will need smaller matrices which is fine; a subset of each
        // spoof'd matrix will be copied over by the deprecation layer
        ComplexMatrixN spoofOps[NUM_QUBITS];
        for (int i=0; i<NUM_QUBITS; i++) {
            spoofOps[i] = createComplexMatrixN(NUM_QUBITS);
            syncCompMatr(spoofOps[i]);
        }
        
        SECTION( "repetition of target" ) {
            
            // make valid targets
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
                
            // duplicate one
            int badInd = GENERATE( range(0,NUM_QUBITS) );
            int copyInd = GENERATE_COPY( filter([=](int i){ return i!=badInd; }, range(0,NUM_QUBITS)) );
            targs[badInd] = targs[copyInd];
            
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, NUM_QUBITS, spoofOps, 1), ContainsSubstring("target qubits") && ContainsSubstring("unique") );
        }
        SECTION( "qubit indices" ) { 
            
            // make valid targets
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
            
            // make one invalid 
            targs[GENERATE( range(0,NUM_QUBITS) )] = GENERATE( -1, NUM_QUBITS );
            
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, NUM_QUBITS, spoofOps, 1), ContainsSubstring("Invalid target qubit") );
        }

        // cannot test this with deprecated API, since 'numOps' informs a copy before validation

            // SECTION( "number of operators" ) {
                
            //     int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            //     int numOps = GENERATE_REF( -1, 0 );
                            
            //     // make valid targets to avoid triggering target validation
            //     vector<int> targs(numTargs);
            //     for (int i=0; i<numTargs; i++)
            //         targs[i] = i;

            //     REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs.data(), numTargs, spoofOps, numOps), ContainsSubstring("must be given a strictly positive number of matrices") );
            // }

        // this validation cannot be performed by the deprecation API which copies
        // over all ops arguments without prior checking them as NULL

            // SECTION( "initialisation of operators" ) {
                
            //     /* compilers don't auto-initialise to NULL; the below circumstance 
            //      * only really occurs when 'malloc' returns NULL in createCompMatr, 
            //      * which actually triggers its own validation. Hence this test is useless 
            //      * currently.
            //      */
                
            //     int numTargs = NUM_QUBITS;
            //     int numOps = (2*numTargs)*(2*numTargs);
                
            //     vector<ComplexMatrixN> ops(numOps);
            //     for (int i=0; i<numOps; i++)
            //         ops[i] = createComplexMatrixN(numTargs);
                
            //     // make one of the max-ops explicitly NULL
            //     ops[GENERATE_COPY( range(0,numTargs) )].cpuElems = NULL;
                
            //     // make valid targets to avoid triggering target validation
            //     vector<int> targs(numTargs);
            //     for (int i=0; i<numTargs; i++)
            //         targs[i] = i;
                
            //     REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps), ContainsSubstring("ComplexMatrixN") && ContainsSubstring("created") );

            //     for (int i=0; i<numOps; i++)
            //         destroyComplexMatrixN(ops[i]);
            // }

        // this validation cannot be performed by the deprecation API which copies
        // over all ops arguments without prior checking them as matching in dimension

            // SECTION( "dimension of operators" ) {
            
            //     // make valid (dimension-wise) max-qubits Kraus map
            //     int numTargs = NUM_QUBITS;
            //     int numOps = (2*numTargs)*(2*numTargs);
            //     vector<ComplexMatrixN> ops(numOps);
            //     for (int i=0; i<numOps; i++)
            //         ops[i] = createComplexMatrixN(numTargs);
                
            //     // make one have wrong-dimensions 
            //     int badInd = GENERATE_COPY( range(0,numTargs) );
            //     destroyComplexMatrixN(ops[badInd]);
            //     ops[badInd] = createComplexMatrixN(numTargs - 1);
                
            //     // make valid targets to avoid triggering target validation
            //     vector<int> targs(numTargs);
            //     for (int i=0; i<numTargs; i++)
            //         targs[i] = i;
                    
            //     REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps), ContainsSubstring("same number of qubits") );
                
            //     for (int i=0; i<numOps; i++)
            //         destroyComplexMatrixN(ops[i]);
            // }

        SECTION( "trace preserving" ) {
            
            int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            int maxNumOps = (2*numTargs) * (2*numTargs);
            int numOps = GENERATE_COPY( 1, 2, maxNumOps );
            
            // generate a valid map
            vector<QMatrix> matrs = getRandomKrausMap(numTargs, numOps);
            vector<ComplexMatrixN> ops(numOps);
            for (int i=0; i<numOps; i++) {
                ops[i] = createComplexMatrixN(numTargs);
                toComplexMatrixN(matrs[i], ops[i]);
            }
            
            // make only one invalid
            ops[GENERATE_COPY( 0, numOps - 1)].cpuElems[0][0] = -123456789;
            
            // make valid targets to avoid triggering target validation
            vector<int> targs(numTargs);
            for (int i=0; i<numTargs; i++)
                targs[i] = i;
            
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps), ContainsSubstring("trace preserving") );

            for (int i=0; i<numOps; i++)
                destroyComplexMatrixN(ops[i]);
        }
        SECTION( "density-matrix" ) {
            
            Qureg statevec = createQureg(NUM_QUBITS);
            
            // make valid targets to avoid triggering target validation
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
                
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(statevec, targs, NUM_QUBITS, spoofOps, 1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(statevec, getQuESTEnv());
            
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
                ops[0].cpuElems[i][i] = 1;
            syncCompMatr(ops[0]);
            
            // fake a smaller qureg 
            qureg.isDistributed = 1;
            qureg.numAmpsPerNode = minAmps - 1;
            qureg.logNumAmpsPerNode = log2(minAmps) - 1;
            REQUIRE_THROWS_WITH( mixMultiQubitKrausMap(qureg, targs, NUM_QUBITS, ops, 1), ContainsSubstring("each node's communication buffer") && ContainsSubstring("cannot simultaneously store") );
            
            destroyComplexMatrixN(ops[0]);
        }

        for (int i=0; i<NUM_QUBITS; i++)
            destroyCompMatr(spoofOps[i]);
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixPauli
 * @ingroup deprecatedtests 
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
            REQUIRE_THROWS_WITH( mixPauli(qureg, target, 0, 0, 0), ContainsSubstring("Invalid target") );
            
        }
        SECTION( "probability" ) {
                
            int target = 0;
            
            // probs clearly must be in [0, 1]
            REQUIRE_THROWS_WITH( mixPauli(qureg, target, -.1,   0,   0), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,   0, -.1,   0), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,   0,   0, -.1), ContainsSubstring("probability is invalid") );
            
            // max single-non-zero-prob is 0.5
            REQUIRE_THROWS_WITH( mixPauli(qureg, target, .6,  0,  0), ContainsSubstring("probabilities exceed that which induce maximal mixing") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,  0, .6,  0), ContainsSubstring("probabilities exceed that which induce maximal mixing") );
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,  0,  0, .6), ContainsSubstring("probabilities exceed that which induce maximal mixing") );
            
            // must satisfy px, py, pz < 1 - px - py - pz
            REQUIRE_THROWS_WITH( mixPauli(qureg, target,  .3,  .3, .3), ContainsSubstring("probabilities exceed that which induce maximal mixing") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_NOTHROW( mixPauli(vec, 0, 0, 0, 0) ); // zero-prob ok in v4
            REQUIRE_THROWS_WITH( mixPauli(vec, 0, 0.1, 0, 0), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            REQUIRE_THROWS_WITH( mixPauli(vec, 0, 0, 0.1, 0), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            REQUIRE_THROWS_WITH( mixPauli(vec, 0, 0, 0, 0.1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixKrausMap
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        int numOps = GENERATE( range(1,5) ); // max 4 inclusive
        vector<QMatrix> matrs = getRandomKrausMap(1, numOps);
        
        vector<ComplexMatrix2> ops(numOps);
        for (int i=0; i<numOps; i++)
            ops[i] = toComplexMatrix2(matrs[i]);

        v3_mixKrausMap(qureg, target, ops.data(), numOps);
        
        // set ref -> K_i ref K_i^dagger
        vector<QMatrix> matrRefs(numOps);
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

        // cannot use NULL because v3 deprecation API copies ops arg before v4 validation
        ComplexMatrix2 spoofOps[1];
        
        SECTION( "number of operators" ) {
            
            int numOps = 0;

            REQUIRE_THROWS_WITH( v3_mixKrausMap(qureg, 0, spoofOps, numOps), ContainsSubstring("must be given a strictly positive number of matrices") );
        }
        SECTION( "trace preserving" ) {
            
            // valid Kraus map
            int numOps = GENERATE( range(1,5) ); // max 4 inclusive
            vector<QMatrix> matrs = getRandomKrausMap(1, numOps);
            vector<ComplexMatrix2> ops(numOps);
            for (int i=0; i<numOps; i++)
                ops[i] = toComplexMatrix2(matrs[i]);
                
            // make invalid
            ops[GENERATE_REF( range(0,numOps) )].real[0][0] = 9999;
            REQUIRE_THROWS_WITH( v3_mixKrausMap(qureg, 0, ops.data(), numOps), ContainsSubstring("trace preserving") );
        }
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( v3_mixKrausMap(qureg, target, spoofOps, 1), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( v3_mixKrausMap(vec, 0, spoofOps, 1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
        SECTION( "operators fit in node" ) {
            
            qureg.isDistributed = 1;
            qureg.numAmpsPerNode = 3; // min 4
            qureg.logNumAmpsPerNode = 1; // min 2

            REQUIRE_THROWS_WITH( v3_mixKrausMap(qureg, 0, spoofOps, 1), ContainsSubstring("each node's communication buffer") && ContainsSubstring("cannot simultaneously store") );
        }        
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixNonTPKrausMap
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixNonTPKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        int numOps = GENERATE( range(1,5) ); // max 4 inclusive
        
        // map consists of unconstrained 2x2 random matrices
        vector<QMatrix> matrs;
        for (int i=0; i<numOps; i++)
            matrs.push_back(getRandomQMatrix(2));
        
        vector<ComplexMatrix2> ops(numOps);
        for (int i=0; i<numOps; i++)
            ops[i] = toComplexMatrix2(matrs[i]);
        mixNonTPKrausMap(qureg, target, ops.data(), numOps);
        
        // set ref -> K_i ref K_i^dagger
        vector<QMatrix> matrRefs(numOps);
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], target, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref, 1E2*REAL_EPS) );
    }
    SECTION( "input validation" ) {

        // v3 deprecation API copies ops list before v4 validation, so we cannot pass ops=NULL
        ComplexMatrix2 spoofOps[1];
        
        SECTION( "number of operators" ) {
            
            int numOps = 0; 
            REQUIRE_THROWS_WITH( mixNonTPKrausMap(qureg, 0, spoofOps, numOps), ContainsSubstring("must be given a strictly positive number of matrices") );
        }
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixNonTPKrausMap(qureg, target, spoofOps, 1), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( mixNonTPKrausMap(vec, 0, spoofOps, 1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
        SECTION( "operators fit in node" ) {
            
            qureg.isDistributed = 1;
            qureg.numAmpsPerNode = 3; // min 4
            qureg.logNumAmpsPerNode = 1; // min 2
            REQUIRE_THROWS_WITH( mixNonTPKrausMap(qureg, 0, spoofOps, 1), ContainsSubstring("each node's communication buffer") && ContainsSubstring("cannot simultaneously store") );
        }        
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixNonTPMultiQubitKrausMap
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixNonTPMultiQubitKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    // figure out max-num (inclusive) targs allowed by hardware backend
    // (each node must contain as 2^(2*numTargs) amps)
    int maxNumTargs = calcLog2(qureg.numAmpsPerNode) / 2;
    
    SECTION( "correctness" ) {
        
        /* note that this function incurs a stack overhead when numTargs < 4,
         * and a heap overhead when numTargs >= 4
         */
         
        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) ); // inclusive upper bound
        
        // previously, we tried every unique set of targets, via:
        //      int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        // alas, this is too slow for CI, so we instead try a fixed number of random sets:
        GENERATE( range(0, 10) );
        vector<int> targs(numTargs);
        setRandomTargets(targs, NUM_QUBITS);
        
        // try the min and max number of operators, and 2 random numbers 
        // (there are way too many to try all!)
        int maxNumOps = (2*numTargs)*(2*numTargs);
        int numOps = GENERATE_COPY( 1, maxNumOps, take(2,random(1,maxNumOps)) );
        
        // use a new random map, of unconstrained random (2^numTargs x 2^numTargs) matrices
        vector<QMatrix> matrs;
        for (int i=0; i<numOps; i++)
            matrs.push_back(getRandomQMatrix(1 << numTargs));
                
        // create map in QuEST datatypes
        vector<ComplexMatrixN> ops(numOps);
        for (int i=0; i<numOps; i++) {
            ops[i] = createComplexMatrixN(numTargs);
            toComplexMatrixN(matrs[i], ops[i]);
        }
                
        mixNonTPMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps);
                
        // set ref -> K_i ref K_i^dagger
        vector<QMatrix> matrRefs(numOps);
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], targs.data(), numTargs, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref, 1E5*REAL_EPS) );
        
        // cleanup QuEST datatypes
        for (int i=0; i<numOps; i++)
            destroyComplexMatrixN(ops[i]);
    }
    SECTION( "input validation" ) {

        // spoof a list of properly-initialised CompMatr to avoid seg-fault in deprecation API.
        // note some functions will need smaller matrices which is fine; a subset of each
        // spoof'd matrix will be copied over by the deprecation layer
        ComplexMatrixN spoofOps[NUM_QUBITS];
        for (int i=0; i<NUM_QUBITS; i++) {
            spoofOps[i] = createComplexMatrixN(NUM_QUBITS);
            syncCompMatr(spoofOps[i]);
        }
        
        SECTION( "repetition of target" ) {
            
            // make valid targets
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
                
            // duplicate one
            int badInd = GENERATE( range(0,NUM_QUBITS) );
            int copyInd = GENERATE_COPY( filter([=](int i){ return i!=badInd; }, range(0,NUM_QUBITS)) );
            targs[badInd] = targs[copyInd];
            
            REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(qureg, targs, NUM_QUBITS, spoofOps, 1), ContainsSubstring("target qubits") && ContainsSubstring("unique") );
        }
        SECTION( "qubit indices" ) { 
            
            // make valid targets
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
            
            // make one invalid 
            targs[GENERATE( range(0,NUM_QUBITS) )] = GENERATE( -1, NUM_QUBITS );
            
            REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(qureg, targs, NUM_QUBITS, spoofOps, 1), ContainsSubstring("Invalid target qubit") );
        }

        // cannot test this with deprecated API, since 'numOps' informs a copy before validation

            // SECTION( "number of operators" ) {
                
            //     int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
            //     int numOps = GENERATE_REF( -1, 0 ); 
                            
            //     // make valid targets to avoid triggering target validation
            //     vector<int> targs(numTargs);
            //     for (int i=0; i<numTargs; i++)
            //         targs[i] = i;
            //     REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(qureg, targs.data(), numTargs, spoofOps, numOps), ContainsSubstring("must be given a strictly positive number of matrices") );
            // }

        // this validation cannot be performed by the deprecation API which copies
        // over all ops arguments without prior checking them as NULL

            // SECTION( "initialisation of operators" ) {
                
            //     /* compilers don't auto-initialise to NULL; the below circumstance 
            //      * only really occurs when 'malloc' returns NULL in createComplexMatrixN, 
            //      * which actually triggers its own validation. Hence this test is useless 
            //      * currently.
            //      */
                
            //     int numTargs = NUM_QUBITS;
            //     int numOps = (2*numTargs)*(2*numTargs);
                
            //     // no need to initialise ops, but set their attribs correct to avoid triggering other validation
            //     vector<ComplexMatrixN> ops(numOps);
            //     for (int i=0; i<numOps; i++)
            //         ops[i].numQubits = numTargs;
                
            //     // make one of the max-ops explicitly NULL
            //     ops[GENERATE_COPY( range(0,numTargs) )].cpuElems = NULL;
                
            //     // make valid targets to avoid triggering target validation
            //     vector<int> targs(numTargs);
            //     for (int i=0; i<numTargs; i++)
            //         targs[i] = i;
                
            //     REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps), ContainsSubstring("ComplexMatrixN") && ContainsSubstring("created") );
            // }

        // this validation cannot be performed by the deprecation API which copies
        // over all ops arguments without prior checking them as matching in dimension
            
            // SECTION( "dimension of operators" ) {
            
            //     // make valid (dimension-wise) max-qubits Kraus map
            //     int numTargs = NUM_QUBITS;
            //     int numOps = (2*numTargs)*(2*numTargs);
            //     vector<ComplexMatrixN> ops(numOps);
            //     for (int i=0; i<numOps; i++)
            //         ops[i] = createComplexMatrixN(numTargs);
                
            //     // make one have wrong-dimensions 
            //     int badInd = GENERATE_COPY( range(0,numTargs) );
            //     destroyComplexMatrixN(ops[badInd]);
            //     ops[badInd] = createComplexMatrixN(numTargs - 1);
                
            //     // make valid targets to avoid triggering target validation
            //     vector<int> targs(numTargs);
            //     for (int i=0; i<numTargs; i++)
            //         targs[i] = i;
                    
            //     REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(qureg, targs.data(), numTargs, ops.data(), numOps), ContainsSubstring("same number of qubits") );
                
            //     for (int i=0; i<numOps; i++)
            //         destroyComplexMatrixN(ops[i]);
            // }

        SECTION( "density-matrix" ) {
            
            Qureg statevec = createQureg(NUM_QUBITS);
            
            // make valid targets to avoid triggering target validation
            int targs[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                targs[i] = i;
                
            REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(statevec, targs, NUM_QUBITS, spoofOps, 1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(statevec, getQuESTEnv());
            
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
                ops[0].cpuElems[i][i] = 1;
            syncCompMatr(ops[0]);
            
            // fake a smaller qureg 
            qureg.isDistributed = 1;
            qureg.numAmpsPerNode = minAmps - 1;
            qureg.logNumAmpsPerNode = log2(minAmps) - 1;
            REQUIRE_THROWS_WITH( mixNonTPMultiQubitKrausMap(qureg, targs, NUM_QUBITS, ops, 1), ContainsSubstring("each node's communication buffer") && ContainsSubstring("cannot simultaneously store") );
            
            destroyComplexMatrixN(ops[0]);
        }

        for (int i=0; i<NUM_QUBITS; i++)
            destroyComplexMatrixN(spoofOps[i]);
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixNonTPTwoQubitKrausMap
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixNonTPTwoQubitKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int targ1 = GENERATE( range(0,NUM_QUBITS) );
        int targ2 = GENERATE_COPY( filter([=](int t){ return t!=targ1; }, range(0,NUM_QUBITS)) );
        int numOps = GENERATE( range(1,17) ); // max 16 inclusive
        
        // map consists of unconstrained 4x4 random matrices
        vector<QMatrix> matrs;
        for (int i=0; i<numOps; i++)
            matrs.push_back(getRandomQMatrix(4));
        
        vector<ComplexMatrix4> ops(numOps);
        for (int i=0; i<numOps; i++)
            ops[i] = toComplexMatrix4(matrs[i]);
        mixNonTPTwoQubitKrausMap(qureg, targ1, targ2, ops.data(), numOps);
        
        // set ref -> K_i ref K_i^dagger
        int targs[2] = {targ1, targ2};
        vector<QMatrix> matrRefs(numOps);
        for (int i=0; i<numOps; i++) {
            matrRefs[i] = ref;
            applyReferenceOp(matrRefs[i], targs, 2, matrs[i]);
        }
        ref = getZeroMatrix(ref.size());
        for (int i=0; i<numOps; i++)
            ref += matrRefs[i];
        
        REQUIRE( areEqual(qureg, ref, 1E4*REAL_EPS) );
    }
    SECTION( "input validation" ) {

        // deprecated v3 API copies ops before v4 validation, so we cannot pass ops=NULL
        ComplexMatrix4 spoofOps[1];

        // cannot test this with deprecated API, since 'numOps' informs a copy before validation
        
            // SECTION( "number of operators" ) {
                
            //     int numOps = GENERATE( -1, 0 );
            //     REQUIRE_THROWS_WITH( mixNonTPTwoQubitKrausMap(qureg, 0,1, spoofOps, numOps), ContainsSubstring("must be given a strictly positive number of matrices") );
            // }
            
        SECTION( "target collision" ) {
            
            int target = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixNonTPTwoQubitKrausMap(qureg, target, target, spoofOps, 1), ContainsSubstring("target qubits") && ContainsSubstring("unique") );
        }
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixNonTPTwoQubitKrausMap(qureg, 0,target, spoofOps, 1), ContainsSubstring("Invalid target qubit") );
            REQUIRE_THROWS_WITH( mixNonTPTwoQubitKrausMap(qureg, target,0, spoofOps, 1), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( mixNonTPTwoQubitKrausMap(vec, 0,1, spoofOps, 1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
        SECTION( "operators fit in node" ) {
            
            qureg.isDistributed = 1;
            qureg.numAmpsPerNode = 15; // min 16
            qureg.logNumAmpsPerNode = 3; // min 4
            REQUIRE_THROWS_WITH( mixNonTPTwoQubitKrausMap(qureg, 0,1, spoofOps, 1), ContainsSubstring("each node's communication buffer") && ContainsSubstring("cannot simultaneously store") );
        }        
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixTwoQubitDephasing
 * @ingroup deprecatedtests 
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
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, 0, targ, 0), ContainsSubstring("Invalid target") );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, targ, 0, 0), ContainsSubstring("Invalid target") );
        }
        SECTION( "target collision" ) {
            
            int targ = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, targ, targ, 0), ContainsSubstring("target") && ContainsSubstring("unique") );
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, 0, 1, -.1), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(qureg, 0, 1, 3/4. + .01), ContainsSubstring("probability") && ContainsSubstring("3/4") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_NOTHROW( mixTwoQubitDephasing(vec, 0, 1, 0) ); // zero-prob ok in v4
            REQUIRE_THROWS_WITH( mixTwoQubitDephasing(vec, 0, 1, 0.1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixTwoQubitDepolarising
 * @ingroup deprecatedtests 
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
                QMatrix op = getKroneckerProduct(paulis[i], paulis[j]);
                applyReferenceOp(term, targs, 2, op);
                ref += (prob/15.) * term;
            }
        }
        
        REQUIRE( areEqual(qureg, ref, 1E4*REAL_EPS) );
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int targ = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, 0, targ, 0), ContainsSubstring("Invalid target") );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, targ, 0, 0), ContainsSubstring("Invalid target") );
        }
        SECTION( "target collision" ) {
            
            int targ = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, targ, targ, 0), ContainsSubstring("target") && ContainsSubstring("unique") );
        }
        SECTION( "probability" ) {
            
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, 0, 1, -.1), ContainsSubstring("probability is invalid") );
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(qureg, 0, 1, 15/16. + .01), ContainsSubstring("probability") && ContainsSubstring("15/16") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_NOTHROW( mixTwoQubitDepolarising(vec, 0, 1, 0) ); // zero prob ok in v4
            REQUIRE_THROWS_WITH( mixTwoQubitDepolarising(vec, 0, 1, 0.1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
    }
    destroyQureg(qureg, getQuESTEnv());
}



/** @sa mixTwoQubitKrausMap
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "mixTwoQubitKrausMap", "[decoherence]" ) {
    
    PREPARE_TEST(qureg, ref);
    
    SECTION( "correctness" ) {
        
        int targ1 = GENERATE( range(0,NUM_QUBITS) );
        int targ2 = GENERATE_COPY( filter([=](int t){ return t!=targ1; }, range(0,NUM_QUBITS)) );
        int numOps = GENERATE( range(1,17) ); // max 16 inclusive
        vector<QMatrix> matrs = getRandomKrausMap(2, numOps);
        
        vector<ComplexMatrix4> ops(numOps);
        for (int i=0; i<numOps; i++)
            ops[i] = toComplexMatrix4(matrs[i]);
        mixTwoQubitKrausMap(qureg, targ1, targ2, ops.data(), numOps);
        
        // set ref -> K_i ref K_i^dagger
        int targs[2] = {targ1, targ2};
        vector<QMatrix> matrRefs(numOps);
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

        ComplexMatrix4 emptyOps[1];
        
        SECTION( "number of operators" ) {
            
            // in v4, a separate createKausMap call is made which throws the below error
            int numOps = 0;
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,1, emptyOps, numOps), ContainsSubstring("must be given a strictly positive number of matrices") );
        }
        SECTION( "trace preserving" ) {
            
            // valid Kraus map
            int numOps = GENERATE( range(1,16) );
            vector<QMatrix> matrs = getRandomKrausMap(2, numOps);
            vector<ComplexMatrix4> ops(numOps);
            for (int i=0; i<numOps; i++)
                ops[i] = toComplexMatrix4(matrs[i]);
                
            // make only one of the ops at a time invalid
            ops[GENERATE_REF( range(0,numOps) )].real[0][0] = 999;
            
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,1, ops.data(), numOps), ContainsSubstring("trace preserving") );
        }
        SECTION( "target collision" ) {
            
            int target = GENERATE( range(0,NUM_QUBITS) );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, target, target, emptyOps, 1), ContainsSubstring("target qubits") && ContainsSubstring("unique") );
        }
        SECTION( "qubit index" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,target, emptyOps, 1), ContainsSubstring("Invalid target qubit") );
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, target,0, emptyOps, 1), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg vec = createQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(vec, 0,1, emptyOps, 1), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec, getQuESTEnv());
        }
        SECTION( "operators fit in node" ) {
            
            qureg.isDistributed = 1;
            qureg.numAmpsPerNode = 15; // min 16
            qureg.logNumAmpsPerNode = 3; // min 4
            REQUIRE_THROWS_WITH( mixTwoQubitKrausMap(qureg, 0,1, emptyOps, 1), ContainsSubstring("each node's communication buffer") && ContainsSubstring("cannot simultaneously store") );
        }        
    }
    destroyQureg(qureg, getQuESTEnv());
}
