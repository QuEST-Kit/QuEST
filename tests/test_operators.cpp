
#include "catch.hpp"
#include "QuEST.h"
#include "utilities.hpp"

/** Prepares the needed data structures for unit testing some operators.
 * This creates a statevector and density matrix of the size NUM_QUBITS,
 * and corresponding QVector and QMatrix instances for analytic comparison.
 */
#define PREPARE_TEST(quregVec, quregMatr, refVec, refMatr) \
    Qureg quregVec = createQureg(NUM_QUBITS, QUEST_ENV); \
    Qureg quregMatr = createDensityQureg(NUM_QUBITS, QUEST_ENV); \
    initDebugState(quregVec); \
    initDebugState(quregMatr); \
    QVector refVec = toQVector(quregVec); \
    QMatrix refMatr = toQMatrix(quregMatr);

/** Destroys the data structures made by PREPARE_TEST */
#define CLEANUP_TEST(quregVec, quregMatr) \
    destroyQureg(quregVec, QUEST_ENV); \
    destroyQureg(quregMatr, QUEST_ENV);

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;



/** @sa applyMatrix2
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyMatrix2", "[operators]" ) {
        
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    // every test will use a unique random matrix
    QMatrix op = getRandomQMatrix(2); // 2-by-2
    ComplexMatrix2 matr = toComplexMatrix2(op); 

    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        
        // reference boilerplate
        int* ctrls = NULL;
        int numCtrls = 0;
        int targs[] = {target};
        int numTargs = 1;
        
        SECTION( "state-vector" ) {
        
            applyMatrix2(quregVec, target, matr);
            applyReferenceMatrix(refVec, ctrls, numCtrls, targs, numTargs, op);

            REQUIRE( areEqual(quregVec, refVec) );
        }
        SECTION( "density-matrix" ) {
        
            applyMatrix2(quregMatr, target, matr);
            applyReferenceMatrix(refMatr, ctrls, numCtrls, targs, numTargs, op);
            
            REQUIRE( areEqual(quregMatr, refMatr, 10*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyMatrix2(quregVec, target, matr), Contains("Invalid target") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyMatrix4
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyMatrix4", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    // in distributed mode, each node must be able to fit all amps modified by matrix 
    REQUIRE( quregVec.numAmpsPerChunk >= 4 );
    
    // every test will use a unique random matrix
    QMatrix op = getRandomQMatrix(4); // 4-by-4
    ComplexMatrix4 matr = toComplexMatrix4(op); 

    SECTION( "correctness" ) {
        
        int targ1 = GENERATE( range(0,NUM_QUBITS) );
        int targ2 = GENERATE_COPY( filter([=](int t){ return t!=targ1; }, range(0,NUM_QUBITS)) );
        
        // reference boilerplate
        int* ctrls = NULL;
        int numCtrls = 0;
        int targs[] = {targ1, targ2};
        int numTargs = 2;
        
        SECTION( "state-vector" ) {
        
            applyMatrix4(quregVec, targ1, targ2, matr);
            applyReferenceMatrix(refVec, ctrls, numCtrls, targs, numTargs, op);
            REQUIRE( areEqual(quregVec, refVec) );
        }
        SECTION( "density-matrix" ) {

            applyMatrix4(quregMatr, targ1, targ2, matr);
            applyReferenceMatrix(refMatr, ctrls, numCtrls, targs, numTargs, op);
            REQUIRE( areEqual(quregMatr, refMatr, 10*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int targ1 = GENERATE( -1, NUM_QUBITS );
            int targ2 = 0;
            REQUIRE_THROWS_WITH( applyMatrix4(quregVec, targ1, targ2, matr), Contains("Invalid target") );
            REQUIRE_THROWS_WITH( applyMatrix4(quregVec, targ2, targ1, matr), Contains("Invalid target") );
        }
        SECTION( "repetition of targets" ) {
            
            int qb = 0;
            REQUIRE_THROWS_WITH( applyMatrix4(quregVec, qb, qb, matr), Contains("target") && Contains("unique") );
        }
        SECTION( "matrix fits in node" ) {
                
            // pretend we have a very limited distributed memory
            quregVec.numAmpsPerChunk = 1;
            REQUIRE_THROWS_WITH( applyMatrix4(quregVec, 0, 1, matr), Contains("targets too many qubits"));
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyMatrixN
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyMatrixN", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    // figure out max-num (inclusive) targs allowed by hardware backend
    int maxNumTargs = calcLog2(quregVec.numAmpsPerChunk);
    
    SECTION( "correctness" ) {
        
        // generate all possible qubit arrangements
        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) ); // inclusive upper bound
        int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        
        // for each qubit arrangement, use a new random matrix
        QMatrix op = getRandomQMatrix(1 << numTargs);
        ComplexMatrixN matr = createComplexMatrixN(numTargs);
        toComplexMatrixN(op, matr);
        
        // reference boilerplate
        int* ctrls = NULL;
        int numCtrls = 0;
    
        SECTION( "state-vector" ) {
            
            applyMatrixN(quregVec, targs, numTargs, matr);
            applyReferenceMatrix(refVec, ctrls, numCtrls, targs, numTargs, op);
            REQUIRE( areEqual(quregVec, refVec) );
        }
        SECTION( "density-matrix" ) {

            applyMatrixN(quregMatr, targs, numTargs, matr);
            applyReferenceMatrix(refMatr, ctrls, numCtrls, targs, numTargs, op);
            REQUIRE( areEqual(quregMatr, refMatr, 100*REAL_EPS) );
        }
        destroyComplexMatrixN(matr);
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of targets" ) {
            
            // there cannot be more targets than qubits in register
            int numTargs = GENERATE( -1, 0, NUM_QUBITS+1 );
            int targs[NUM_QUBITS+1]; // prevents seg-fault if validation doesn't trigger
            ComplexMatrixN matr = createComplexMatrixN(NUM_QUBITS+1); // prevent seg-fault
            
            REQUIRE_THROWS_WITH( applyMatrixN(quregVec, targs, numTargs, matr), Contains("Invalid number of target"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "repetition in targets" ) {
            
            int numTargs = 3;
            int targs[] = {1,2,2};
            ComplexMatrixN matr = createComplexMatrixN(numTargs); // prevents seg-fault if validation doesn't trigger
            
            REQUIRE_THROWS_WITH( applyMatrixN(quregVec, targs, numTargs, matr), Contains("target") && Contains("unique"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "qubit indices" ) {
            
            int numTargs = 3;
            int targs[] = {1,2,3};
            ComplexMatrixN matr = createComplexMatrixN(numTargs); // prevents seg-fault if validation doesn't trigger
            
            int inv = GENERATE( -1, NUM_QUBITS );
            targs[GENERATE_COPY( range(0,numTargs) )] = inv; // make invalid target
            REQUIRE_THROWS_WITH( applyMatrixN(quregVec, targs, numTargs, matr), Contains("Invalid target") );
            
            destroyComplexMatrixN(matr);
        }
        SECTION( "matrix creation" ) {
            
            int numTargs = 3;
            int targs[] = {1,2,3};
            
            /* compilers don't auto-initialise to NULL; the below circumstance 
             * only really occurs when 'malloc' returns NULL in createComplexMatrixN, 
             * which actually triggers its own validation. Hence this test is useless 
             * currently.
             */
            ComplexMatrixN matr;
            matr.real = NULL;
            matr.imag = NULL; 
            REQUIRE_THROWS_WITH( applyMatrixN(quregVec, targs, numTargs, matr), Contains("created") );
        }
        SECTION( "matrix dimensions" ) {
            
            int targs[2] = {1,2};
            ComplexMatrixN matr = createComplexMatrixN(3); // intentionally wrong size
            
            REQUIRE_THROWS_WITH( applyMatrixN(quregVec, targs, 2, matr), Contains("matrix size"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "matrix fits in node" ) {
                
            // pretend we have a very limited distributed memory (judged by matr size)
            quregVec.numAmpsPerChunk = 1;
            int qb[] = {1,2};
            ComplexMatrixN matr = createComplexMatrixN(2); // prevents seg-fault if validation doesn't trigger
            REQUIRE_THROWS_WITH( applyMatrixN(quregVec, qb, 2, matr), Contains("targets too many qubits"));
            destroyComplexMatrixN(matr);
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyMultiControlledMatrixN
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyMultiControlledMatrixN", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    // figure out max-num targs (inclusive) allowed by hardware backend
    int maxNumTargs = calcLog2(quregVec.numAmpsPerChunk);
    if (maxNumTargs >= NUM_QUBITS)
        maxNumTargs = NUM_QUBITS - 1; // leave room for min-number of control qubits
        
    SECTION( "correctness" ) {
        
        // try all possible numbers of targets and controls
        int numTargs = GENERATE_COPY( range(1,maxNumTargs+1) );
        int maxNumCtrls = NUM_QUBITS - numTargs;
        int numCtrls = GENERATE_COPY( range(1,maxNumCtrls+1) );
        
        // generate all possible valid qubit arrangements
        int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        int* ctrls = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numCtrls, targs, numTargs) );
        
        // for each qubit arrangement, use a new random unitary
        QMatrix op = getRandomQMatrix(1 << numTargs);
        ComplexMatrixN matr = createComplexMatrixN(numTargs);
        toComplexMatrixN(op, matr);
    
        SECTION( "state-vector" ) {
            
            applyMultiControlledMatrixN(quregVec, ctrls, numCtrls, targs, numTargs, matr);
            applyReferenceMatrix(refVec, ctrls, numCtrls, targs, numTargs, op);
            REQUIRE( areEqual(quregVec, refVec) );
        }
        SECTION( "density-matrix" ) {

            applyMultiControlledMatrixN(quregMatr, ctrls, numCtrls, targs, numTargs, matr);
            applyReferenceMatrix(refMatr, ctrls, numCtrls, targs, numTargs, op);
            REQUIRE( areEqual(quregMatr, refMatr, 100*REAL_EPS) );
        }
        destroyComplexMatrixN(matr);
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of targets" ) {
            
            // there cannot be more targets than qubits in register
            // (numTargs=NUM_QUBITS is caught elsewhere, because that implies ctrls are invalid)
            int numTargs = GENERATE( -1, 0, NUM_QUBITS+1 );
            int targs[NUM_QUBITS+1]; // prevents seg-fault if validation doesn't trigger
            int ctrls[] = {0};
            ComplexMatrixN matr = createComplexMatrixN(NUM_QUBITS+1); // prevent seg-fault
            toComplexMatrixN(getRandomQMatrix( 1 << (NUM_QUBITS+1)), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 1, targs, numTargs, matr), Contains("Invalid number of target"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "repetition in targets" ) {
            
            int ctrls[] = {0};
            int numTargs = 3;
            int targs[] = {1,2,2};
            ComplexMatrixN matr = createComplexMatrixN(numTargs); // prevents seg-fault if validation doesn't trigger
            toComplexMatrixN(getRandomQMatrix(1 << numTargs), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 1, targs, numTargs, matr), Contains("target") && Contains("unique"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "number of controls" ) {
            
            int numCtrls = GENERATE( -1, 0, NUM_QUBITS, NUM_QUBITS+1 );
            int ctrls[NUM_QUBITS+1]; // avoids seg-fault if validation not triggered
            int targs[1] = {0};
            ComplexMatrixN matr = createComplexMatrixN(1);
            toComplexMatrixN(getRandomQMatrix(1 << 1), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, numCtrls, targs, 1, matr), Contains("Invalid number of control"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "repetition in controls" ) {
            
            int ctrls[] = {0,1,1};
            int targs[] = {3};
            ComplexMatrixN matr = createComplexMatrixN(1);
            toComplexMatrixN(getRandomQMatrix(1 << 1), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 3, targs, 1, matr), Contains("control") && Contains("unique"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "control and target collision" ) {
            
            int ctrls[] = {0,1,2};
            int targs[] = {3,1,4};
            ComplexMatrixN matr = createComplexMatrixN(3);
            toComplexMatrixN(getRandomQMatrix(1 << 3), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 3, targs, 3, matr), Contains("Control") && Contains("target") && Contains("disjoint"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "qubit indices" ) {
            
            // valid inds
            int numQb = 2;
            int qb1[2] = {0,1};
            int qb2[2] = {2,3};
            ComplexMatrixN matr = createComplexMatrixN(numQb);
            toComplexMatrixN(getRandomQMatrix(1 << numQb), matr);
            
            // make qb1 invalid
            int inv = GENERATE( -1, NUM_QUBITS );
            qb1[GENERATE_COPY(range(0,numQb))] = inv;
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, qb1, numQb, qb2, numQb, matr), Contains("Invalid control") );
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, qb2, numQb, qb1, numQb, matr), Contains("Invalid target") );
            destroyComplexMatrixN(matr);
        }
        SECTION( "matrix creation" ) {
            
            int ctrls[1] = {0};
            int targs[3] = {1,2,3};
            
            /* compilers don't auto-initialise to NULL; the below circumstance 
             * only really occurs when 'malloc' returns NULL in createComplexMatrixN, 
             * which actually triggers its own validation. Hence this test is useless 
             * currently.
             */
            ComplexMatrixN matr;
            matr.real = NULL;
            matr.imag = NULL; 
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 1, targs, 3, matr), Contains("created") );
        }
        SECTION( "matrix dimensions" ) {
            
            int ctrls[1] = {0};
            int targs[2] = {1,2};
            ComplexMatrixN matr = createComplexMatrixN(3); // intentionally wrong size
            toComplexMatrixN(getRandomQMatrix(1 << 3), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 1, targs, 2, matr), Contains("matrix size"));
            destroyComplexMatrixN(matr);
        }
        SECTION( "matrix fits in node" ) {
                
            // pretend we have a very limited distributed memory (judged by matr size)
            quregVec.numAmpsPerChunk = 1;
            int ctrls[1] = {0};
            int targs[2] = {1,2};
            ComplexMatrixN matr = createComplexMatrixN(2);
            toComplexMatrixN(getRandomQMatrix(1 << 2), matr);
            
            REQUIRE_THROWS_WITH( applyMultiControlledMatrixN(quregVec, ctrls, 1, targs, 2, matr), Contains("targets too many qubits"));
            destroyComplexMatrixN(matr);
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyPauliHamil
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyPauliHamil", "[operators]" ) {
    
    Qureg vecIn = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg vecOut = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg matIn = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    Qureg matOut = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
    initDebugState(vecIn);
    initDebugState(matIn);

    SECTION( "correctness" ) {
        
        /* it's too expensive to try every possible Pauli configuration, so
         * we'll try 10 random codes, and for each, random coefficients
         */
        GENERATE( range(0,10) );
        
        int numTerms = GENERATE( 1, 2, 10, 15 );
        PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);
        setRandomPauliSum(hamil);
        QMatrix refHamil = toQMatrix(hamil);
        
        SECTION( "state-vector" ) {
            
            QVector vecRef = toQVector(vecIn);
            applyPauliHamil(vecIn, hamil, vecOut);
            
            // ensure vecIn barely changes under precision
            REQUIRE( areEqual(vecIn, vecRef) );
            
            // ensure vecOut changed correctly 
            REQUIRE( areEqual(vecOut, refHamil * vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            QMatrix matRef = toQMatrix(matIn);
            applyPauliHamil(matIn, hamil, matOut);
            
            // ensure matIn barely changes under precision
            REQUIRE( areEqual(matIn, matRef) );
            
            // ensure matOut changed correctly 
            REQUIRE( areEqual(matOut, refHamil * matRef, 1E2*REAL_EPS) );
        }
        
        destroyPauliHamil(hamil);
    }
    SECTION( "input validation" ) {
        
        SECTION( "pauli codes" ) {
            
            int numTerms = 3;
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);

            // make one pauli code wrong
            hamil.pauliCodes[GENERATE_COPY( range(0,numTerms*NUM_QUBITS) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( applyPauliHamil(vecIn, hamil, vecOut), Contains("Invalid Pauli code") );
            
            destroyPauliHamil(hamil);
        }
        SECTION( "qureg dimensions" ) {
            
            Qureg badVec = createQureg(NUM_QUBITS+1, QUEST_ENV);
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 1);

            REQUIRE_THROWS_WITH( applyPauliHamil(vecIn, hamil, badVec), Contains("Dimensions of the qubit registers don't match") );
            
            destroyQureg(badVec, QUEST_ENV);
            destroyPauliHamil(hamil);
        }
        SECTION( "qureg types" ) {
            
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 1);
            
            REQUIRE_THROWS_WITH( applyPauliHamil(vecIn, hamil, matOut), Contains("Registers must both be state-vectors or both be density matrices") );
            
            destroyPauliHamil(hamil);
        }
        SECTION( "matching hamiltonian qubits" ) {
            
            PauliHamil hamil = createPauliHamil(NUM_QUBITS + 1, 1);
            
            REQUIRE_THROWS_WITH( applyPauliHamil(vecIn, hamil, vecOut), Contains("same number of qubits") );
            REQUIRE_THROWS_WITH( applyPauliHamil(matIn, hamil, matOut), Contains("same number of qubits") );
            
            destroyPauliHamil(hamil);
        }
    }
    destroyQureg(vecIn, QUEST_ENV);
    destroyQureg(vecOut, QUEST_ENV);
    destroyQureg(matIn, QUEST_ENV);
    destroyQureg(matOut, QUEST_ENV);
}



/** @sa applyPauliSum
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyPauliSum", "[operators]" ) {
    
    Qureg vecIn = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg vecOut = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg matIn = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    Qureg matOut = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
    initDebugState(vecIn);
    initDebugState(matIn);

    SECTION( "correctness" ) {
        
        /* it's too expensive to try ALL Pauli sequences, via 
         *      pauliOpType* paulis = GENERATE_COPY( pauliseqs(numPaulis) );.
         * Furthermore, take(10, pauliseqs(numTargs)) will try the same pauli codes.
         * Hence, we instead opt to repeatedly randomly generate pauliseqs
         */
        GENERATE( range(0,10) ); // gen 10 random pauli-codes
        
        int numTerms = GENERATE( 1, 2, 10, 15);
        int numPaulis = numTerms * NUM_QUBITS;
        
        // each test will use random coefficients
        qreal coeffs[numTerms];
        pauliOpType paulis[numPaulis];
        setRandomPauliSum(coeffs, paulis, NUM_QUBITS, numTerms);
        QMatrix pauliSum = toQMatrix(coeffs, paulis, NUM_QUBITS, numTerms);
        
        SECTION( "state-vector" ) {
            
            QVector vecRef = toQVector(vecIn);
            applyPauliSum(vecIn, paulis, coeffs, numTerms, vecOut);
            
            // ensure vecIn barely changes under precision
            REQUIRE( areEqual(vecIn, vecRef) );
            
            // ensure vecOut changed correctly 
            REQUIRE( areEqual(vecOut, pauliSum * vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            QMatrix matRef = toQMatrix(matIn);
            applyPauliSum(matIn, paulis, coeffs, numTerms, matOut);
            
            // ensure matIn barely changes under precision
            REQUIRE( areEqual(matIn, matRef) );
            
            // ensure matOut changed correctly 
            REQUIRE( areEqual(matOut, pauliSum * matRef, 1E2*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of terms" ) {
            
            int numTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( applyPauliSum(vecIn, NULL, NULL, numTerms, vecOut), Contains("Invalid number of terms") );
        }
        SECTION( "pauli codes" ) {
            
            // valid codes
            int numTerms = 3;
            int numPaulis = numTerms*NUM_QUBITS;
            pauliOpType paulis[numPaulis];
            for (int i=0; i<numPaulis; i++)
                paulis[i] = PAULI_I;
            
            // make one invalid 
            paulis[GENERATE_COPY( range(0,numPaulis) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( applyPauliSum(vecIn, paulis, NULL, numTerms, vecOut), Contains("Invalid Pauli code") );
        }
        SECTION( "qureg dimensions" ) {
            
            Qureg badVec = createQureg(NUM_QUBITS+1, QUEST_ENV);
            pauliOpType paulis[NUM_QUBITS];
            REQUIRE_THROWS_WITH( applyPauliSum(vecIn, paulis, NULL, 1, badVec), Contains("Dimensions of the qubit registers don't match") );
            destroyQureg(badVec, QUEST_ENV);
        }
        SECTION( "qureg types" ) {
            
            pauliOpType paulis[NUM_QUBITS];
            REQUIRE_THROWS_WITH( applyPauliSum(vecIn, paulis, NULL, 1, matOut), Contains("Registers must both be state-vectors or both be density matrices") );
        }
    }
    destroyQureg(vecIn, QUEST_ENV);
    destroyQureg(vecOut, QUEST_ENV);
    destroyQureg(matIn, QUEST_ENV);
    destroyQureg(matOut, QUEST_ENV);
}



/** @sa applyTrotterCircuit
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyTrotterCircuit", "[operators]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    initDebugState(vec);
    initDebugState(mat);
    
    Qureg vecRef = createCloneQureg(vec, QUEST_ENV);
    Qureg matRef = createCloneQureg(mat, QUEST_ENV);

    SECTION( "correctness" ) {
    
        SECTION( "one term" ) {
            
            // a Hamiltonian with one term has an exact (trivial) Trotterisation
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 1);
            
            // H = coeff X Y Z (on qubits 0,1,2)
            qreal coeff = getRandomReal(-5, 5);
            int numTargs = 3;
            int targs[] =  {0, 1, 2};
            pauliOpType codes[] = {PAULI_X, PAULI_Y, PAULI_Z};
            hamil.termCoeffs[0] = coeff;
            for (int i=0; i<numTargs; i++)
                hamil.pauliCodes[targs[i]] = codes[i];
                
            // time can be negative
            qreal time = getRandomReal(-2,2);
            
            // by commutation, all reps & orders yield the same total unitary
            int reps = GENERATE( range(1,5) );
            
            // applyTrotter(t, 1, 1) = exp(-i t coeff paulis)
            // multiRotatePauli(a) = exp(- i a / 2 paulis)

            SECTION( "state-vector" ) {
                
                int order = GENERATE( 1, 2, 4 );
                
                applyTrotterCircuit(vec, hamil, time, order, reps);
                multiRotatePauli(vecRef, targs, codes, numTargs, 2*time*coeff);
                REQUIRE( areEqual(vec, vecRef) );
            }
            SECTION( "density-matrix" ) {
                
                int order = GENERATE( 1, 2 ); // precision bites density-matrices quickly
                
                applyTrotterCircuit(mat, hamil, time, order, reps);
                multiRotatePauli(matRef, targs, codes, numTargs, 2*time*coeff);
                REQUIRE( areEqual(mat, matRef, 1E2*REAL_EPS) );
            }
                    
            destroyPauliHamil(hamil);
        }
        SECTION( "commuting terms" ) {
            
            // a Hamiltonian of commuting terms, Trotterises exactly
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 2);
            
            // H = c0 X Y I + c1 X I Z (on qubits 0,1,2)
            int targs[] = {0, 1, 2};
            hamil.pauliCodes[0] = PAULI_X;
            hamil.pauliCodes[0 + NUM_QUBITS] = PAULI_X;
            hamil.pauliCodes[1] = PAULI_Y;
            hamil.pauliCodes[1 + NUM_QUBITS] = PAULI_I;
            hamil.pauliCodes[2] = PAULI_I;
            hamil.pauliCodes[2 + NUM_QUBITS] = PAULI_Z;
            for (int i=0; i<hamil.numSumTerms; i++)
                hamil.termCoeffs[i] = getRandomReal(-5,5); 
            
            // time can be negative
            qreal time = getRandomReal(-2,2);
            
            // applyTrotter(t, 1, 1) = exp(-i t c0 paulis0) exp(-i t c1 paulis1)
            // multiRotatePauli(a) = exp(- i a / 2 paulis)
            
            SECTION( "state-vector" ) {
                
                int reps = GENERATE( range(1,5) );
                int order = GENERATE( 1, 2, 4 );
                
                applyTrotterCircuit(vec, hamil, time, order, reps);
                multiRotatePauli(vecRef, targs, hamil.pauliCodes, 3, 2*time*hamil.termCoeffs[0]);
                multiRotatePauli(vecRef, targs, &(hamil.pauliCodes[NUM_QUBITS]), 3, 2*time*hamil.termCoeffs[1]);
                REQUIRE( areEqual(vec, vecRef, 10*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                int reps = GENERATE( range(1,5) );
                int order = GENERATE( 1, 2 ); // precision hurts density matrices quickly
                
                applyTrotterCircuit(mat, hamil, time, order, reps);
                multiRotatePauli(matRef, targs, hamil.pauliCodes, 3, 2*time*hamil.termCoeffs[0]);
                multiRotatePauli(matRef, targs, &(hamil.pauliCodes[NUM_QUBITS]), 3, 2*time*hamil.termCoeffs[1]);
                REQUIRE( areEqual(mat, matRef, 1E2*REAL_EPS) );
            }

            destroyPauliHamil(hamil);
        }
        SECTION( "general" ) {
            
            /* We'll consider an analytic time-evolved state, so that we can avoid 
             * comparing applyTrotterCircuit to other numerical approximations.
             * We can construct such a state, by using a Hamiltonian with known
             * analytic eigenvalues, and hence a known period. Time evolution of the 
             * period will just yield the input state.
             *
             * E.g. H = pi sqrt(2) X Y Z X Y + pi Y Z X Y Z + pi Z X Y Z X
             * has (degenerate) eigenvalues +- 2 pi, so the period
             * of the Hamiltonian is t=1. 
             */
             
            // hardcoded 5 qubits here in the Pauli codes
            REQUIRE( NUM_QUBITS == 5 );
            
            // (verbose for C++ compatibility)
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 3);
            qreal coeffs[] = {M_PI * sqrt(2.0), M_PI, M_PI};
            enum pauliOpType codes[] = {
                PAULI_X, PAULI_Y, PAULI_Z, PAULI_X, PAULI_Y,
                PAULI_Y, PAULI_Z, PAULI_X, PAULI_Y, PAULI_Z,
                PAULI_Z, PAULI_X, PAULI_Y, PAULI_Z, PAULI_X};
            initPauliHamil(hamil, coeffs, codes);
                    
            // evolving to t=1 should leave the input state unchanged
            qreal time = 1;

            // since unnormalised (initDebugState), max fid is 728359.8336
            qreal fidNorm = 728359.8336;
            
            SECTION( "absolute" ) {
                
                // such a high order and reps should yield precise solution
                int order = 4;
                int reps = 20;
                applyTrotterCircuit(vec, hamil, time, 4, 20);
                qreal fid = calcFidelity(vec, vecRef) / fidNorm;
                
                REQUIRE( fid == Approx(1).epsilon(1E-8) );
            }
            SECTION( "repetitions scaling" ) {
                
                // exclude order 1; too few reps for monotonic increase of accuracy
                int order = GENERATE( 2, 4, 6 ); 
                
                // accuracy should increase with increasing repetitions
                int reps[] = {1, 5, 10};
                
                qreal prevFid = 0;
                for (int i=0; i<3; i++) {
                    initDebugState(vec);
                    applyTrotterCircuit(vec, hamil, time, order, reps[i]);
                    qreal fid = calcFidelity(vec, vecRef) / fidNorm;

                    REQUIRE( fid >= prevFid );
                    prevFid = fid;
                }
            }
            SECTION( "order scaling" ) {
                
                // exclude order 1; too few reps for monotonic increase of accuracy
                int reps = GENERATE( 5, 10 ); 
                
                // accuracy should increase with increasing repetitions
                int orders[] = {1, 2, 4, 6};
                
                qreal prevFid = 0;
                for (int i=0; i<4; i++) {
                    initDebugState(vec);
                    applyTrotterCircuit(vec, hamil, time, orders[i], reps);
                    qreal fid = calcFidelity(vec, vecRef) / fidNorm;

                    REQUIRE( fid >= prevFid );
                    prevFid = fid;
                }
            }
            
            destroyPauliHamil(hamil);
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "repetitions" ) {
            
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 1);
            int reps = GENERATE( -1, 0 );
            
            REQUIRE_THROWS_WITH( applyTrotterCircuit(vec, hamil, 1, 1, reps), Contains("repetitions must be >=1") );
            
            destroyPauliHamil(hamil);
        }
        SECTION( "order" ) {
            
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 1);
            int order = GENERATE( -1, 0, 3, 5, 7 );
            
            REQUIRE_THROWS_WITH( applyTrotterCircuit(vec, hamil, 1, order, 1), Contains("order must be 1, or an even number") );
            
            destroyPauliHamil(hamil);
        }
        SECTION( "pauli codes" ) {
            
            int numTerms = 3;
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);

            // make one pauli code wrong
            hamil.pauliCodes[GENERATE_COPY( range(0,numTerms*NUM_QUBITS) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( applyTrotterCircuit(vec, hamil, 1, 1, 1), Contains("Invalid Pauli code") );
            
            destroyPauliHamil(hamil);
        }
        SECTION( "matching hamiltonian qubits" ) {
            
            PauliHamil hamil = createPauliHamil(NUM_QUBITS + 1, 1);
            
            REQUIRE_THROWS_WITH( applyTrotterCircuit(vec, hamil, 1, 1, 1), Contains("same number of qubits") );
            REQUIRE_THROWS_WITH( applyTrotterCircuit(mat, hamil, 1, 1, 1), Contains("same number of qubits") );
            
            destroyPauliHamil(hamil);
        }
    }
    
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
    destroyQureg(vecRef, QUEST_ENV);
    destroyQureg(matRef, QUEST_ENV);
}

