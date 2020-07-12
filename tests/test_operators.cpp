
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
        
        int numTerms = GENERATE( 1, 2, 10, 15);
        int numPaulis = numTerms * NUM_QUBITS;
        
        // each test will use random coefficients
        qreal coeffs[numTerms];
        for (int i=0; i<numTerms; i++)
            coeffs[i] = getRandomReal(-5, 5);
            
        /* it's too expensive to try ALL Pauli sequences, via 
         *      pauliOpType* paulis = GENERATE_COPY( pauliseqs(numPaulis) );.
         * Furthermore, take(10, pauliseqs(numTargs)) will try the same pauli codes.
         * Hence, we instead opt to repeatedly randomly generate pauliseqs
         */
        GENERATE( range(0,10) ); // gen 10 random pauli-codes
        pauliOpType paulis[numPaulis];
        for (int i=0; i<numPaulis; i++)
            paulis[i] = (pauliOpType) getRandomInt(0,4);
            
        // build the resulting Pauli sum; sum_j coeffs[j] * paulis[n*j : (n+1)j]
        QMatrix iMatr{{1,0},{0,1}};
        QMatrix xMatr{{0,1},{1,0}};
        QMatrix yMatr{{0,-1i},{1i,0}};
        QMatrix zMatr{{1,0},{0,-1}};
        QMatrix pauliSum = getZeroMatrix(1<<NUM_QUBITS);
        
        for (int t=0; t<numTerms; t++) {
            QMatrix pauliProd = QMatrix{{1}};
            
            for (int q=0; q<NUM_QUBITS; q++) {
                int i = q + t*NUM_QUBITS;
                
                QMatrix fac;
                if (paulis[i] == PAULI_I) fac = iMatr;
                if (paulis[i] == PAULI_X) fac = xMatr;
                if (paulis[i] == PAULI_Y) fac = yMatr;
                if (paulis[i] == PAULI_Z) fac = zMatr;
                pauliProd = getKroneckerProduct(fac, pauliProd);
            }
            pauliSum += coeffs[t] * pauliProd;
        }
        
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