
#include "catch.hpp"
#include "QuEST.h"
#include "QuEST_test_utils.hpp"

/** The default number of qubits in the registers created for unit testing 
 * (both statevectors and density matrices). Creation of non-NUM_QUBITS sized 
 * Quregs should be justified in a comment. 
 * Note that the smaller this number is, the fewer nodes can be employed in 
 * distribution testing, since each node must contain at least one amplitude.
 * Furthermore, the larger this number is, the greater the deviation of correct 
 * results from their expected value, due to numerical error; this is especially 
 * apparent for density matrices.
 */
#define NUM_QUBITS 5

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;



TEST_CASE( "applyPauliSum", "[operators]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vecIn = createQureg(NUM_QUBITS, env);
    Qureg vecOut = createQureg(NUM_QUBITS, env);
    Qureg matIn = createDensityQureg(NUM_QUBITS, env);
    Qureg matOut = createDensityQureg(NUM_QUBITS, env);
    
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
            
            Qureg badVec = createQureg(NUM_QUBITS+1, env);
            pauliOpType paulis[NUM_QUBITS];
            REQUIRE_THROWS_WITH( applyPauliSum(vecIn, paulis, NULL, 1, badVec), Contains("Dimensions of the qubit registers don't match") );
            destroyQureg(badVec, env);
        }
        SECTION( "qureg types" ) {
            
            pauliOpType paulis[NUM_QUBITS];
            REQUIRE_THROWS_WITH( applyPauliSum(vecIn, paulis, NULL, 1, matOut), Contains("Registers must both be state-vectors or both be density matrices") );
        }
    }
    destroyQureg(vecIn, env);
    destroyQureg(vecOut, env);
    destroyQureg(matIn, env);
    destroyQureg(matOut, env);
    destroyQuESTEnv(env);
}