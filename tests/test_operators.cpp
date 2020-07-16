
#include "catch.hpp"
#include "QuEST.h"
#include "utilities.hpp"

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;



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