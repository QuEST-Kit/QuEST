
#include "catch.hpp"
#include "QuEST.h"
#include "utilities.hpp"

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;



/** @sa fromComplex
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "fromComplex", "[data_structures]" ) {
    
    Complex a = {.real=.5, .imag=-.2};
    qcomp b = fromComplex(a);
    
    REQUIRE( a.real == real(b) );
    REQUIRE( a.imag == imag(b) );
}



/** @sa getStaticComplexMatrixN
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "getStaticComplexMatrixN", "[data_structures]" ) {
    
    /* use of this function is illegal in C++ */
    SUCCEED( );
}



/** @sa toComplex
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "toComplex", "[data_structures]" ) {
    
    qcomp a = .5 - .2i;
    Complex b = toComplex(a);
    
    REQUIRE( real(a) == b.real );
    REQUIRE( imag(a) == b.imag );
}



/** @sa createCloneQureg
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createCloneQureg", "[data_structures]" ) {
        
    SECTION( "state-vector" ) {
        
        Qureg a = createQureg(NUM_QUBITS, QUEST_ENV);    
        Qureg b = createCloneQureg(a, QUEST_ENV);
        
        // check properties are the same
        REQUIRE( b.isDensityMatrix == a.isDensityMatrix );
        REQUIRE( b.numQubitsRepresented == a.numQubitsRepresented );
        REQUIRE( b.numQubitsInStateVec == a.numQubitsInStateVec );
        REQUIRE( b.numAmpsPerChunk == a.numAmpsPerChunk );
        REQUIRE( b.numAmpsTotal == a.numAmpsTotal );
        
        // check state-vector is the same (works for GPU and distributed)
        REQUIRE( areEqual(a, b) );  
        
        destroyQureg(a, QUEST_ENV);
        destroyQureg(b, QUEST_ENV);
    }
    SECTION( "density-matrix" ) {
        
        Qureg a = createDensityQureg(NUM_QUBITS, QUEST_ENV);
        Qureg b = createCloneQureg(a, QUEST_ENV);
        
        // check properties are the same
        REQUIRE( b.isDensityMatrix == a.isDensityMatrix );
        REQUIRE( b.numQubitsRepresented == a.numQubitsRepresented );
        REQUIRE( b.numQubitsInStateVec == a.numQubitsInStateVec );
        REQUIRE( b.numAmpsPerChunk == a.numAmpsPerChunk );
        REQUIRE( b.numAmpsTotal == a.numAmpsTotal );
        
        // check state-vector is the same (works for GPU and distributed)
        REQUIRE( areEqual(a, b) );  
        
        destroyQureg(a, QUEST_ENV);
        destroyQureg(b, QUEST_ENV);
    }
}



/** @sa createComplexMatrixN
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createComplexMatrixN", "[data_structures]" ) {
    
    SECTION( "correctness" ) {
        
        int numQb = GENERATE( range(1,10+1) );
        ComplexMatrixN m = createComplexMatrixN(numQb);
        
        // ensure elems are created and initialised to 0
        REQUIRE( areEqual(toQMatrix(m), getZeroMatrix(1<<numQb)) );
        
        destroyComplexMatrixN(m);
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of qubits" ) {
            
            int numQb = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( createComplexMatrixN(numQb), Contains("Invalid number of qubits") );
        }
    }
}



/** @sa createDensityQureg
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createDensityQureg", "[data_structures]" ) {
        
    // must be at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks) - 1; // density matrix has 2*numQb in state-vec
    if (minNumQb <= 0)
        minNumQb = 1;
    
    SECTION( "correctness" ) {
        
        // try 10 valid number of qubits
        int numQb = GENERATE_COPY( range(minNumQb, minNumQb+10) );
        Qureg reg = createDensityQureg(numQb, QUEST_ENV);
        
        // ensure elems (CPU and/or GPU) are created, and reg begins in |0><0|
        QMatrix ref = getZeroMatrix(1<<numQb);
        ref[0][0] = 1; // |0><0|
        REQUIRE( areEqual(reg, ref) );
        
        destroyQureg(reg, QUEST_ENV);
    }
    SECTION( "input validation") {
        
        SECTION( "number of qubits" ) {
            
            int numQb = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( createDensityQureg(numQb, QUEST_ENV), Contains("Invalid number of qubits") );
        }
        SECTION( "number of amplitudes" ) {
            
            // use local QuESTEnv to safely modify
            QuESTEnv env = QUEST_ENV;
            
            // too many amplitudes to store in type
            int maxQb = (int) calcLog2(SIZE_MAX) / 2;
            REQUIRE_THROWS_WITH( createDensityQureg(maxQb+1, env), Contains("Too many qubits") && Contains("size_t type") );
            
            /* n-qubit density matrix contains 2^(2n) amplitudes 
             * so can be spread between at most 2^(2n) ranks
             */
            int minQb = GENERATE_COPY( range(3,maxQb) );
            env.numRanks = (int) pow(2, 2*minQb);
            int numQb = GENERATE_COPY( range(1,minQb) );
            REQUIRE_THROWS_WITH( createDensityQureg(numQb, env), Contains("Too few qubits") );
        }
        SECTION( "available memory" ) {
            
            /* there is no reliable way to force the malloc statements to
             * fail, and hence trigger the matrixInit validation */
            SUCCEED( );
        }
    }
}



/** @sa createPauliHamil
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createPauliHamil", "[data_structures]" ) {

    SECTION( "correctness" ) {

        int numQb = GENERATE( range(1,5) );
        int numTerms = GENERATE( range(1,5) );
        PauliHamil hamil = createPauliHamil(numQb, numTerms);
        
        // check fields are correct
        REQUIRE( hamil.numQubits == numQb );
        REQUIRE( hamil.numSumTerms == numTerms );
        
        // check all Pauli codes are identity
        int numPaulis = numQb * numTerms;
        for (int i=0; i<numPaulis; i++) {
            REQUIRE( hamil.pauliCodes[i] == PAULI_I );
        }
            
        // check all term coefficients can be written to (no seg fault)
        for (int j=0; j<numTerms; j++) {
            hamil.termCoeffs[j] = 1;
            REQUIRE( hamil.termCoeffs[j] == 1 );
        }
        
        destroyPauliHamil(hamil);
    }
    SECTION( "input validation") {

        SECTION( "number of qubits" ) {

            int numQb = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( createPauliHamil(numQb, 1), Contains("The number of qubits and terms in the PauliHamil must be strictly positive.") );
        }
        SECTION( "number of terms" ) {

            int numTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( createPauliHamil(1, numTerms), Contains("The number of qubits and terms in the PauliHamil must be strictly positive.") );
        }
    }
}



/** @sa createPauliHamilFromFile
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createPauliHamilFromFile", "[data_structures]" ) {

    // a file created & populated during the test, and deleted afterward
    char fn[] = "temp_test_output_file.txt";

    SECTION( "correctness" ) {
        
        SECTION( "general" ) {
            
            // for several sizes...
            int numQb = GENERATE( range(1,6) );
            int numTerms = GENERATE( range(1,6) );
            int numPaulis = numQb*numTerms;
            
            // create a PauliHamil with random elements
            qreal coeffs[numTerms];
            enum pauliOpType paulis[numPaulis];
            setRandomPauliSum(coeffs, paulis, numQb, numTerms);
            
            // write the Hamiltonian to file (with trailing whitespace, and trailing newline)
            FILE* file = fopen(fn, "w");
            int i=0;
            for (int n=0; n<numTerms; n++) {
                fprintf(file, "%.20f ", coeffs[n]);
                for (int q=0; q<numQb; q++)
                    fprintf(file, "%d ", (int) paulis[i++]);
                fprintf(file, "\n");
            }
            fprintf(file, "\n");
            fclose(file);
            
            // load the file as a PauliHamil
            PauliHamil hamil = createPauliHamilFromFile(fn);

            // check fields agree
            REQUIRE( hamil.numQubits == numQb );
            REQUIRE( hamil.numSumTerms == numTerms );
            
            // check elements agree
            i=0;
            for (int n=0; n<numTerms; n++) {
                REQUIRE( hamil.termCoeffs[n] == coeffs[n] );
                for (int q=0; q<numQb; q++) {
                    REQUIRE( hamil.pauliCodes[i] == paulis[i] );
                    i++;
                }
            }
            
            destroyPauliHamil(hamil);
        }
        SECTION( "edge cases" ) {
            
            SECTION( "no trailing newline or space" ) {
                
                FILE* file = fopen(fn, "w");
                fprintf(file, ".1 1 0 1");
                fclose(file);
                PauliHamil hamil = createPauliHamilFromFile(fn);
                
                REQUIRE( hamil.numSumTerms == 1 );
                REQUIRE( hamil.numQubits == 3 );
                destroyPauliHamil(hamil); 
            }
            SECTION( "trailing newlines" ) {
                
                FILE* file = fopen(fn, "w");
                fprintf(file, ".1 1 0 1\n\n\n");
                fclose(file);
                PauliHamil hamil = createPauliHamilFromFile(fn);
                
                REQUIRE( hamil.numSumTerms == 1 );
                REQUIRE( hamil.numQubits == 3 );
                destroyPauliHamil(hamil); 
            }
            SECTION( "trailing spaces" ) {
                
                FILE* file = fopen(fn, "w");
                fprintf(file, ".1 1 0 1    ");
                fclose(file);
                PauliHamil hamil = createPauliHamilFromFile(fn);
                
                REQUIRE( hamil.numSumTerms == 1 );
                REQUIRE( hamil.numQubits == 3 );
                destroyPauliHamil(hamil); 
            }
        }
    }
    SECTION( "input validation") {

        SECTION( "number of qubits" ) {

            FILE* file = fopen(fn, "w");
            fprintf(file, ".1 ");
            fclose(file);
            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("The number of qubits") && Contains("strictly positive"));
        }
        SECTION( "coefficient type" ) {

            FILE* file = fopen(fn, "w");
            fprintf(file, "notanumber 1 2 3");
            fclose(file);
            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("Failed to parse") && Contains("coefficient"));
        }
        SECTION( "pauli code" ) {

            // invalid int
            FILE* file = fopen(fn, "w");
            fprintf(file, ".1 1 2 4");
            fclose(file);
            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("invalid pauli code"));
            
            // invalid type
            file = fopen(fn, "w");
            fprintf(file, ".1 1 2 notanumber");
            fclose(file);
            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("Failed to parse the next expected Pauli code"));
        }
    }
    
    // delete the test file
    remove(fn);
}



/** @sa createQuESTEnv
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createQuESTEnv", "[data_structures]" ) {
    
    /* there is no meaningful way to test this */
    SUCCEED( );
}



/** @sa createQureg
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createQureg", "[data_structures]" ) {
        
    // must be at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;
    
    SECTION( "correctness" ) {
        
        // try 10 valid number of qubits
        int numQb = GENERATE_COPY( range(minNumQb, minNumQb+10) );
        Qureg reg = createQureg(numQb, QUEST_ENV);
        
        // ensure elems (CPU and/or GPU) are created, and reg begins in |0>
        QVector ref = QVector(1<<numQb);
        ref[0] = 1; // |0>
        REQUIRE( areEqual(reg, ref) );
        
        destroyQureg(reg, QUEST_ENV);
    }
    SECTION( "input validation") {
        
        SECTION( "number of qubits" ) {
            
            int numQb = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( createQureg(numQb, QUEST_ENV), Contains("Invalid number of qubits") );
        }
        SECTION( "number of amplitudes" ) {
            
            // use local QuESTEnv to safely modify
            QuESTEnv env = QUEST_ENV;
            
            // too many amplitudes to store in type
            int maxQb = (int) calcLog2(SIZE_MAX);
            REQUIRE_THROWS_WITH( createQureg(maxQb+1, env), Contains("Too many qubits") && Contains("size_t type") );
            
            // too few amplitudes to distribute
            int minQb = GENERATE_COPY( range(2,maxQb) );
            env.numRanks = (int) pow(2, minQb);
            int numQb = GENERATE_COPY( range(1,minQb) );
            REQUIRE_THROWS_WITH( createQureg(numQb, env), Contains("Too few qubits") );
        }
        SECTION( "available memory" ) {
            
            /* there is no reliable way to force the malloc statements to
             * fail, and hence trigger the matrixInit validation */
            SUCCEED( );
        }
    }
}



/** @sa destroyComplexMatrixN
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "destroyComplexMatrixN", "[data_structures]" ) {
    
    SECTION( "correctness" ) {
        
        /* there is no meaningful way to test this */
        SUCCEED( );
    }
    SECTION( "input validation" ) {
        
        SECTION( "matrix not created" ) {
            
            /* this is an artificial test case since nothing in the QuEST API 
             * automatically sets un-initialised ComplexMatrixN fields to 
             * the NULL pointer.
             */
             ComplexMatrixN m;
             m.real = NULL;
             
             /* the error message is also somewhat unrelated, but oh well 
              */
             REQUIRE_THROWS_WITH( destroyComplexMatrixN(m), Contains("The ComplexMatrixN was not successfully created") );
        }
    }
}



/** @sa destroyPauliHamil
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "destroyPauliHamil", "[data_structures]" ) {
    
    /* there is no meaningful way to test this.
     * We e.g. cannot check that the pointers are NULL because 
     * they are not updated; this function passes the struct by value,
     * not by reference. We also cannot reliably monitor the 
     * memory used in the heap at runtime.
     */
    SUCCEED( );
}



/** @sa destroyQuESTEnv
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "destroyQuESTEnv", "[data_structures]" ) {

    /* there is no meaningful way to test this */
    SUCCEED( );
}



/** @sa destroyQureg
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "destroyQureg", "[data_structures]" ) {
    
    /* there is no meaningful way to test this.
     * We e.g. cannot check that the pointers are NULL because 
     * they are not updated; this function passes the struct by value,
     * not by reference. We also cannot reliably monitor the 
     * memory used in the heap at runtime.
     */
    SUCCEED( );
}



/** @sa initComplexMatrixN
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initComplexMatrixN", "[data_structures]" ) {
    
    /* use of this function is illegal in C++ */
    SUCCEED( );
}



/** @sa initPauliHamil
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initPauliHamil", "[data_structures]" ) {
    
    SECTION( "correctness" ) {
        
        PauliHamil hamil = createPauliHamil(3, 2);
        
        qreal coeffs[] = {-5, 5};
        enum pauliOpType codes[] = {
            PAULI_X, PAULI_Y, PAULI_Z,
            PAULI_Z, PAULI_Y, PAULI_X};
        initPauliHamil(hamil, coeffs, codes);
        
        // check everything written correctly
        for (int t=0; t<2; t++) {
            REQUIRE( coeffs[t] == hamil.termCoeffs[t] );
            for (int q=0; q<3; q++) {
                int ind = 3*t+q;
                REQUIRE( codes[ind] == hamil.pauliCodes[ind] );
            }
        }
            
        destroyPauliHamil(hamil);
    }
    SECTION( "input validation" ) {
        
        SECTION( "parameters" ) {
            
            // parameters checked before codes, so safe to leave un-initialised
            qreal coeffs[1];
            enum pauliOpType codes[1];
            PauliHamil hamil;
            
            hamil.numQubits = GENERATE( -1, 0 );
            hamil.numSumTerms = 1;
            REQUIRE_THROWS_WITH( initPauliHamil(hamil, coeffs, codes), Contains("number of qubits") && Contains("strictly positive") );
            
            hamil.numQubits = 1;
            hamil.numSumTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( initPauliHamil(hamil, coeffs, codes), Contains("terms") && Contains("strictly positive") );
        }
        SECTION( "Pauli codes" ) {
        
            int numQb = 3;
            int numTerms = 2;
            int numCodes = numQb * numTerms;
            qreal coeffs[numTerms];
            enum pauliOpType codes[numCodes];
            
            // make only one code invalid
            for (int i=0; i<numCodes; i++)
                codes[i] = PAULI_I;
            codes[GENERATE_COPY( range(0,numCodes) )] = (pauliOpType) GENERATE( -1, 4 );
            
            PauliHamil hamil = createPauliHamil(numQb, numTerms);
            REQUIRE_THROWS_WITH( initPauliHamil(hamil, coeffs, codes), Contains("Invalid Pauli code") );
            destroyPauliHamil(hamil);
        }
    }
}