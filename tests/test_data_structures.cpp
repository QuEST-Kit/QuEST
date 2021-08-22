
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
    
    qcomp a = qcomp(.5,-.2);
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



/** @sa createDiagonalOp
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createDiagonalOp", "[data_structures]" ) {
    
    // must be at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;
        
    SECTION( "correctness" ) {
        
        // try 10 valid number of qubits
        int numQb = GENERATE_COPY( range(minNumQb, minNumQb+10) );
        DiagonalOp op = createDiagonalOp(numQb, QUEST_ENV);
        
        // check properties are correct
        REQUIRE( op.numQubits == numQb );
        REQUIRE( op.chunkId == QUEST_ENV.rank );
        REQUIRE( op.numChunks == QUEST_ENV.numRanks );
        REQUIRE( op.numElemsPerChunk == (1LL << numQb) / QUEST_ENV.numRanks );
        REQUIRE( op.real != NULL );
        REQUIRE( op.imag != NULL );
        
        // check all elements in CPU are zero 
        REQUIRE( areEqual(toQVector(op), QVector(1LL << numQb)) );
        
        // (no concise way to check this for GPU)
        
        destroyDiagonalOp(op, QUEST_ENV);
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of qubits" ) {
            
            int numQb = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( createDiagonalOp(numQb, QUEST_ENV), Contains("Invalid number of qubits") );
        }
        SECTION( "number of elements" ) {
            
            // use local QuESTEnv to safely modify
            QuESTEnv env = QUEST_ENV;
            
            // too many amplitudes to store in type
            int maxQb = (int) calcLog2(SIZE_MAX);
            REQUIRE_THROWS_WITH( createDiagonalOp(maxQb+1, env), Contains("Too many qubits") && Contains("size_t type") );
            
            // too few amplitudes to distribute
            int minQb = GENERATE_COPY( range(2,maxQb) );
            env.numRanks = (int) pow(2, minQb);
            int numQb = GENERATE_COPY( range(1,minQb) );
            REQUIRE_THROWS_WITH( createDiagonalOp(numQb, env), Contains("Too few qubits") && Contains("distributed"));
        }
        SECTION( "available memory" ) {
            
            /* there is no reliable way to force the malloc statements to
             * fail, and hence trigger the diagonalOpInit validation */
            SUCCEED( );
        }
    }
}



/** @sa createDiagonalOpFromPauliHamilFile
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "createDiagonalOpFromPauliHamilFile", "[data_structures]" ) {

    // files created & populated during the test, and deleted afterward
    char fnPrefix[] = "temp_createDiagonalOpFromPauliHamilFile";
    char fn[100];
    
    // each test uses a unique filename (managed by master node), to avoid file IO locks
    // (it's safe for sub-test to overwrite this, since after each test, all files
    //  with prefix fnPrefix are deleted)
    setUniqueFilename(fn, fnPrefix);
    
    // diagonal op must have at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;

    SECTION( "correctness" ) {

        SECTION( "general" ) {
            
            // try several Pauli Hamiltonian sizes
            int numQb = GENERATE_COPY( range(minNumQb, 6+minNumQb) );
            int numTerms = GENERATE_COPY( 1, minNumQb, 10*minNumQb );
            
            // create a PauliHamil with random elements
            PauliHamil hamil = createPauliHamil(numQb, numTerms);
            setRandomDiagPauliHamil(hamil);
            
            // write the Hamiltonian to file (with trailing whitespace, and trailing newline)
            if (QUEST_ENV.rank == 0) { 
                FILE* file = fopen(fn, "w");
                int i=0;
                for (int n=0; n<numTerms; n++) {
                    fprintf(file, REAL_STRING_FORMAT, hamil.termCoeffs[n]);
                    fprintf(file, " ");
                    for (int q=0; q<numQb; q++)
                        fprintf(file, "%d ", (int) hamil.pauliCodes[i++]);
                    fprintf(file, "\n");
                }
                fprintf(file, "\n");
                fclose(file);
            }
            syncQuESTEnv(QUEST_ENV);
            
            // load the file as a diagonal operator, and compare
            DiagonalOp op = createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV);
            REQUIRE( areEqual(toQMatrix(op), toQMatrix(hamil)) );

            destroyPauliHamil(hamil);
            destroyDiagonalOp(op, QUEST_ENV);
        }
        SECTION( "edge cases" ) {
            
            // prepare a valid single-term diagonal Pauli Hamiltonian
            qreal coeffs[] = {.1};
            pauliOpType codes[minNumQb];
            for (int q=0; q<minNumQb; q++)
                codes[q] = (q%2)? PAULI_I : PAULI_Z; 
                
            QMatrix ref = toQMatrix(coeffs, codes, minNumQb, 1);

            // prepare basic encoding string 
            string line = to_string(coeffs[0]) + " ";
            for (int q=0; q<minNumQb; q++)
                line += to_string(codes[q]) + ((q<minNumQb-1)? " ":"");
            
            SECTION( "no trailing newline or space" ) {
                
                writeToFileSynch(fn, line);
                
                DiagonalOp op = createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV);
                REQUIRE( areEqual(ref, toQMatrix(op)) );
                
                destroyDiagonalOp(op, QUEST_ENV);
            }
            SECTION( "trailing newlines" ) {
                
                writeToFileSynch(fn, line + "\n\n\n");
                
                DiagonalOp op = createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV);
                REQUIRE( areEqual(ref, toQMatrix(op)) );
                
                destroyDiagonalOp(op, QUEST_ENV);
            }
            SECTION( "trailing spaces" ) {
                
                writeToFileSynch(fn, line + "    ");
                
                DiagonalOp op = createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV);
                REQUIRE( areEqual(ref, toQMatrix(op)) );
                
                destroyDiagonalOp(op, QUEST_ENV);
            }
        }
    }
    SECTION( "input validation") {

        SECTION( "number of qubits" ) {
            
            writeToFileSynch(fn, ".1 "); // 0 qubits
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("The number of qubits") && Contains("strictly positive"));
        }
        SECTION( "number of elements" ) {
            
            // too many amplitudes to store in type
            int maxQb = (int) calcLog2(SIZE_MAX);
            
            // encode one more qubit than legal to file
            string line = ".1 ";
            for (int q=0; q<(maxQb+1); q++)
                line += "3 "; // trailing space ok
            writeToFileSynch(fn, line);
            
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("Too many qubits") && Contains("size_t type") );
            
            // use local QuESTEnv to safely modify
            QuESTEnv env = QUEST_ENV;
            
            // too few elements to distribute
            int minQb = GENERATE_COPY( range(2,maxQb) );
            env.numRanks = (int) pow(2, minQb);
            int numQb = GENERATE_COPY( range(1,minQb) );
            
            line = ".1 ";
            for (int q=0; q<numQb; q++)
                line += "3 "; // trailing space ok
            setUniqueFilename(fn, fnPrefix);
            writeToFileSynch(fn, line);
            
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, env), Contains("Too few qubits") && Contains("distributed") );
        }
        SECTION( "coefficient type" ) {

            writeToFileSynch(fn, "notanumber 1 2 3");
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("Failed to parse") && Contains("coefficient"));
        }
        SECTION( "pauli code" ) {

            writeToFileSynch(fn, ".1 0 3 2");  // final is invalid Y
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("contained operators other than PAULI_Z and PAULI_I"));
            
            setUniqueFilename(fn, fnPrefix);
            writeToFileSynch(fn, ".1 0 1 3");  // second is invalid X
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("contained operators other than PAULI_Z and PAULI_I"));
            
            setUniqueFilename(fn, fnPrefix);
            writeToFileSynch(fn, ".1 0 1 4");  // final is invalid Pauli code
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("invalid pauli code"));
            
            setUniqueFilename(fn, fnPrefix);
            writeToFileSynch(fn, ".1 3 0 notanumber");  // final is invalid type
            REQUIRE_THROWS_WITH( createDiagonalOpFromPauliHamilFile(fn, QUEST_ENV), Contains("Failed to parse the next expected Pauli code"));
        }
    }
    
    // delete all files created above
    deleteFilesWithPrefixSynch(fnPrefix);
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
    
    // files created & populated during the test, and deleted afterward
    char fnPrefix[] = "temp_createPauliHamilFromFile";
    char fn[100];
    
    // each test uses a unique filename (managed by master node), to avoid file IO locks
    // (it's safe for sub-test to overwrite this, since after each test, all files
    //  with prefix fnPrefix are deleted)
    setUniqueFilename(fn, fnPrefix);

    SECTION( "correctness" ) {
        
        SECTION( "general" ) {
            
            // for several sizes...
            int numQb = GENERATE( 1, 5, 10, 15 );
            int numTerms = GENERATE( 1, 10, 30 );
            int numPaulis = numQb*numTerms;
            
            // create a PauliHamil with random elements
            qreal coeffs[numTerms];
            enum pauliOpType paulis[numPaulis];
            setRandomPauliSum(coeffs, paulis, numQb, numTerms);
            
            // write the Hamiltonian to file (with trailing whitespace, and trailing newline)
            if (QUEST_ENV.rank == 0) {
                FILE* file = fopen(fn, "w");
                int i=0;
                for (int n=0; n<numTerms; n++) {
                    fprintf(file, REAL_STRING_FORMAT, coeffs[n]);
                    fprintf(file, " ");
                    for (int q=0; q<numQb; q++)
                        fprintf(file, "%d ", (int) paulis[i++]);
                    fprintf(file, "\n");
                }
                fprintf(file, "\n");
                fclose(file);
            }
            syncQuESTEnv(QUEST_ENV);
            
            // load the file as a PauliHamil
            PauliHamil hamil = createPauliHamilFromFile(fn);

            // check fields agree
            REQUIRE( hamil.numQubits == numQb );
            REQUIRE( hamil.numSumTerms == numTerms );
            
            // check elements agree
            int j=0;
            for (int n=0; n<numTerms; n++) {
                REQUIRE( absReal(hamil.termCoeffs[n] - coeffs[n]) <= REAL_EPS );
                for (int q=0; q<numQb; q++) {
                    REQUIRE( hamil.pauliCodes[j] == paulis[j] );
                    j++;
                }
            }
            
            destroyPauliHamil(hamil);
        }
        SECTION( "edge cases" ) {
            
            SECTION( "no trailing newline or space" ) {
                
                writeToFileSynch(fn, ".1 1 0 1");
        
                PauliHamil hamil = createPauliHamilFromFile(fn);
                REQUIRE( hamil.numSumTerms == 1 );
                REQUIRE( hamil.numQubits == 3 );
                
                destroyPauliHamil(hamil); 
            }
            SECTION( "trailing newlines" ) {
                
                writeToFileSynch(fn, ".1 1 0 1\n\n\n");

                PauliHamil hamil = createPauliHamilFromFile(fn);
                REQUIRE( hamil.numSumTerms == 1 );
                REQUIRE( hamil.numQubits == 3 );
                
                destroyPauliHamil(hamil); 
            }
            SECTION( "trailing spaces" ) {
                
                writeToFileSynch(fn, ".1 1 0 1    ");
                
                PauliHamil hamil = createPauliHamilFromFile(fn);
                REQUIRE( hamil.numSumTerms == 1 );
                REQUIRE( hamil.numQubits == 3 );
                
                destroyPauliHamil(hamil); 
            }
        }
    }
    SECTION( "input validation") {

        SECTION( "number of qubits" ) {

            writeToFileSynch(fn, ".1 ");

            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("The number of qubits") && Contains("strictly positive"));
        }
        SECTION( "coefficient type" ) {

            writeToFileSynch(fn, "notanumber 1 2 3");

            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("Failed to parse") && Contains("coefficient"));
        }
        SECTION( "pauli code" ) {
            
            writeToFileSynch(fn, ".1 1 2 4"); // invalid int

            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("invalid pauli code"));
            
            setUniqueFilename(fn, fnPrefix);
            writeToFileSynch(fn, ".1 1 2 notanumber"); // invalid type
            
            REQUIRE_THROWS_WITH( createPauliHamilFromFile(fn), Contains("Failed to parse the next expected Pauli code"));
        }
    }
    
    // cleanup temp files
    deleteFilesWithPrefixSynch(fn);
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



/** @sa destroyDiagonalOp
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "destroyDiagonalOp", "[data_structures]" ) {

    /* there is no meaningful way to test this */
    SUCCEED( );
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



/** @sa initDiagonalOp
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initDiagonalOp", "[data_structures]" ) {
    
    // must be at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;
    
    SECTION( "correctness" ) {
        
        // try 10 valid number of qubits
        int numQb = GENERATE_COPY( range(minNumQb, minNumQb+10) );
        DiagonalOp op = createDiagonalOp(numQb, QUEST_ENV);
        
        long long int len = (1LL << numQb);
        qreal reals[len];
        qreal imags[len];
        long long int n;
        for (n=0; n<len; n++) {
            reals[n] = (qreal)    n;
            imags[n] = (qreal) -2*n; // (n - 2n i)
        }
        initDiagonalOp(op, reals, imags);
        
        // check that op.real and op.imag were modified...
        REQUIRE( areEqual(toQVector(op), reals, imags) );
        
        // and also that GPU real and imag were modified
        // via if it modifies an all-unity state-vec correctly 
        Qureg qureg = createQureg(numQb, QUEST_ENV);
        for (long long int i=0; i<qureg.numAmpsPerChunk; i++) {
            qureg.stateVec.real[i] = 1;
            qureg.stateVec.imag[i] = 1;
        }
        copyStateToGPU(qureg);
        
        QVector prodRef = toQMatrix(op) * toQVector(qureg);
        
        // (n - 2n i) * (1 + 1i) = 3n - n*i
        applyDiagonalOp(qureg, op);
        copyStateFromGPU(qureg);
        QVector result = toQVector(qureg);    
        REQUIRE( areEqual(prodRef, result) );
        
        destroyQureg(qureg, QUEST_ENV);
        destroyDiagonalOp(op, QUEST_ENV);
    }
}



/** @sa initDiagonalOpFromPauliHamil
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initDiagonalOpFromPauliHamil", "[data_structures]" ) {
    
    // distributed diagonal op must contain at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;
    
    SECTION( "correctness" ) {    

        // try (at most) 10 valid number of qubits (even for validation)
        int numQb = GENERATE_COPY( range(minNumQb, min(10,minNumQb+10)) );
        DiagonalOp op = createDiagonalOp(numQb, QUEST_ENV);
        
        // try several sized random all-Z Hamiltonians
        int numTerms = GENERATE_COPY( 1, numQb, 5*numQb );
        PauliHamil hamil = createPauliHamil(numQb, numTerms);
        setRandomDiagPauliHamil(hamil);

        initDiagonalOpFromPauliHamil(op, hamil);
        REQUIRE( areEqual(toQMatrix(op), toQMatrix(hamil)) );

        destroyPauliHamil(hamil);
        destroyDiagonalOp(op, QUEST_ENV);
    }
    SECTION( "input validation" ) {
        
        SECTION( "hamiltonian parameters" ) {
            
            DiagonalOp op = createDiagonalOp(minNumQb, QUEST_ENV);
            PauliHamil hamil;
            
            hamil.numQubits = GENERATE( -1, 0 );
            hamil.numSumTerms = 1;
            REQUIRE_THROWS_WITH( initDiagonalOpFromPauliHamil(op, hamil), Contains("number of qubits") && Contains("strictly positive") );
            
            hamil.numQubits = minNumQb;
            hamil.numSumTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( initDiagonalOpFromPauliHamil(op, hamil), Contains("terms") && Contains("strictly positive") );
            
            destroyDiagonalOp(op, QUEST_ENV);
        }
        SECTION( "mismatching dimensions" ) {
            
            int numQbA = minNumQb+1;
            int numQbB = GENERATE_COPY( numQbA-1, numQbA+1 );
            
            DiagonalOp op = createDiagonalOp(numQbA, QUEST_ENV);
            PauliHamil hamil = createPauliHamil(numQbB, 1);
            
            REQUIRE_THROWS_WITH( initDiagonalOpFromPauliHamil(op, hamil), Contains("Pauli Hamiltonian and diagonal operator have different, incompatible dimensions") );
            
            destroyDiagonalOp(op, QUEST_ENV);
            destroyPauliHamil(hamil);
        }
        SECTION( "pauli codes" ) {
        
            DiagonalOp op = createDiagonalOp(minNumQb, QUEST_ENV);
            PauliHamil hamil = createPauliHamil(minNumQb, 5);  // all I
            
            // make only one code invalid
            int numCodes = minNumQb * hamil.numSumTerms;
            int ind = GENERATE_COPY( range(0,numCodes) );
            hamil.pauliCodes[ind] = GENERATE( PAULI_X, PAULI_Y );
            
            REQUIRE_THROWS_WITH( initDiagonalOpFromPauliHamil(op, hamil), Contains("contained operators other than PAULI_Z and PAULI_I") );
            
            destroyDiagonalOp(op, QUEST_ENV);
            destroyPauliHamil(hamil);
        }
    }
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



/** @sa setDiagonalOpElems
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "setDiagonalOpElems", "[data_structures]" ) {
    
    // must be at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;
    
    // try 10 valid number of qubits (even for validation)
    int numQb = GENERATE_COPY( range(minNumQb, minNumQb+10) );
    DiagonalOp op = createDiagonalOp(numQb, QUEST_ENV);
    
    SECTION( "correctness" ) {
    
        // make entire array on every node
        long long int len = (1LL << numQb);
        qreal reals[len];
        qreal imags[len];
        long long int n;
        for (n=0; n<len; n++) {
            reals[n] = (qreal)    n;
            imags[n] = (qreal) -2*n; // (n - 2n i)
        }
        
        // set one value at a time (only relevant nodes will update)
        for (n=0; n<len; n++)
            setDiagonalOpElems(op, n, &reals[n], &imags[n], 1);
        
        // check op.real and op.imag updated correctly 
        REQUIRE( areEqual(toQVector(op), reals, imags) );
        
        // no check that GPU values updated (occurs in initDiagonalOp)
    }
    SECTION( "input validation" ) {
        
        long long int maxInd = (1LL << numQb);
        qreal *reals;
        qreal *imags;
        
        SECTION( "start index" ) {
            
            int startInd = GENERATE_COPY( -1, maxInd );
            int numAmps = 1;
            REQUIRE_THROWS_WITH( setDiagonalOpElems(op, startInd, reals, imags, numAmps), Contains("Invalid element index") );
        }
        
        SECTION( "number of elements" ) {
            
            // independent
            int startInd = 0;
            int numAmps = GENERATE_COPY( -1, maxInd+1 );
            REQUIRE_THROWS_WITH( setDiagonalOpElems(op, startInd, reals, imags, numAmps), Contains("Invalid number of elements") );

            // invalid considering start-index
            startInd = maxInd - 1;
            numAmps = 2;
            REQUIRE_THROWS_WITH( setDiagonalOpElems(op, startInd, reals, imags, numAmps), Contains("More elements given than exist") );
        }
    }
    
    destroyDiagonalOp(op, QUEST_ENV);
}



/** @sa syncDiagonalOp
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "syncDiagonalOp", "[data_structures]" ) {
    
    // must be at least one amplitude per node
    int minNumQb = calcLog2(QUEST_ENV.numRanks);
    if (minNumQb == 0)
        minNumQb = 1;
        
    SECTION( "correctness" ) {
        
        // try 10 valid number of qubits
        int numQb = GENERATE_COPY( range(minNumQb, minNumQb+10) );
        DiagonalOp op = createDiagonalOp(numQb, QUEST_ENV);
        
        // check that changes get sync'd to the GPU...
        long long int n;
        for (n=0; n<op.numElemsPerChunk; n++) {
            op.real[n] = (qreal)    n;
            op.imag[n] = (qreal) -2*n; // (n - 2n i)
        }
        syncDiagonalOp(op);
        
        // via if it modifies an all-unity state-vec correctly 
        Qureg qureg = createQureg(numQb, QUEST_ENV);
        for (long long int i=0; i<qureg.numAmpsPerChunk; i++) {
            qureg.stateVec.real[i] = 1;
            qureg.stateVec.imag[i] = 1;
        }
        copyStateToGPU(qureg);
        
        // (n - 2n i) * (1 + 1i) = 3n - n*i
        applyDiagonalOp(qureg, op);
        copyStateFromGPU(qureg);
        for (n=0; n<qureg.numAmpsPerChunk; n++) {
            REQUIRE( qureg.stateVec.real[n] == 3*n );
            REQUIRE( qureg.stateVec.imag[n] == -n );
        }
        
        destroyQureg(qureg, QUEST_ENV);
        destroyDiagonalOp(op, QUEST_ENV);
    }
}