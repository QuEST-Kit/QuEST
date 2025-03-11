/** @file
 * Ported tests of the deprecated QuEST v3 interface,
 * unit testing the "state initialisation" module.
 * 
 * This file should be excluded from doxygen parsing so 
 * as not to conflict with the doc of the v4 unit tests.
 * 
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#define INCLUDE_DEPRECATED_FUNCTIONS 1
#define DISABLE_DEPRECATION_WARNINGS 1
#include "quest/include/quest.h"

#include "test_utilities.hpp"
    
/* allows concise use of ContainsSubstring in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::ContainsSubstring;
    
    

/** @sa cloneQureg
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "cloneQureg", "[state_initialisations]" ) {
    
    Qureg vec1 = createForcedQureg(NUM_QUBITS);
    Qureg mat1 = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            Qureg vec2 = createForcedQureg(NUM_QUBITS);
            
            // make sure states start differently
            initDebugState(vec1);    
            initBlankState(vec2);
            REQUIRE( !areEqual(vec1, vec2) );
            
            // make sure vec2 is changed
            QVector copy1 = toQVector(vec1);
            cloneQureg(vec2, vec1);
            REQUIRE( areEqual(vec1, vec2) );
            
            // make sure vec1 unaffected
            REQUIRE( areEqual(vec1, copy1) );
            
            destroyQureg(vec2);
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat2 = createForcedDensityQureg(NUM_QUBITS);

            // make sure states start differently
            initDebugState(mat1);
            initBlankState(mat2);
            REQUIRE( !areEqual(mat1, mat2) );
            
            // make sure vec2 is changed
            QMatrix copy1 = toQMatrix(mat1);
            cloneQureg(mat2, mat1);
            REQUIRE( areEqual(mat1, mat2) );
            
            // make sure vec1 unaffected
            REQUIRE( areEqual(mat1, copy1) );
            
            destroyQureg(mat2);
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qureg type" ) {
            
            REQUIRE_THROWS_WITH( cloneQureg(mat1, vec1), ContainsSubstring("must both be statevectors, or both be density matrices") );
            REQUIRE_THROWS_WITH( cloneQureg(vec1, mat1), ContainsSubstring("must both be statevectors, or both be density matrices") );
        }
        SECTION( "qureg dimensions" ) {
            
            Qureg vec3 = createForcedQureg(vec1.numQubits + 1);
            Qureg mat3 = createForcedDensityQureg(mat1.numQubits + 1);
            
            REQUIRE_THROWS_WITH( cloneQureg(vec1, vec3), ContainsSubstring("different number of qubits") );
            REQUIRE_THROWS_WITH( cloneQureg(mat1, mat3), ContainsSubstring("different number of qubits") );
            
            destroyQureg(vec3);
            destroyQureg(mat3);
        }
    }
    destroyQureg(vec1);
    destroyQureg(mat1);
}



/** @sa initBlankState
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "initBlankState", "[state_initialisations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
        
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initBlankState(vec);
            REQUIRE( areEqual(vec, QVector(1<<NUM_QUBITS)) );
        }
        SECTION( "density-matrix" ) {
            
            initBlankState(mat);
            REQUIRE( areEqual(mat, getZeroMatrix(1<<NUM_QUBITS)) );
        }
    }
    SECTION( "input validation" ) {
        
        // no user validation
        SUCCEED( );
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa initClassicalState
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "initClassicalState", "[state_initialisations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        int numInds = (1<<NUM_QUBITS);
        int ind = GENERATE_COPY( range(0,numInds) );
        
        SECTION( "state-vector" ) {
            
            initClassicalState(vec, ind);
            QVector vecRef = QVector(1<<NUM_QUBITS);
            vecRef[ind] = 1;
            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            initClassicalState(mat, ind);
            QMatrix matRef = getZeroMatrix(1<<NUM_QUBITS);
            matRef[ind][ind] = 1;
            REQUIRE( areEqual(mat, matRef) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, (1<<NUM_QUBITS) );
            REQUIRE_THROWS_WITH( initClassicalState(vec, ind), ContainsSubstring("Basis state index") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa initPlusState
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "initPlusState", "[state_initialisations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            // |+> = 1/sqrt(N^2) sum_i |i>
            //     = 1/sqrt(N^2) {1, ..., 1}
            initPlusState(vec);
            QVector vecRef = QVector(1<<NUM_QUBITS);
            for (size_t i=0; i<vecRef.size(); i++)
                vecRef[i] = 1./sqrt(pow(2,NUM_QUBITS));
            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            // |+><+| = 1/sqrt(N^2) sum_i |i> 1/sqrt(N^2) sum_j <j|
            //        = 1/(N^2) sum_{ij} |i><j|
            //        = 1/(N^2) {{1, ..., 1}, ...}
            initPlusState(mat);
            QMatrix matRef = getZeroMatrix(1<<NUM_QUBITS);
            for (size_t i=0; i<matRef.size(); i++)
                for (size_t j=0; j<matRef.size(); j++)
                    matRef[i][j] = 1./pow(2, NUM_QUBITS);
            REQUIRE( areEqual(mat, matRef) );
        }
    }
    SECTION( "input validation" ) {
        
        // no user validation
        SUCCEED( );
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa initPureState
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "initPureState", "[state_initialisations]" ) {
    
    Qureg vec1 = createForcedQureg(NUM_QUBITS);
    Qureg mat1 = createForcedDensityQureg(NUM_QUBITS);
        
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            /* state-vector version just performs cloneQureg */
            
            Qureg vec2 = createForcedQureg(NUM_QUBITS);
            
            // make sure states start differently
            initDebugState(vec1);    
            initBlankState(vec2);
            REQUIRE( !areEqual(vec1, vec2) );
            
            // make sure vec2 is overwritten with vec1
            QVector copy1 = toQVector(vec1);
            initPureState(vec2, vec1);
            REQUIRE( areEqual(vec1, vec2) );
            
            // make sure vec1 was not modified 
            REQUIRE( areEqual(vec1, copy1) );
            
            destroyQureg(vec2);
        }
        SECTION( "density-matrix" ) {
            
            /* density matrix version initialises matrix in |pure><pure| */
            
            initDebugState(vec1); // |vec1> = sum_i a_i |i>
            QVector copy1 = toQVector(vec1);
            
            // make sure mat1 is modified correctly
            initBlankState(mat1); 
            initPureState(mat1, vec1); // mat1 = |vec1><vec1| = sum_{ij} a_i a_j* |i><j|
            
            QMatrix matRef = getZeroMatrix(1<<NUM_QUBITS);
            for (size_t i=0; i<matRef.size(); i++)
                for (size_t j=0; j<matRef.size(); j++)
                    matRef[i][j] = copy1[i] * conj(copy1[j]);
            REQUIRE( areEqual(mat1, matRef) );
            
            // make sure vec1 was not modified
            REQUIRE( areEqual(vec1, copy1) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qureg types" ) {
            
            // density matrix as second arg is illegal (regardless of first arg)
            REQUIRE_THROWS_WITH( initPureState(vec1, mat1), ContainsSubstring("(the second argument) must be a statevector") );
            REQUIRE_THROWS_WITH( initPureState(mat1, mat1), ContainsSubstring("(the second argument) must be a statevector") );
        }
        SECTION( "qureg dimensions" ) {
            
            Qureg vec2 = createForcedQureg(NUM_QUBITS + 1);
            REQUIRE_THROWS_WITH( initPureState(vec1, vec2), ContainsSubstring("differing number of qubits") );
            REQUIRE_THROWS_WITH( initPureState(mat1, vec2), ContainsSubstring("differing number of qubits") );
            destroyQureg(vec2);
        }
    }
    destroyQureg(vec1);
    destroyQureg(mat1);
}



/** @sa initArbitraryPureState
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "initArbitraryPureState", "[state_initialisations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);

    // create random (unnormalised) vector
    QVector vecRef = getRandomQVector(1<<NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {

            initArbitraryPureState(vec, vecRef.data());

            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {

            initArbitraryPureState(mat, vecRef.data());

            // matRef = |vecRef><vecRef|
            QMatrix matRef = getZeroMatrix(vecRef.size());
            for (size_t r=0; r<vecRef.size(); r++)
                for (size_t c=0; c<vecRef.size(); c++)
                    matRef[r][c] = vecRef[r] * conj(vecRef[c]);

            REQUIRE( areEqual(mat, matRef) );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa initZeroState
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "initZeroState", "[state_initialisations]" ) {

    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initBlankState(vec);
            initZeroState(vec);
            
            QVector refVec = QVector(vec.numAmps);
            refVec[0] = 1;
            REQUIRE( areEqual(vec, refVec) );
        }
        SECTION( "density-matrix" ) {
            
            initBlankState(mat);
            initZeroState(mat);
            
            QMatrix refMat = getZeroMatrix(1<<mat.numQubits);
            refMat[0][0] = 1;
            REQUIRE( areEqual(mat, refMat) );
        }
    }
    SECTION( "input validation" ) {
        
        // no input validation 
        SUCCEED( );
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa setAmps
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "setQuregAmps", "[state_initialisations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    
    int maxInd = vec.numAmps;
    QVector amps = getRandomQVector(maxInd);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            // all valid number of amplitudes and offsets
            int startInd = GENERATE_COPY( range(0,maxInd) );
            int numAmps = GENERATE_COPY( range(0,1+maxInd-startInd) ); // upper-bound allows all amps specified
            
            // check both specified and un-specified amplitudes are correctly handled
            initDebugState(vec);
            QVector vecRef = toQVector(vec);

            setQuregAmps(vec, startInd, amps.data(), numAmps);
            for (int i=0; i<numAmps; i++)
                vecRef[startInd+i] = amps[i];
                
            REQUIRE( areEqual(vec, vecRef) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "start index" ) {
            
            int startInd = GENERATE_COPY( -1, maxInd );
            int numAmps = 0;
            REQUIRE_THROWS_WITH( setQuregAmps(vec, startInd, amps.data(), numAmps), ContainsSubstring("starting basis state index") );
        }
        
        SECTION( "number of amplitudes" ) {
            
            // independent
            int startInd = 0;
            int numAmps = GENERATE_COPY( -1, maxInd+1 );
            REQUIRE_THROWS_WITH( setQuregAmps(vec, startInd, amps.data(), numAmps), ContainsSubstring("number of amplitudes") );

            // invalid considering start-index
            startInd = maxInd - 1;
            numAmps = 2;
            REQUIRE_THROWS_WITH( setQuregAmps(vec, startInd, amps.data(), numAmps), ContainsSubstring("implies an end index") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createForcedDensityQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( setQuregAmps(mat, 0, amps.data(), 0), ContainsSubstring("Expected a statevector Qureg but received a density matrix") );
            destroyQureg(mat);
        }
    }
    destroyQureg(vec);
}



// setDensityAmps removed because replacement setDensityQuregAmps
// has an incompatible signature (accepting qcomps directly). Note
// I do not know why we updated setQuregAmps above to the v4 API :)

    // /** @sa setDensityAmps
    //  * @ingroup deprecatedtests 
    //  * @author Tyson Jones 
    //  */
    // TEST_CASE( "setDensityAmps", "[state_initialisations]" ) {
        
    //     Qureg matr = createForcedDensityQureg(NUM_QUBITS);
    //     QVector amps = getRandomQVector(matr.numAmps);
    //     int maxInd = (1 << NUM_QUBITS);
        
    //     SECTION( "correctness" ) {
            
    //         SECTION( "density-matrix" ) {
                
    //             // try all valid number of amplitudes and offsets
    //             int startRow = GENERATE_COPY( range(0,maxInd) );
    //             int startCol = GENERATE_COPY( range(0,maxInd) );
                
    //             // determine the max number of amps that can be passed from the given start indices
    //             int numPriorAmps = startRow + startCol*(1 << matr.numQubits);
    //             int maxNumAmps = matr.numAmps - numPriorAmps;

    //             // previously, we tried all possible number of amps, like so:
    //             //      int numAmps = GENERATE_COPY( range(0,maxNumAmps) ); // upper-bound allows all amps specified
    //             // but this is too many and causes a timeout on Ubuntu. 
    //             // So we instead randomly choose the number of amps, and repeat 10 times.
    //             int numAmps = getRandomInt(0, maxNumAmps);
    //             GENERATE_COPY( range(0,10) );
                
    //             // check both specified and un-specified amplitudes are correctly handled
    //             initDebugState(matr);
    //             QMatrix matrRef = toQMatrix(matr);
                
    //             setDensityQuregAmps(matr, startRow, startCol, reals, imags, numAmps);
                
    //             int r=startRow;
    //             int c=startCol;
    //             for (int i=0; i<numAmps; i++) {
    //                 qcomp amp = reals[i] + imags[i] * (qcomp) 1i;
    //                 matrRef[r][c] = amp;
                    
    //                 r++;
    //                 if (r >= maxInd ) {
    //                     r=0;
    //                     c++;
    //                 }
    //             }
                    
    //             REQUIRE( areEqual(matr, matrRef) );
    //         }
    //     }
    //     SECTION( "input validation" ) {
            
    //         SECTION( "start index" ) {
                
    //             int badInd = GENERATE_COPY( -1, maxInd );
    //             int numAmps = 0;
    //             REQUIRE_THROWS_WITH( setDensityAmps(matr, badInd, 0, reals, imags, numAmps), ContainsSubstring("Invalid amplitude index") );
    //             REQUIRE_THROWS_WITH( setDensityAmps(matr, 0, badInd, reals, imags, numAmps), ContainsSubstring("Invalid amplitude index") );
    //         }
    //         SECTION( "number of amplitudes" ) {
                
    //             // independent
    //             int numAmps = GENERATE_COPY( -1, matr.numAmps+1 );
    //             REQUIRE_THROWS_WITH( setDensityAmps(matr, 0, 0, reals, imags, numAmps), ContainsSubstring("Invalid number of amplitudes") );

    //             // invalid considering start-index
    //             REQUIRE_THROWS_WITH( setDensityAmps(matr, maxInd-1, maxInd-1, reals, imags, 2), ContainsSubstring("More amplitudes given than exist") );
    //             REQUIRE_THROWS_WITH( setDensityAmps(matr, maxInd-1, maxInd-2, reals, imags, maxInd+2), ContainsSubstring("More amplitudes given than exist") );
    //         }
    //         SECTION( "state-vector" ) {
                
    //             Qureg vec = createForcedQureg(NUM_QUBITS);
    //             REQUIRE_THROWS_WITH( setDensityAmps(vec, 0, 0, reals, imags, 0), ContainsSubstring("valid only for density matrices") );
    //             destroyQureg(vec);
    //         }
    //     }
    //     destroyQureg(matr);
    // }



// setQuregToPauliHamil removed because replacement setQuregToPauliStrSum
// accepts new type PauliStrSum which has an initialiser incompatible
// with the deprecated PauliHamil type

    // /** @sa setQuregToPauliHamil
    //  * @ingroup deprecatedtests 
    //  * @author Tyson Jones 
    //  */
    // TEST_CASE( "setQuregToPauliHamil", "[state_initialisations]" ) {
        
    //     Qureg rho = createForcedDensityQureg(NUM_QUBITS);
        
    //     SECTION( "correctness" ) {
            
    //         // too expensive to enumerate all paulis, so try 10x with random ones
    //         GENERATE( range(0,10) );
            
    //         int numTerms = GENERATE( 1, 2, 10, 15, 100, 1000 );
    //         PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);
    //         setRandomPauliSum(hamil);

    //         setQuregToPauliHamil(rho, hamil);
    //         REQUIRE( areEqual(rho, toQMatrix(hamil)) );
            
    //         destroyPauliHamil(hamil);
    //     }
    //     SECTION( "input validation" ) {
            
    //         SECTION( "density-matrix" ) {
                
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS, 5);
    //             Qureg vec = createForcedQureg(NUM_QUBITS);
                
    //             REQUIRE_THROWS_WITH( setQuregToPauliHamil(vec, hamil), ContainsSubstring("density matrices") );
                
    //             destroyQureg(vec);
    //             destroyPauliHamil(hamil);
    //         }
    //         SECTION( "pauli codes" ) {
                
    //             int numTerms = 3;
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);
                
    //             // make one pauli code wrong
    //             hamil.pauliCodes[GENERATE_COPY( range(0,numTerms*NUM_QUBITS) )] = (pauliOpType) GENERATE( -1, 4 );
    //             REQUIRE_THROWS_WITH( setQuregToPauliHamil(rho, hamil), ContainsSubstring("Invalid Pauli code") );
                
    //             destroyPauliHamil(hamil);
    //         }
    //         SECTION( "matching hamiltonian qubits" ) {
                
    //             int numTerms = 1;
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS + 1, numTerms);
                
    //             REQUIRE_THROWS_WITH( setQuregToPauliHamil(rho, hamil), ContainsSubstring("same number of qubits") );
                
    //             destroyPauliHamil(hamil);
    //         }
    //     }
    //     destroyQureg(rho);
    // }



/** @sa setWeightedQureg
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "setWeightedQureg", "[state_initialisations]" ) {
        
    SECTION( "correctness" ) {
        
        // repeat each test below 10 times 
        GENERATE( range(0,10) );
        
        /* note tolerance in areEqual increases with tests, since 
         * small differences propogate in vecC which is not re-initialised 
         */
        
        SECTION( "state-vector" ) {
            
            // make three random vectors
            Qureg vecA = createForcedQureg(NUM_QUBITS);
            Qureg vecB = createForcedQureg(NUM_QUBITS);
            Qureg vecC = createForcedQureg(NUM_QUBITS);
            for (int j=0; j<vecA.numAmpsPerNode; j++) {
                vecA.cpuAmps[j] = getRandomComplex();
                vecB.cpuAmps[j] = getRandomComplex();
                vecC.cpuAmps[j] = getRandomComplex();
            }
            copyStateToGPU(vecA); copyStateToGPU(vecB); copyStateToGPU(vecC);
            QVector refA = toQVector(vecA);
            QVector refB = toQVector(vecB);
            QVector refC = toQVector(vecC);
            QVector refOut;
            
            // get three random factors
            qcomp numA = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            qcomp numB = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            qcomp numC = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            Complex facA; facA.real = real(numA); facA.imag = imag(numA);
            Complex facB; facB.real = real(numB); facB.imag = imag(numB);
            Complex facC; facC.real = real(numC); facC.imag = imag(numC);
            
            // check out-qureg is correct, when all quregs are unique... 
            setWeightedQureg(facA, vecA, facB, vecB, facC, vecC);
            refOut = numA*refA + numB*refB + numC*refC;
            REQUIRE( areEqual(vecC, refOut) );
            
            // ... and that other qureg's aren't modified
            REQUIRE( areEqual(vecA, refA) );
            REQUIRE( areEqual(vecB, refB) );
            
            // check quregOut correct, when it's also qureg2
            refC = toQVector(vecC);
            setWeightedQureg(facB, vecB, facC, vecC, facA, vecC);
            refOut = numB*refB + numC*refC + numA*refC;
            REQUIRE( areEqual(vecC, refOut, 10*REAL_EPS) );
            
            // ... and that the remaining qureg is not modified
            REQUIRE( areEqual(vecB, refB) );
            
            // check quregOut correct, when it's also qureg1
            refC = toQVector(vecC);
            setWeightedQureg(facC, vecC, facB, vecB, facA, vecC);
            refOut = numC*refC + numB*refB + numA*refC;
            REQUIRE( areEqual(vecC, refOut, 10*REAL_EPS) );
            
            // ... and that the remaining qureg is not modified
            REQUIRE( areEqual(vecB, refB) );
            
            // check quregOut is correct when it's both input quregs
            refC = toQVector(vecC);
            setWeightedQureg(facA, vecC, facB, vecC, facC, vecC);
            refOut = numA*refC + numB*refC + numC*refC;
            REQUIRE( areEqual(vecC, refOut, 1E3*REAL_EPS) );
        
            // cleanup
            destroyQureg(vecA);
            destroyQureg(vecB);
            destroyQureg(vecC);
        }

        // v4 does not permit superposing of density matrices

            // SECTION( "density-matrix" ) {
                
            //     // make three random matrices
            //     Qureg matA = createForcedDensityQureg(NUM_QUBITS);
            //     Qureg matB = createForcedDensityQureg(NUM_QUBITS);
            //     Qureg matC = createForcedDensityQureg(NUM_QUBITS);
            //     for (int j=0; j<matA.numAmpsPerNode; j++) {
            //         matA.cpuAmps[j] = getRandomComplex();
            //         matB.cpuAmps[j] = getRandomComplex();
            //         matC.cpuAmps[j] = getRandomComplex();
            //     }
            //     copyStateToGPU(matA); copyStateToGPU(matB); copyStateToGPU(matC);
            //     QMatrix refA = toQMatrix(matA);
            //     QMatrix refB = toQMatrix(matB);
            //     QMatrix refC = toQMatrix(matC);
            //     QMatrix refOut;
                
            //     // get three random factors
            //     qcomp numA = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            //     qcomp numB = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            //     qcomp numC = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            //     Complex facA; facA.real = real(numA); facA.imag = imag(numA);
            //     Complex facB; facB.real = real(numB); facB.imag = imag(numB);
            //     Complex facC; facC.real = real(numC); facC.imag = imag(numC);
                
            //     // check out-qureg is correct, when all quregs are unique... 
            //     setWeightedQureg(facA, matA, facB, matB, facC, matC);
            //     refOut = numA*refA + numB*refB + numC*refC;
            //     REQUIRE( areEqual(matC, refOut) );
                
            //     // ... and that other qureg's aren't modified
            //     REQUIRE( areEqual(matA, refA) );
            //     REQUIRE( areEqual(matB, refB) );
                
            //     // check quregOut correct, when it's also qureg2
            //     refC = toQMatrix(matC);
            //     setWeightedQureg(facB, matB, facC, matC, facA, matC);
            //     refOut = numB*refB + numC*refC + numA*refC;
            //     REQUIRE( areEqual(matC, refOut, 10*REAL_EPS) );
                
            //     // ... and that the remaining qureg is not modified
            //     REQUIRE( areEqual(matB, refB) );
                
            //     // check quregOut correct, when it's also qureg1
            //     refC = toQMatrix(matC);
            //     setWeightedQureg(facC, matC, facB, matB, facA, matC);
            //     refOut = numC*refC + numB*refB + numA*refC;
            //     REQUIRE( areEqual(matC, refOut, 1E2*REAL_EPS) );
                
            //     // ... and that the remaining qureg is not modified
            //     REQUIRE( areEqual(matB, refB) );
                
            //     // check quregOut is correct when it's both input quregs
            //     refC = toQMatrix(matC);
            //     setWeightedQureg(facA, matC, facB, matC, facC, matC);
            //     refOut = numA*refC + numB*refC + numC*refC;
            //     REQUIRE( areEqual(matC, refOut, 1E3*REAL_EPS) );
            
            //     // cleanup
            //     destroyQureg(matA);
            //     destroyQureg(matB);
            //     destroyQureg(matC);
            // }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qureg types" ) {
            
            Qureg vec = createForcedQureg(NUM_QUBITS);
            Qureg mat = createForcedDensityQureg(NUM_QUBITS);
            Complex f; f.real = 0; f.imag = 0;
            
            // two state-vecs, one density-matrix
            REQUIRE_THROWS_WITH( setWeightedQureg(f, mat, f, vec, f, vec), ContainsSubstring("Cannot superpose a density matrix. All quregs must be statevectors") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vec, f, mat, f, vec), ContainsSubstring("Cannot superpose a density matrix. All quregs must be statevectors") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vec, f, vec, f, mat), ContainsSubstring("Cannot superpose a density matrix. All quregs must be statevectors") );

            // one state-vec, two density-matrices
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vec, f, mat, f, mat), ContainsSubstring("Cannot superpose a density matrix. All quregs must be statevectors") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, mat, f, vec, f, mat), ContainsSubstring("Cannot superpose a density matrix. All quregs must be statevectors") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, mat, f, mat, f, vec), ContainsSubstring("Cannot superpose a density matrix. All quregs must be statevectors") );
        
            destroyQureg(vec);
            destroyQureg(mat);
        } 
        SECTION( "qureg dimensions" ) {
            
            Qureg vecA = createForcedQureg(NUM_QUBITS);
            Qureg vecB = createForcedQureg(NUM_QUBITS + 1);
            Qureg matA = createForcedDensityQureg(NUM_QUBITS);
            Qureg matB = createForcedDensityQureg(NUM_QUBITS + 1);
            Complex f; f.real = 0; f.imag = 0;
            
            // state-vecs
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vecA, f, vecB, f, vecB), ContainsSubstring("differing numbers of qubits") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vecB, f, vecA, f, vecB), ContainsSubstring("differing numbers of qubits") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vecB, f, vecB, f, vecA), ContainsSubstring("differing numbers of qubits") );
            
            // v4 does not permit superposing density matrices

                // // density-matrices
                // REQUIRE_THROWS_WITH( setWeightedQureg(f, matA, f, matB, f, matB), ContainsSubstring("differing numbers of qubits") );
                // REQUIRE_THROWS_WITH( setWeightedQureg(f, matB, f, matA, f, matB), ContainsSubstring("differing numbers of qubits") );
                // REQUIRE_THROWS_WITH( setWeightedQureg(f, matB, f, matB, f, matA), ContainsSubstring("differing numbers of qubits") );
                
            destroyQureg(vecA);
            destroyQureg(vecB);
            destroyQureg(matA);
            destroyQureg(matB);
        }
    }
}
