
#include "catch.hpp"
#include "QuEST.h"
#include "utilities.hpp"
    
/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;
    
    

/** @sa cloneQureg
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "cloneQureg", "[state_initialisations]" ) {
    
    Qureg vec1 = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat1 = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            Qureg vec2 = createQureg(NUM_QUBITS, QUEST_ENV);
            
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
            
            destroyQureg(vec2, QUEST_ENV);
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat2 = createDensityQureg(NUM_QUBITS, QUEST_ENV);

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
            
            destroyQureg(mat2, QUEST_ENV);
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qureg type" ) {
            
            REQUIRE_THROWS_WITH( cloneQureg(mat1, vec1),  Contains("both be state-vectors") && Contains("density matrices") );
            REQUIRE_THROWS_WITH( cloneQureg(vec1, mat1),  Contains("both be state-vectors") && Contains("density matrices") );
        }
        SECTION( "qureg dimensions" ) {
            
            Qureg vec3 = createQureg(vec1.numQubitsRepresented + 1, QUEST_ENV);
            Qureg mat3 = createDensityQureg(mat1.numQubitsRepresented + 1, QUEST_ENV);
            
            REQUIRE_THROWS_WITH( cloneQureg(vec1, vec3), Contains("Dimensions") && Contains("don't match") );
            REQUIRE_THROWS_WITH( cloneQureg(mat1, mat3), Contains("Dimensions") && Contains("don't match") );
            
            destroyQureg(vec3, QUEST_ENV);
            destroyQureg(mat3, QUEST_ENV);
        }
    }
    destroyQureg(vec1, QUEST_ENV);
    destroyQureg(mat1, QUEST_ENV);
}



/** @sa initBlankState
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initBlankState", "[state_initialisations]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
        
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
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
}



/** @sa initClassicalState
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initClassicalState", "[state_initialisations]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
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
            REQUIRE_THROWS_WITH( initClassicalState(vec, ind), Contains("Invalid state index") );
        }
    }
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
}



/** @sa initPlusState
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initPlusState", "[state_initialisations]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
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
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
}



/** @sa initPureState
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initPureState", "[state_initialisations]" ) {
    
    Qureg vec1 = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat1 = createDensityQureg(NUM_QUBITS, QUEST_ENV);
        
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            /* state-vector version just performs cloneQureg */
            
            Qureg vec2 = createQureg(NUM_QUBITS, QUEST_ENV);
            
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
            
            destroyQureg(vec2, QUEST_ENV);
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
            REQUIRE_THROWS_WITH( initPureState(vec1, mat1), Contains("Second argument must be a state-vector") );
            REQUIRE_THROWS_WITH( initPureState(mat1, mat1), Contains("Second argument must be a state-vector") );
        }
        SECTION( "qureg dimensions" ) {
            
            Qureg vec2 = createQureg(NUM_QUBITS + 1, QUEST_ENV);
            REQUIRE_THROWS_WITH( initPureState(vec1, vec2), Contains("Dimensions") && Contains("don't match") );
            REQUIRE_THROWS_WITH( initPureState(mat1, vec2), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(vec2, QUEST_ENV);
        }
    }
    destroyQureg(vec1, QUEST_ENV);
    destroyQureg(mat1, QUEST_ENV);
}



/** @sa initStateFromAmps
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initStateFromAmps", "[state_initialisations]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            // create random (unnormalised) vector
            QVector vecRef = getRandomQVector(1<<NUM_QUBITS);
            
            qreal ampsRe[vec.numAmpsTotal];
            qreal ampsIm[vec.numAmpsTotal];
            for (size_t i=0; i<vecRef.size(); i++) {
                ampsRe[i] = real(vecRef[i]);
                ampsIm[i] = imag(vecRef[i]);
            }
            
            initStateFromAmps(vec, ampsRe, ampsIm);
            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            // create random (unnormalised) matrix
            QMatrix matRef = getRandomQMatrix(1<<NUM_QUBITS);
            
            qreal ampsRe[mat.numAmpsTotal];
            qreal ampsIm[mat.numAmpsTotal];
            
            // populate column-wise 
            long long int i=0;
            for (size_t c=0; c<matRef.size(); c++) {
                for (size_t r=0; r<matRef.size(); r++) {
                    ampsRe[i] = real(matRef[r][c]);
                    ampsIm[i] = imag(matRef[r][c]);
                    i++;
                }
            }
    
            initStateFromAmps(mat, ampsRe, ampsIm);
            REQUIRE( areEqual(mat, matRef) );
        }
    }
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
}



/** @sa initZeroState
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "initZeroState", "[state_initialisations]" ) {

    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initBlankState(vec);
            initZeroState(vec);
            
            QVector refVec = QVector(vec.numAmpsTotal);
            refVec[0] = 1;
            REQUIRE( areEqual(vec, refVec) );
        }
        SECTION( "density-matrix" ) {
            
            initBlankState(mat);
            initZeroState(mat);
            
            QMatrix refMat = getZeroMatrix(1<<mat.numQubitsRepresented);
            refMat[0][0] = 1;
            REQUIRE( areEqual(mat, refMat) );
        }
    }
    SECTION( "input validation" ) {
        
        // no input validation 
        SUCCEED( );
    }
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
}



/** @sa setAmps
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "setAmps", "[state_initialisations]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    
    int maxInd = vec.numAmpsTotal;
    qreal reals[maxInd];
    qreal imags[maxInd];
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            // all valid number of amplitudes and offsets
            int startInd = GENERATE_COPY( range(0,maxInd) );
            int numAmps = GENERATE_COPY( range(0,1+maxInd-startInd) ); // upper-bound allows all amps specified
            
            // generate random amplitudes
            for (int i=0; i<numAmps; i++) {
                reals[i] = getRandomReal(-5,5);
                imags[i] = getRandomReal(-5,5);
            }
            
            // check both specified and un-specified amplitudes are correctly handled
            initDebugState(vec);
            QVector vecRef = toQVector(vec);
            
            setAmps(vec, startInd, reals, imags, numAmps);
            for (int i=0; i<numAmps; i++)
                vecRef[startInd+i] = reals[i] + imags[i] * (qcomp) 1i;
                
            REQUIRE( areEqual(vec, vecRef) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "start index" ) {
            
            int startInd = GENERATE_COPY( -1, maxInd );
            int numAmps = 0;
            REQUIRE_THROWS_WITH( setAmps(vec, startInd, reals, imags, numAmps), Contains("Invalid amplitude index") );
        }
        
        SECTION( "number of amplitudes" ) {
            
            // independent
            int startInd = 0;
            int numAmps = GENERATE_COPY( -1, maxInd+1 );
            REQUIRE_THROWS_WITH( setAmps(vec, startInd, reals, imags, numAmps), Contains("Invalid number of amplitudes") );

            // invalid considering start-index
            startInd = maxInd - 1;
            numAmps = 2;
            REQUIRE_THROWS_WITH( setAmps(vec, startInd, reals, imags, numAmps), Contains("More amplitudes given than exist") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
            REQUIRE_THROWS_WITH( setAmps(mat, 0, reals, imags, 0), Contains("valid only for state-vectors") );
            destroyQureg(mat, QUEST_ENV);
        }
    }
    destroyQureg(vec, QUEST_ENV);
}



/** @sa setWeightedQureg
 * @ingroup unittest 
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
            Qureg vecA = createQureg(NUM_QUBITS, QUEST_ENV);
            Qureg vecB = createQureg(NUM_QUBITS, QUEST_ENV);
            Qureg vecC = createQureg(NUM_QUBITS, QUEST_ENV);
            for (int j=0; j<vecA.numAmpsPerChunk; j++) {
                vecA.stateVec.real[j] = getRandomReal(-5,5); vecA.stateVec.imag[j] = getRandomReal(-5,5);
                vecB.stateVec.real[j] = getRandomReal(-5,5); vecB.stateVec.imag[j] = getRandomReal(-5,5);
                vecC.stateVec.real[j] = getRandomReal(-5,5); vecC.stateVec.imag[j] = getRandomReal(-5,5);
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
            Complex facA = toComplex(numA);
            Complex facB = toComplex(numB);
            Complex facC = toComplex(numC);
            
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
            destroyQureg(vecA, QUEST_ENV);
            destroyQureg(vecB, QUEST_ENV);
            destroyQureg(vecC, QUEST_ENV);
        }
        SECTION( "density-matrix" ) {
            
            // make three random matrices
            Qureg matA = createDensityQureg(NUM_QUBITS, QUEST_ENV);
            Qureg matB = createDensityQureg(NUM_QUBITS, QUEST_ENV);
            Qureg matC = createDensityQureg(NUM_QUBITS, QUEST_ENV);
            for (int j=0; j<matA.numAmpsPerChunk; j++) {
                matA.stateVec.real[j] = getRandomReal(-5,5); matA.stateVec.imag[j] = getRandomReal(-5,5);
                matB.stateVec.real[j] = getRandomReal(-5,5); matB.stateVec.imag[j] = getRandomReal(-5,5);
                matC.stateVec.real[j] = getRandomReal(-5,5); matC.stateVec.imag[j] = getRandomReal(-5,5);
            }
            copyStateToGPU(matA); copyStateToGPU(matB); copyStateToGPU(matC);
            QMatrix refA = toQMatrix(matA);
            QMatrix refB = toQMatrix(matB);
            QMatrix refC = toQMatrix(matC);
            QMatrix refOut;
            
            // get three random factors
            qcomp numA = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            qcomp numB = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            qcomp numC = getRandomReal(-5,5) + getRandomReal(-5,5) * (qcomp) 1i;
            Complex facA = toComplex(numA);
            Complex facB = toComplex(numB);
            Complex facC = toComplex(numC);
            
            // check out-qureg is correct, when all quregs are unique... 
            setWeightedQureg(facA, matA, facB, matB, facC, matC);
            refOut = numA*refA + numB*refB + numC*refC;
            REQUIRE( areEqual(matC, refOut) );
            
            // ... and that other qureg's aren't modified
            REQUIRE( areEqual(matA, refA) );
            REQUIRE( areEqual(matB, refB) );
            
            // check quregOut correct, when it's also qureg2
            refC = toQMatrix(matC);
            setWeightedQureg(facB, matB, facC, matC, facA, matC);
            refOut = numB*refB + numC*refC + numA*refC;
            REQUIRE( areEqual(matC, refOut, 10*REAL_EPS) );
            
            // ... and that the remaining qureg is not modified
            REQUIRE( areEqual(matB, refB) );
            
            // check quregOut correct, when it's also qureg1
            refC = toQMatrix(matC);
            setWeightedQureg(facC, matC, facB, matB, facA, matC);
            refOut = numC*refC + numB*refB + numA*refC;
            REQUIRE( areEqual(matC, refOut, 1E2*REAL_EPS) );
            
            // ... and that the remaining qureg is not modified
            REQUIRE( areEqual(matB, refB) );
            
            // check quregOut is correct when it's both input quregs
            refC = toQMatrix(matC);
            setWeightedQureg(facA, matC, facB, matC, facC, matC);
            refOut = numA*refC + numB*refC + numC*refC;
            REQUIRE( areEqual(matC, refOut, 1E3*REAL_EPS) );
        
            // cleanup
            destroyQureg(matA, QUEST_ENV);
            destroyQureg(matB, QUEST_ENV);
            destroyQureg(matC, QUEST_ENV);
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qureg types" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
            Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
            Complex f = {.real=0, .imag=0};
            
            // two state-vecs, one density-matrix
            REQUIRE_THROWS_WITH( setWeightedQureg(f, mat, f, vec, f, vec), Contains("state-vectors or") && Contains("density matrices") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vec, f, mat, f, vec), Contains("state-vectors or") && Contains("density matrices") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vec, f, vec, f, mat), Contains("state-vectors or") && Contains("density matrices") );

            // one state-vec, two density-matrices
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vec, f, mat, f, mat), Contains("state-vectors or") && Contains("density matrices") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, mat, f, vec, f, mat), Contains("state-vectors or") && Contains("density matrices") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, mat, f, mat, f, vec), Contains("state-vectors or") && Contains("density matrices") );
        
            destroyQureg(vec, QUEST_ENV);
            destroyQureg(mat, QUEST_ENV);
        } 
        SECTION( "qureg dimensions" ) {
            
            Qureg vecA = createQureg(NUM_QUBITS, QUEST_ENV);
            Qureg vecB = createQureg(NUM_QUBITS + 1, QUEST_ENV);
            Qureg matA = createDensityQureg(NUM_QUBITS, QUEST_ENV);
            Qureg matB = createDensityQureg(NUM_QUBITS + 1, QUEST_ENV);
            Complex f = {.real=0, .imag=0};
            
            // state-vecs
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vecA, f, vecB, f, vecB), Contains("Dimensions") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vecB, f, vecA, f, vecB), Contains("Dimensions") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, vecB, f, vecB, f, vecA), Contains("Dimensions") );
            
            // density-matrices
            REQUIRE_THROWS_WITH( setWeightedQureg(f, matA, f, matB, f, matB), Contains("Dimensions") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, matB, f, matA, f, matB), Contains("Dimensions") );
            REQUIRE_THROWS_WITH( setWeightedQureg(f, matB, f, matB, f, matA), Contains("Dimensions") );
            
            destroyQureg(vecA, QUEST_ENV);
            destroyQureg(vecB, QUEST_ENV);
            destroyQureg(matA, QUEST_ENV);
            destroyQureg(matB, QUEST_ENV);
        }
    }
}

