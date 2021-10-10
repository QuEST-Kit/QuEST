
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



/** @sa applyDiagonalOp
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyDiagonalOp", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    SECTION( "correctness" ) {
        
        // try 10 random operators 
        GENERATE( range(0,10) );
        
        // make a totally random (non-Hermitian) diagonal oeprator
        DiagonalOp op = createDiagonalOp(NUM_QUBITS, QUEST_ENV);
        for (long long int i=0; i<op.numElemsPerChunk; i++) {
            op.real[i] = getRandomReal(-5, 5);
            op.imag[i] = getRandomReal(-5, 5);
        }
        syncDiagonalOp(op);
        
        SECTION( "state-vector" ) {
            
            QVector ref = toQMatrix(op) * refVec;
            applyDiagonalOp(quregVec, op);
            REQUIRE( areEqual(quregVec, ref) );
        }
        SECTION( "density-matrix" ) {
            
            QMatrix ref = toQMatrix(op) * refMatr;
            applyDiagonalOp(quregMatr, op);
            REQUIRE( areEqual(quregMatr, ref, 100*REAL_EPS) );
        }
        
        destroyDiagonalOp(op, QUEST_ENV);
    }
    SECTION( "input validation" ) {
        
        SECTION( "mismatching size" ) {
            
            DiagonalOp op = createDiagonalOp(NUM_QUBITS + 1, QUEST_ENV);
            
            REQUIRE_THROWS_WITH( applyDiagonalOp(quregVec, op), Contains("equal number of qubits"));
            REQUIRE_THROWS_WITH( applyDiagonalOp(quregMatr, op), Contains("equal number of qubits"));
            
            destroyDiagonalOp(op, QUEST_ENV);
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}

/** @sa applyFullQFT
 * @ingroup unittest
 * @author Tyson Jones 
 */
TEST_CASE( "applyFullQFT", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    SECTION( "correctness" ) {
        
        // try 10 times for each kind of input
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            SECTION( "normalised" ) {
                
                refVec = getRandomStateVector(NUM_QUBITS);
                toQureg(quregVec, refVec);
                refVec = getDFT(refVec);
                
                applyFullQFT(quregVec);
                REQUIRE( areEqual(quregVec, refVec) );
            }
            SECTION( "unnormalised" ) {
                
                refVec = getRandomQVector(1 << NUM_QUBITS);
                toQureg(quregVec, refVec);
                refVec = getDFT(refVec);
                
                applyFullQFT(quregVec);
                REQUIRE( areEqual(quregVec, refVec) );
            }
        }
        SECTION( "density-matrix" ) {
            
            SECTION( "pure" ) {
                
                /* a pure density matrix should be mapped to a pure state 
                 * corresponding to the state-vector DFT
                 */

                refVec = getRandomStateVector(NUM_QUBITS);
                refMatr = getPureDensityMatrix(refVec);
                
                toQureg(quregMatr, refMatr);
                applyFullQFT(quregMatr);
                
                refVec = getDFT(refVec);
                refMatr = getPureDensityMatrix(refVec);
                
                REQUIRE( areEqual(quregMatr, refMatr) );
            }
            SECTION( "mixed" ) {
                
                /* a mixed density matrix, conceptualised as a mixture of orthogonal
                 * state-vectors, should be mapped to an equally weighted mixture 
                 * of DFTs of each state-vector (because QFT is unitary and hence 
                 * maintains state orthogonality)
                 */
                
                int numStates = (1 << NUM_QUBITS)/4; // quarter as many states as possible
                std::vector<QVector> states = getRandomOrthonormalVectors(NUM_QUBITS, numStates);
                std::vector<qreal> probs = getRandomProbabilities(numStates);
                
                // set qureg to random mixture 
                refMatr = getMixedDensityMatrix(probs, states);
                toQureg(quregMatr, refMatr);
                
                // apply QFT to mixture
                applyFullQFT(quregMatr);
                
                // compute dft of mixture, via dft of each state
                refMatr = getZeroMatrix(1 << NUM_QUBITS);
                for (int i=0; i<numStates; i++)
                    refMatr += probs[i] * getPureDensityMatrix(getDFT(states[i]));
                
                REQUIRE( areEqual(quregMatr, refMatr) );
            }
            SECTION( "unnormalised" ) {
                
                /* repeat method above, except that we use unnormalised vectors, 
                 * and mix them with arbitrary complex numbers instead of probabilities,
                 * yielding an unnormalised density matrix 
                 */
                
                int numVecs = (1 << NUM_QUBITS)/4; // quarter as many states as possible
                std::vector<QVector> vecs;
                std::vector<qcomp> coeffs;
                for (int i=0; i<numVecs; i++) {
                    vecs.push_back(getRandomQVector(1 << NUM_QUBITS));
                    coeffs.push_back(getRandomComplex());
                }
                
                // produce unnormalised matrix
                refMatr = getZeroMatrix(1 << NUM_QUBITS);
                for (int i=0; i<numVecs; i++)
                    refMatr += coeffs[i] * getPureDensityMatrix(vecs[i]);
                    
                toQureg(quregMatr, refMatr);
                applyFullQFT(quregMatr);
                
                // compute target matrix via dft of each unnormalised vector 
                refMatr = getZeroMatrix(1 << NUM_QUBITS);
                for (int i=0; i<numVecs; i++)
                    refMatr += coeffs[i] * getPureDensityMatrix(getDFT(vecs[i]));
                    
                REQUIRE( areEqual(quregMatr, refMatr) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SUCCEED( );
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



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



/** @sa applyMultiVarPhaseFunc
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyMultiVarPhaseFunc", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encodings
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every possible number of registers 
        // (between #qubits containing 1, and 1 containing #qubits)
        int numRegs;
        int maxNumRegs = 0;
        if (encoding == UNSIGNED)
            maxNumRegs = NUM_QUBITS;
        if (encoding == TWOS_COMPLEMENT)
            maxNumRegs = NUM_QUBITS/2;  // floors
        numRegs = GENERATE_COPY( range(1, maxNumRegs+1) );
        
        // try every possible total number of involed qubits
        int totalNumQubits;
        int minTotalQubits = 0;
        if (encoding == UNSIGNED)
            // each register must contain at least 1 qubit
            minTotalQubits = numRegs;
        if (encoding == TWOS_COMPLEMENT)
            // each register must contain at least 2 qubits
            minTotalQubits = 2*numRegs;
        totalNumQubits = GENERATE_COPY( range(minTotalQubits,NUM_QUBITS+1) );
                        
        // try every qubits subset and ordering 
        int* regs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), totalNumQubits) );
        
        // assign each sub-reg its minimum length
        int unallocQubits = totalNumQubits;
        int numQubitsPerReg[numRegs];
        for (int i=0; i<numRegs; i++) 
            if (encoding == UNSIGNED) {
                numQubitsPerReg[i] = 1;
                unallocQubits -= 1;
            }
            else if (encoding == TWOS_COMPLEMENT) {
                numQubitsPerReg[i] = 2;
                unallocQubits -= 2;
            }
        // and randomly allocate the remaining qubits between the registers
        while (unallocQubits > 0) {
            numQubitsPerReg[getRandomInt(0,numRegs)] += 1;
            unallocQubits--;
        }
        

        // each register gets a random number of terms in the full phase function
        int numTermsPerReg[numRegs];
        int numTermsTotal = 0;
        for (int r=0; r<numRegs; r++) {
            int numTerms = getRandomInt(1,5);
            numTermsPerReg[r] = numTerms;
            numTermsTotal += numTerms;
        }
            
        // populate the multi-var phase function with random but POSITIVE-power terms,
        // which must further be integers in two's complement
        int flatInd = 0;
        qreal coeffs[numTermsTotal];
        qreal expons[numTermsTotal];
        for (int r=0; r<numRegs; r++) {
            for (int t=0; t<numTermsPerReg[r]; t++) {
                coeffs[flatInd] = getRandomReal(-10,10);
                if (encoding == TWOS_COMPLEMENT)
                    expons[flatInd] = getRandomInt(0, 3+1);
                else if (encoding == UNSIGNED)
                    expons[flatInd] = getRandomReal(0, 3);
                    
                flatInd++;
            }
        }
        
        /* To perform this calculation more distinctly from the QuEST 
         * core method, we can exploit that 
         * exp(i (f1[x1] + f2[x2]))|x2>|x1> = exp(i f2[x2])|x2> (x) exp(i f1[x1])|x1> 
         * and so compute a separate diagonal matrix for each sub-register with 
         * phases determined only by its own variable, and Kronecker multiply 
         * them together. The target qubits of the Kronecker product is then 
         * just the full register, stored in *regs. 
         */
        QMatrix allRegMatr{{1}};
        int startInd = 0;
        
        for (int r=0; r<numRegs; r++) {
            
            QMatrix singleRegMatr = getZeroMatrix( 1 << numQubitsPerReg[r] );
            
            for (size_t i=0; i<singleRegMatr.size(); i++) {
                
                long long int ind = 0;
                if (encoding == UNSIGNED)
                    ind = i;
                if (encoding == TWOS_COMPLEMENT)
                    ind = getTwosComplement(i, numQubitsPerReg[r]);
                
                qreal phase = 0;
                for (int t=0; t<numTermsPerReg[r]; t++)
                    phase += coeffs[t+startInd] * pow(ind, expons[t+startInd]);
                    
                singleRegMatr[i][i] = expI(phase);
            }                
            allRegMatr = getKroneckerProduct(singleRegMatr, allRegMatr);
            startInd += numTermsPerReg[r];
        }
        
        SECTION( "state-vector" ) {
                
            applyMultiVarPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, coeffs, expons, numTermsPerReg);
            applyReferenceOp(refVec, regs, totalNumQubits, allRegMatr);
            REQUIRE( areEqual(quregVec, refVec, 1E4*REAL_EPS) );
        }
        SECTION( "density-matrix" ) {
            
            applyMultiVarPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, coeffs, expons, numTermsPerReg);
            applyReferenceOp(refMatr, regs, totalNumQubits, allRegMatr);
            REQUIRE( areEqual(quregMatr, refMatr, 1E6*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        int numRegs = 2;
        int numQubitsPerReg[] = {2,3};
        int qubits[] = {0,1,2,3,4};
        
        SECTION( "number of registers" ) {
            
            numRegs = GENERATE_COPY( -1, 0, 1+MAX_NUM_REGS_APPLY_ARBITRARY_PHASE );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL), Contains("Invalid number of qubit subregisters") );
        }
        SECTION( "number of qubits" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0, 1+NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            qubits[GENERATE(2,3,4)] = qubits[1];
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {
            
            qubits[GENERATE(range(0,NUM_QUBITS))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL), Contains("Invalid qubit index") );
        }
        SECTION( "number of terms" ) {

            int numTermsPerReg[] = {3, 3};
            
            numTermsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, numTermsPerReg), Contains("Invalid number of terms in the phase function") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, enc, NULL, NULL, NULL), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = 1;
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, NULL, NULL, NULL), Contains("A sub-register contained too few qubits to employ TWOS_COMPLEMENT encoding") );
        }
        SECTION( "fractional exponent" ) {
            
            int numTermsPerReg[] = {3, 3};
            qreal coeffs[] = {0,0,0,  0,0,0};
            qreal expos[] =  {1,2,3,  1,2,3};
            
            expos[GENERATE(range(0,6))] = GENERATE( 0.5, 1.999, 5.0001 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, coeffs, expos, numTermsPerReg), Contains("The phase function contained a fractional exponent, which is illegal in TWOS_COMPLEMENT") );
            
            // ensure fractional exponents are valid in unsigned mode however
            REQUIRE_NOTHROW( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, coeffs, expos, numTermsPerReg) );

        }
        SECTION( "negative exponent" ) {

            int numTermsPerReg[] = {3, 3};
            qreal coeffs[] = {0,0,0,  0,0,0};
            qreal expos[] =  {1,2,3,  1,2,3};
            
            expos[GENERATE(range(0,6))] = GENERATE( -1, -2, -2.5 );
            enum bitEncoding enc = GENERATE( UNSIGNED, TWOS_COMPLEMENT );
    
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFunc(quregVec, qubits, numQubitsPerReg, numRegs, enc, coeffs, expos, numTermsPerReg), Contains("The phase function contained an illegal negative exponent") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyMultiVarPhaseFuncOverrides
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyMultiVarPhaseFuncOverrides", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encodings
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every possible number of registers 
        // (between #qubits containing 1, and 1 containing #qubits)
        int numRegs;
        int maxNumRegs = 0;
        if (encoding == UNSIGNED)
            maxNumRegs = NUM_QUBITS;
        if (encoding == TWOS_COMPLEMENT)
            maxNumRegs = NUM_QUBITS/2;  // floors
        numRegs = GENERATE_COPY( range(1, maxNumRegs+1) );
        
        // try every possible total number of involed qubits
        int totalNumQubits;
        int minTotalQubits = 0;
        if (encoding == UNSIGNED)
            // each register must contain at least 1 qubit
            minTotalQubits = numRegs;
        if (encoding == TWOS_COMPLEMENT)
            // each register must contain at least 2 qubits
            minTotalQubits = 2*numRegs;
        totalNumQubits = GENERATE_COPY( range(minTotalQubits,NUM_QUBITS+1) );
                        
        // try every qubits subset and ordering 
        int* regs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), totalNumQubits) );
        
        // assign each sub-reg its minimum length
        int unallocQubits = totalNumQubits;
        int numQubitsPerReg[numRegs];
        for (int i=0; i<numRegs; i++) 
            if (encoding == UNSIGNED) {
                numQubitsPerReg[i] = 1;
                unallocQubits -= 1;
            }
            else if (encoding == TWOS_COMPLEMENT) {
                numQubitsPerReg[i] = 2;
                unallocQubits -= 2;
            }
        // and randomly allocate the remaining qubits between the registers
        while (unallocQubits > 0) {
            numQubitsPerReg[getRandomInt(0,numRegs)] += 1;
            unallocQubits--;
        }
        

        // each register gets a random number of terms in the full phase function
        int numTermsPerReg[numRegs];
        int numTermsTotal = 0;
        for (int r=0; r<numRegs; r++) {
            int numTerms = getRandomInt(1,5);
            numTermsPerReg[r] = numTerms;
            numTermsTotal += numTerms;
        }
            
        // populate the multi-var phase function with random but POSITIVE-power terms,
        // which must further be integers in two's complement
        int flatInd = 0;
        qreal coeffs[numTermsTotal];
        qreal expons[numTermsTotal];
        for (int r=0; r<numRegs; r++) {
            for (int t=0; t<numTermsPerReg[r]; t++) {
                coeffs[flatInd] = getRandomReal(-10,10);
                if (encoding == TWOS_COMPLEMENT)
                    expons[flatInd] = getRandomInt(0, 3+1);
                else if (encoding == UNSIGNED)
                    expons[flatInd] = getRandomReal(0, 3);
                    
                flatInd++;
            }
        }
        
        
        // choose a random number of overrides (even overriding every amplitude)
        int numOverrides = getRandomInt(0, (1<<totalNumQubits) + 1);
        
        // randomise each override index (uniqueness isn't checked)
        long long int overrideInds[numOverrides*numRegs];
        flatInd = 0;
        for (int v=0; v<numOverrides; v++) {
            for (int r=0; r<numRegs; r++) {
                if (encoding == UNSIGNED)
                    overrideInds[flatInd] = getRandomInt(0, 1<<numQubitsPerReg[r]);
                else if (encoding == TWOS_COMPLEMENT) 
                    overrideInds[flatInd] = getRandomInt(-(1<<(numQubitsPerReg[r]-1)), (1<<(numQubitsPerReg[r]-1))-1);
                flatInd++;
            }
        }

        // override to a random phase
        qreal overridePhases[numOverrides];
        for (int v=0; v<numOverrides; v++)
            overridePhases[v] = getRandomReal(-4, 4); // periodic in [-pi, pi]
            
        /* To perform this calculation more distinctly from the QuEST 
         * core method, we can exploit that 
         * exp(i (f1[x1] + f2[x2]))|x2>|x1> = exp(i f2[x2])|x2> (x) exp(i f1[x1])|x1> 
         * and so compute a separate diagonal matrix for each sub-register with 
         * phases determined only by its own variable, and Kronecker multiply 
         * them together. The target qubits of the Kronecker product is then 
         * just the full register, stored in *regs. We do this, then iterate 
         * the list of overrides directly to overwrite the ultimate diagonal 
         * entries
         */
        QMatrix allRegMatr{{1}};
        int startInd = 0;
        for (int r=0; r<numRegs; r++) {
            
            QMatrix singleRegMatr = getZeroMatrix( 1 << numQubitsPerReg[r] );
            
            for (size_t i=0; i<singleRegMatr.size(); i++) {
                
                long long int ind = 0;
                if (encoding == UNSIGNED)
                    ind = i;
                if (encoding == TWOS_COMPLEMENT)
                    ind = getTwosComplement(i, numQubitsPerReg[r]);
                
                qreal phase = 0;
                for (int t=0; t<numTermsPerReg[r]; t++)
                    phase += coeffs[t+startInd] * pow(ind, expons[t+startInd]);
                    
                singleRegMatr[i][i] = expI(phase);
            }
            allRegMatr = getKroneckerProduct(singleRegMatr, allRegMatr);
            startInd += numTermsPerReg[r];
        }    
        setDiagMatrixOverrides(allRegMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
        
        SECTION( "state-vector" ) {
                    
            applyMultiVarPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, coeffs, expons, numTermsPerReg, overrideInds, overridePhases, numOverrides);
            applyReferenceOp(refVec, regs, totalNumQubits, allRegMatr);
            REQUIRE( areEqual(quregVec, refVec, 1E4*REAL_EPS) );
        }
        SECTION( "density-matrix" ) {
            
            applyMultiVarPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, coeffs, expons, numTermsPerReg, overrideInds, overridePhases, numOverrides);
            applyReferenceOp(refMatr, regs, totalNumQubits, allRegMatr);
            REQUIRE( areEqual(quregMatr, refMatr, 1E6*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        int numRegs = 2;
        int numQubitsPerReg[] = {2,3};
        int qubits[] = {0,1,2,3,4};
        
        SECTION( "number of registers" ) {
            
            numRegs = GENERATE_COPY( -1, 0, 1+MAX_NUM_REGS_APPLY_ARBITRARY_PHASE );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL, NULL, NULL, 0), Contains("Invalid number of qubit subregisters") );
        }
        SECTION( "number of qubits" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0, 1+NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL, NULL, NULL, 0), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            qubits[GENERATE(2,3,4)] = qubits[1];
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL, NULL, NULL, 0), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {
            
            qubits[GENERATE(range(0,NUM_QUBITS))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, NULL, NULL, NULL, 0), Contains("Invalid qubit index") );
        }
        SECTION( "number of terms" ) {

            int numTermsPerReg[] = {3, 3};
            
            numTermsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, NULL, NULL, numTermsPerReg, NULL, NULL, 0), Contains("Invalid number of terms in the phase function") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, enc, NULL, NULL, NULL, NULL, NULL, 0), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = 1;
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, NULL, NULL, NULL, NULL, NULL, 0), Contains("A sub-register contained too few qubits to employ TWOS_COMPLEMENT encoding") );
        }
        SECTION( "fractional exponent" ) {
            
            int numTermsPerReg[] = {3, 3};
            qreal coeffs[] = {0,0,0,  0,0,0};
            qreal expos[] =  {1,2,3,  1,2,3};
            
            expos[GENERATE(range(0,6))] = GENERATE( 0.5, 1.999, 5.0001 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, coeffs, expos, numTermsPerReg, NULL, NULL, 0), Contains("The phase function contained a fractional exponent, which is illegal in TWOS_COMPLEMENT") );
            
            // ensure fractional exponents are valid in unsigned mode however
            REQUIRE_NOTHROW( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, coeffs, expos, numTermsPerReg, NULL, NULL, 0) );
        }
        SECTION( "negative exponent" ) {

            int numTermsPerReg[] = {3, 3};
            qreal coeffs[] = {0,0,0,  0,0,0};
            qreal expos[] =  {1,2,3,  1,2,3};
            
            expos[GENERATE(range(0,6))] = GENERATE( -1, -2, -2.5 );
            enum bitEncoding enc = GENERATE( UNSIGNED, TWOS_COMPLEMENT );
    
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, enc, coeffs, expos, numTermsPerReg, NULL, NULL, 0), Contains("The phase function contained an illegal negative exponent") );
        }
        SECTION( "number of overrides" ) {
            
            int numTermsPerReg[] = {3, 3};
            qreal coeffs[] = {0,0,0,  0,0,0};
            qreal expos[] =  {1,2,3,  1,2,3};
   
            int numOverrides = -1;
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, UNSIGNED, coeffs, expos, numTermsPerReg, NULL, NULL, numOverrides), Contains("Invalid number of phase function overrides specified") );
        }
        SECTION( "override indices" ) {
            
            int numTermsPerReg[] = {3, 3};
            qreal coeffs[] = {0,0,0,  0,0,0};
            qreal expos[] =  {1,2,3,  1,2,3};
            
            // numQubitsPerReg = {2, 3}
            int numOverrides = 3;
            long long int overrideInds[] = {0,0, 0,0, 0,0}; // repetition not checked
            qreal overridePhases[]       = {.1,  .1,  .1};
            
            // first element of overrideInds coordinate is a 2 qubit register
            enum bitEncoding enc = UNSIGNED;
            int badInd = GENERATE(0, 2, 4);
            overrideInds[badInd] = GENERATE( -1, (1<<2) );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, enc, coeffs, expos, numTermsPerReg, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            overrideInds[badInd] = 0;
            
            // second element of overrideInds coordinate is a 3 qubit register
            badInd += 1;
            overrideInds[badInd] = GENERATE( -1, (1<<3) );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, enc, coeffs, expos, numTermsPerReg, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            overrideInds[badInd] = 0;
            badInd -= 1;
            
            enc = TWOS_COMPLEMENT;
            int minInd = -(1<<(numQubitsPerReg[0]-1));
            int maxInd = (1<<(numQubitsPerReg[0]-1)) - 1;
            overrideInds[badInd] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, enc, coeffs, expos, numTermsPerReg, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
            overrideInds[badInd] = 0;
            
            badInd++;
            minInd = -(1<<(numQubitsPerReg[1]-1));
            maxInd = (1<<(numQubitsPerReg[1]-1)) -1;
            overrideInds[badInd] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyMultiVarPhaseFuncOverrides(quregVec, qubits, numQubitsPerReg, numRegs, enc, coeffs, expos, numTermsPerReg, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyNamedPhaseFunc
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyNamedPhaseFunc", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encoding
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every possible number of registers 
        // (between #qubits containing 1, and 1 containing #qubits)
        int numRegs;
        int maxNumRegs = 0;
        if (encoding == UNSIGNED)
            maxNumRegs = NUM_QUBITS;
        if (encoding == TWOS_COMPLEMENT)
            maxNumRegs = NUM_QUBITS/2;  // floors
        numRegs = GENERATE_COPY( range(1, maxNumRegs+1) );
        
        // try every possible total number of involved qubits
        int totalNumQubits;
        int minTotalQubits = 0;
        if (encoding == UNSIGNED)
            // each register must contain at least 1 qubit
            minTotalQubits = numRegs;
        if (encoding == TWOS_COMPLEMENT)
            // each register must contain at least 2 qubits
            minTotalQubits = 2*numRegs;
        totalNumQubits = GENERATE_COPY( range(minTotalQubits,NUM_QUBITS+1) );
                        
        // try every qubits subset and ordering 
        int* regs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), totalNumQubits) );
        
        // assign each sub-reg its minimum length
        int unallocQubits = totalNumQubits;
        int numQubitsPerReg[numRegs];
        for (int i=0; i<numRegs; i++) 
            if (encoding == UNSIGNED) {
                numQubitsPerReg[i] = 1;
                unallocQubits -= 1;
            }
            else if (encoding == TWOS_COMPLEMENT) {
                numQubitsPerReg[i] = 2;
                unallocQubits -= 2;
            }
        // and randomly allocate the remaining qubits between the registers
        while (unallocQubits > 0) {
            numQubitsPerReg[getRandomInt(0,numRegs)] += 1;
            unallocQubits--;
        }
        
        // for reference, determine the values corresponding to each register for all basis states
        qreal regVals[1<<totalNumQubits][numRegs];
        for (long long int i=0; i<(1<<totalNumQubits); i++) {
            
            long long int bits = i;
            for (int r=0; r<numRegs; r++) {            
                regVals[i][r] = bits % (1 << numQubitsPerReg[r]);
                bits = bits >> numQubitsPerReg[r];
                
                if (encoding == TWOS_COMPLEMENT)
                    regVals[i][r] = getTwosComplement(regVals[i][r], numQubitsPerReg[r]);
            }
        }
        
        /* the reference diagonal matrix which assumes the qubits are
         * contiguous and strictly increasing between the registers, and hence 
         * only depends on the number of qubits in each register.
         */
        QMatrix diagMatr = getZeroMatrix(1 << totalNumQubits);
            
        SECTION( "NORM" ) {
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
            
            SECTION( "state-vector" ) {
                
                applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, NORM);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, NORM);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 142*REAL_EPS) );
            }
        }
        SECTION( "PRODUCT" ) {
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                diagMatr[i][i] = expI(phase);
            }
            
            SECTION( "state-vector" ) {
                
                applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, PRODUCT);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, PRODUCT);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 142*REAL_EPS) );
            }
        }
        SECTION( "DISTANCE" ) {
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, DISTANCE);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, DISTANCE);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 142*REAL_EPS) );
                }
            }
        }
    }
    SECTION( "input validation" ) {
        
        int numRegs = 2;
        int numQubitsPerReg[] = {2,3};
        int regs[] = {0,1,2,3,4};
        
        SECTION( "number of registers" ) {
            
            numRegs = GENERATE_COPY( -1, 0, 1+MAX_NUM_REGS_APPLY_ARBITRARY_PHASE );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM), Contains("Invalid number of qubit subregisters") );
        }
        SECTION( "number of qubits" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0, 1+NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            regs[GENERATE(2,3,4)] = regs[1];
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {

            regs[GENERATE(range(0,NUM_QUBITS))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM), Contains("Invalid qubit index") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = 1;
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, NORM), Contains("A sub-register contained too few qubits to employ TWOS_COMPLEMENT encoding") );
        }
        SECTION( "phase function name" ) {
            
            enum phaseFunc func = (enum phaseFunc) GENERATE( -1, 14 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func), Contains("Invalid named phase function") );
        }
        SECTION( "phase function parameters" ) {

            enum phaseFunc func = GENERATE( SCALED_NORM, INVERSE_NORM, SCALED_INVERSE_NORM, SCALED_INVERSE_SHIFTED_NORM, SCALED_PRODUCT, INVERSE_PRODUCT, SCALED_INVERSE_PRODUCT, SCALED_DISTANCE, INVERSE_DISTANCE, SCALED_INVERSE_DISTANCE, SCALED_INVERSE_SHIFTED_DISTANCE );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func), Contains("Invalid number of parameters") );
        }
        SECTION( "distance pair registers" ) {
            
            int numQb[] = {1,1,1,1,1};
            int qb[] = {0,1,2,3,4};
            
            numRegs = GENERATE( 1, 3, 5 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFunc(quregVec, qb, numQb, numRegs, UNSIGNED, DISTANCE), Contains("Phase functions DISTANCE") && Contains("even number of sub-registers") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyNamedPhaseFuncOverrides
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyNamedPhaseFuncOverrides", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encoding
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every possible number of registers 
        // (between #qubits containing 1, and 1 containing #qubits)
        int numRegs;
        int maxNumRegs = 0;
        if (encoding == UNSIGNED)
            maxNumRegs = NUM_QUBITS;
        if (encoding == TWOS_COMPLEMENT)
            maxNumRegs = NUM_QUBITS/2;  // floors
        numRegs = GENERATE_COPY( range(1, maxNumRegs+1) );
        
        // try every possible total number of involved qubits
        int totalNumQubits;
        int minTotalQubits = 0;
        if (encoding == UNSIGNED)
            // each register must contain at least 1 qubit
            minTotalQubits = numRegs;
        if (encoding == TWOS_COMPLEMENT)
            // each register must contain at least 2 qubits
            minTotalQubits = 2*numRegs;
        totalNumQubits = GENERATE_COPY( range(minTotalQubits,NUM_QUBITS+1) );
                        
        // try every qubits subset and ordering 
        int* regs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), totalNumQubits) );
        
        // assign each sub-reg its minimum length
        int unallocQubits = totalNumQubits;
        int numQubitsPerReg[numRegs];
        for (int i=0; i<numRegs; i++) 
            if (encoding == UNSIGNED) {
                numQubitsPerReg[i] = 1;
                unallocQubits -= 1;
            }
            else if (encoding == TWOS_COMPLEMENT) {
                numQubitsPerReg[i] = 2;
                unallocQubits -= 2;
            }
        // and randomly allocate the remaining qubits between the registers
        while (unallocQubits > 0) {
            numQubitsPerReg[getRandomInt(0,numRegs)] += 1;
            unallocQubits--;
        }
        
        
        // choose a random number of overrides (even overriding every amplitude)
        int numOverrides = getRandomInt(0, (1<<totalNumQubits) + 1);
        
        // randomise each override index (uniqueness isn't checked)
        long long int overrideInds[numOverrides*numRegs];
        int flatInd = 0;
        for (int v=0; v<numOverrides; v++) {
            for (int r=0; r<numRegs; r++) {
                if (encoding == UNSIGNED)
                    overrideInds[flatInd] = getRandomInt(0, 1<<numQubitsPerReg[r]);
                else if (encoding == TWOS_COMPLEMENT) 
                    overrideInds[flatInd] = getRandomInt(-(1<<(numQubitsPerReg[r]-1)), (1<<(numQubitsPerReg[r]-1))-1);
                flatInd++;
            }
        }

        // override to a random phase
        qreal overridePhases[numOverrides];
        for (int v=0; v<numOverrides; v++)
            overridePhases[v] = getRandomReal(-4, 4); // periodic in [-pi, pi]
            
            
        // determine the values corresponding to each register for all basis states
        qreal regVals[1<<totalNumQubits][numRegs];
        for (long long int i=0; i<(1<<totalNumQubits); i++) {
            
            long long int bits = i;
            for (int r=0; r<numRegs; r++) {            
                regVals[i][r] = bits % (1 << numQubitsPerReg[r]);
                bits = bits >> numQubitsPerReg[r];
                
                if (encoding == TWOS_COMPLEMENT)
                    regVals[i][r] = getTwosComplement(regVals[i][r], numQubitsPerReg[r]);
            }
        }
        
        /* a reference diagonal matrix which assumes the qubits are
         * contiguous and strictly increasing between the registers, and hence 
         * only depends on the number of qubits in each register.
         */
        QMatrix diagMatr = getZeroMatrix(1 << totalNumQubits);
            
        SECTION( "NORM" ) {
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            
            SECTION( "state-vector" ) {
                
                applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, NORM, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, NORM, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E4*REAL_EPS) );
            }
        }
        SECTION( "PRODUCT" ) {
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                diagMatr[i][i] = expI(phase);
            }
            
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            
            SECTION( "state-vector" ) {
                
                applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, PRODUCT, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, PRODUCT, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E4*REAL_EPS) );
            }
        }
        SECTION( "DISTANCE" ) {
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
                
                setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, DISTANCE, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, DISTANCE, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E4*REAL_EPS) );
                }
            }
        }
    }
    SECTION( "input validation" ) {
        
        int numRegs = 2;
        int numQubitsPerReg[] = {2,3};
        int regs[] = {0,1,2,3,4};
        
        SECTION( "number of registers" ) {
            
            numRegs = GENERATE_COPY( -1, 0, 1+MAX_NUM_REGS_APPLY_ARBITRARY_PHASE );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, NULL, 0), Contains("Invalid number of qubit subregisters") );
        }
        SECTION( "number of qubits" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0, 1+NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, NULL, 0), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            regs[GENERATE(2,3,4)] = regs[1];
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, NULL, 0), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {

            regs[GENERATE(range(0,NUM_QUBITS))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, NULL, 0), Contains("Invalid qubit index") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, NULL, 0), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = 1;
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, NORM, NULL, NULL, 0), Contains("A sub-register contained too few qubits to employ TWOS_COMPLEMENT encoding") );
        }
        SECTION( "phase function name" ) {
            
            enum phaseFunc func = (enum phaseFunc) GENERATE( -1, 14 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, NULL, NULL, 0), Contains("Invalid named phase function") );
        }
        SECTION( "phase function parameters" ) {

            enum phaseFunc func = GENERATE( SCALED_NORM, INVERSE_NORM, SCALED_INVERSE_NORM, SCALED_INVERSE_SHIFTED_NORM, SCALED_PRODUCT, INVERSE_PRODUCT, SCALED_INVERSE_PRODUCT, SCALED_DISTANCE, INVERSE_DISTANCE, SCALED_INVERSE_DISTANCE, SCALED_INVERSE_SHIFTED_DISTANCE );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, NULL, NULL, 0), Contains("Invalid number of parameters") );
        }
        SECTION( "distance pair registers" ) {
            
            int numQb[] = {1,1,1,1,1};
            int qb[] = {0,1,2,3,4};
            
            numRegs = GENERATE( 1, 3, 5 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, qb, numQb, numRegs, UNSIGNED, DISTANCE, NULL, NULL, 0), Contains("Phase functions DISTANCE") && Contains("even number of sub-registers") );
        }
        SECTION( "number of overrides" ) {
   
            int numOverrides = -1;
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, NULL, numOverrides), Contains("Invalid number of phase function overrides specified") );
        }
        SECTION( "override indices" ) {
            
            // numQubitsPerReg = {2, 3}
            int numOverrides = 3;
            long long int overrideInds[] = {0,0, 0,0, 0,0}; // repetition not checked
            qreal overridePhases[]       = {.1,  .1,  .1};
            
            // first element of overrideInds coordinate is a 2 qubit register
            enum bitEncoding enc = UNSIGNED;
            int badInd = GENERATE(0, 2, 4);
            overrideInds[badInd] = GENERATE( -1, (1<<2) );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            overrideInds[badInd] = 0;
            
            // second element of overrideInds coordinate is a 3 qubit register
            badInd += 1;
            overrideInds[badInd] = GENERATE( -1, (1<<3) );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            overrideInds[badInd] = 0;
            badInd -= 1;
            
            enc = TWOS_COMPLEMENT;
            int minInd = -(1<<(numQubitsPerReg[0]-1));
            int maxInd = (1<<(numQubitsPerReg[0]-1)) - 1;
            overrideInds[badInd] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
            overrideInds[badInd] = 0;
            
            badInd++;
            minInd = -(1<<(numQubitsPerReg[1]-1));
            maxInd = (1<<(numQubitsPerReg[1]-1)) -1;
            overrideInds[badInd] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyParamNamedPhaseFunc
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyParamNamedPhaseFunc", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encoding
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every possible number of registers 
        // (between #qubits containing 1, and 1 containing #qubits)
        int numRegs;
        int maxNumRegs = 0;
        if (encoding == UNSIGNED)
            maxNumRegs = NUM_QUBITS;
        if (encoding == TWOS_COMPLEMENT)
            maxNumRegs = NUM_QUBITS/2;  // floors
        numRegs = GENERATE_COPY( range(1, maxNumRegs+1) );
        
        // try every possible total number of involved qubits
        int totalNumQubits;
        int minTotalQubits = 0;
        if (encoding == UNSIGNED)
            // each register must contain at least 1 qubit
            minTotalQubits = numRegs;
        if (encoding == TWOS_COMPLEMENT)
            // each register must contain at least 2 qubits
            minTotalQubits = 2*numRegs;
        totalNumQubits = GENERATE_COPY( range(minTotalQubits,NUM_QUBITS+1) );
                        
        // try every qubits subset and ordering 
        int* regs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), totalNumQubits) );
        
        // assign each sub-reg its minimum length
        int unallocQubits = totalNumQubits;
        int numQubitsPerReg[numRegs];
        for (int i=0; i<numRegs; i++) 
            if (encoding == UNSIGNED) {
                numQubitsPerReg[i] = 1;
                unallocQubits -= 1;
            }
            else if (encoding == TWOS_COMPLEMENT) {
                numQubitsPerReg[i] = 2;
                unallocQubits -= 2;
            }
        // and randomly allocate the remaining qubits between the registers
        while (unallocQubits > 0) {
            numQubitsPerReg[getRandomInt(0,numRegs)] += 1;
            unallocQubits--;
        }
            
        /* produce a reference diagonal matrix which assumes the qubits are
         * contiguous and strictly increasing between the registers, and hence 
         * only depends on the number of qubits in each register.
         */
        QMatrix diagMatr = getZeroMatrix(1 << totalNumQubits);
        
        // determine the values corresponding to each register for all basis states
        qreal regVals[1<<totalNumQubits][numRegs];
        for (long long int i=0; i<(1<<totalNumQubits); i++) {
            
            long long int bits = i;
            for (int r=0; r<numRegs; r++) {            
                regVals[i][r] = bits % (1 << numQubitsPerReg[r]);
                bits = bits >> numQubitsPerReg[r];
                
                if (encoding == TWOS_COMPLEMENT)
                    regVals[i][r] = getTwosComplement(regVals[i][r], numQubitsPerReg[r]);
            }
        }
        
        SECTION( "INVERSE_NORM" ) {
            
            enum phaseFunc func = INVERSE_NORM;
            qreal divPhase = getRandomReal(-4, 4);
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = (phase == 0.)? divPhase : 1/sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
                            
            qreal params[] = {divPhase};
            int numParams = 1;
            
            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }

        }
        SECTION( "SCALED_NORM" ) {
            
            enum phaseFunc func = SCALED_NORM;
            qreal coeff = getRandomReal(-10, 10);
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = coeff * sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
                            
            qreal params[] = {coeff};
            int numParams = 1;

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_INVERSE_NORM" ) {
            
            enum phaseFunc func = SCALED_INVERSE_NORM;
            qreal coeff = getRandomReal(-10, 10);
            qreal divPhase = getRandomReal(-4, 4);
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = (phase == 0.)? divPhase : coeff/sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
                            
            qreal params[] = {coeff, divPhase};
            int numParams = 2;

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_INVERSE_SHIFTED_NORM" ) {
            
            enum phaseFunc func = SCALED_INVERSE_SHIFTED_NORM;
            int numParams = 2 + numRegs;
            qreal params[numParams];
            params[0] = getRandomReal(-10, 10); // scaling
            params[1] = getRandomReal(-4, 4); // divergence override
            for (int r=0; r<numRegs; r++)
                params[2+r] = getRandomReal(-8, 8); // shifts
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r] - params[2+r], 2);
                phase = sqrt(phase);
                phase = (phase <= REAL_EPS)? params[1] : params[0]/phase;
                diagMatr[i][i] = expI(phase);
            }
            
            SECTION( "state-vector" ) {

                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "INVERSE_PRODUCT" ) {
            
            enum phaseFunc func = INVERSE_PRODUCT;
            qreal divPhase = getRandomReal(-4, 4);
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                phase = (phase == 0.)? divPhase : 1. / phase;
                diagMatr[i][i] = expI(phase);
            }
                            
            qreal params[] = {divPhase};
            int numParams = 1;

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_PRODUCT" ) {
            
            enum phaseFunc func = SCALED_PRODUCT;
            qreal coeff = getRandomReal(-10, 10);
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                diagMatr[i][i] = expI(coeff * phase);
            }
                            
            qreal params[] = {coeff};
            int numParams = 1;

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_INVERSE_PRODUCT" ) {
            
            enum phaseFunc func = SCALED_INVERSE_PRODUCT;
            qreal coeff = getRandomReal(-10, 10);
            qreal divPhase = getRandomReal(-4, 4);
            qreal params[] = {coeff, divPhase};
            int numParams = 2;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                phase = (phase == 0)? divPhase : coeff / phase;
                diagMatr[i][i] = expI(phase);
            }

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "INVERSE_DISTANCE" ) {
            
            enum phaseFunc func = INVERSE_DISTANCE;
            qreal divPhase = getRandomReal( -4, 4 );
            qreal params[] = {divPhase};
            int numParams = 1;
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = (phase == 0.)? divPhase : 1./sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
        SECTION( "SCALED_DISTANCE" ) {
            
            enum phaseFunc func = SCALED_DISTANCE;
            qreal coeff = getRandomReal( -10, 10 );
            qreal params[] = {coeff};
            int numParams = 1;
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = coeff * sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
        SECTION( "SCALED_INVERSE_DISTANCE" ) {
            
            enum phaseFunc func = SCALED_INVERSE_DISTANCE;
            qreal coeff = getRandomReal( -10, 10 );
            qreal divPhase = getRandomReal( -4, 4 );
            qreal params[] = {coeff, divPhase};
            int numParams = 2;
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = (phase == 0.)? divPhase : coeff/sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
        SECTION( "SCALED_INVERSE_SHIFTED_DISTANCE" ) {
            
            enum phaseFunc func = SCALED_INVERSE_SHIFTED_DISTANCE;
            int numParams = 2 + numRegs/2;
            qreal params[numParams];
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {

                params[0] = getRandomReal( -10, 10 ); // scaling
                params[1] = getRandomReal( -4, 4 ); // divergence override
                for (int r=0; r<numRegs/2; r++)
                    params[2+r] = getRandomReal( -8, 8 ); // shifts
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r]-regVals[i][r+1]-params[2+r/2], 2);
                    phase = sqrt(phase);
                    phase = (phase <= REAL_EPS)? params[1] : params[0]/phase;
                    diagMatr[i][i] = expI(phase);
                }
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFunc(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
    }
    SECTION( "input validation" ) {
        
        int numRegs = 2;
        int numQubitsPerReg[] = {2,3};
        int regs[] = {0,1,2,3,4};
        
        SECTION( "number of registers" ) {
            
            numRegs = GENERATE_COPY( -1, 0, 1+MAX_NUM_REGS_APPLY_ARBITRARY_PHASE );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0), Contains("Invalid number of qubit subregisters") );
        }
        SECTION( "number of qubits" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0, 1+NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            regs[GENERATE(2,3,4)] = regs[1];
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {

            regs[GENERATE(range(0,NUM_QUBITS))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0), Contains("Invalid qubit index") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, 0), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = 1;
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, NORM, NULL, 0), Contains("A sub-register contained too few qubits to employ TWOS_COMPLEMENT encoding") );
        }
        SECTION( "phase function name" ) {
            
            enum phaseFunc func = (enum phaseFunc) GENERATE( -1, 14 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, NULL, 0), Contains("Invalid named phase function") );
        }
        SECTION( "phase function parameters" ) {
            
            qreal params[numRegs + 3];
            for (int r=0; r<numRegs + 3; r++)
                params[r] = 0;
            
            SECTION( "no parameter functions" ) {
                
                enum phaseFunc func = GENERATE( NORM, PRODUCT, DISTANCE  );
                int numParams = GENERATE( -1, 1, 2 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams), Contains("Invalid number of parameters") );
            }
            SECTION( "single parameter functions" ) {
                
                enum phaseFunc func = GENERATE( SCALED_NORM, INVERSE_NORM, SCALED_PRODUCT, INVERSE_PRODUCT, SCALED_DISTANCE, INVERSE_DISTANCE );
                int numParams = GENERATE( -1, 0, 2 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams), Contains("Invalid number of parameters") );
            }
            SECTION( "two parameter functions" ) {    
                
                enum phaseFunc func = GENERATE( SCALED_INVERSE_NORM, SCALED_INVERSE_PRODUCT, SCALED_INVERSE_DISTANCE );
                int numParams = GENERATE( 0, 1, 3 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams), Contains("Invalid number of parameters") );
            }
            SECTION( "shifted distance" ) {
                
                if (numRegs%2 == 0) {
                    enum phaseFunc func = SCALED_INVERSE_SHIFTED_DISTANCE;
                    int numParams = GENERATE_COPY( 0, 1, numRegs/2 - 1, numRegs/2, numRegs/2 + 1, numRegs/2 + 3 );
                    REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams), Contains("Invalid number of parameters") );
                }
            }
            SECTION( "shifted norm" ) {
                
                enum phaseFunc func = SCALED_INVERSE_SHIFTED_NORM;
                int numParams = GENERATE_COPY( 0, 1, numRegs-1, numRegs, numRegs+1, numRegs+3 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams), Contains("Invalid number of parameters") );
            }
        }
        SECTION( "distance pair registers" ) {
            
            int numQb[] = {1,1,1,1,1};
            int qb[] = {0,1,2,3,4};
            
            numRegs = GENERATE( 1, 3, 5 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFunc(quregVec, qb, numQb, numRegs, UNSIGNED, DISTANCE, NULL, 0), Contains("Phase functions DISTANCE") && Contains("even number of sub-registers") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyParamNamedPhaseFuncOverrides
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyParamNamedPhaseFuncOverrides", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encoding
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every possible number of registers 
        // (between #qubits containing 1, and 1 containing #qubits)
        int numRegs;
        int maxNumRegs = 0;
        if (encoding == UNSIGNED)
            maxNumRegs = NUM_QUBITS;
        if (encoding == TWOS_COMPLEMENT)
            maxNumRegs = NUM_QUBITS/2;  // floors
        numRegs = GENERATE_COPY( range(1, maxNumRegs+1) );
        
        // try every possible total number of involved qubits
        int totalNumQubits;
        int minTotalQubits = 0;
        if (encoding == UNSIGNED)
            // each register must contain at least 1 qubit
            minTotalQubits = numRegs;
        if (encoding == TWOS_COMPLEMENT)
            // each register must contain at least 2 qubits
            minTotalQubits = 2*numRegs;
        totalNumQubits = GENERATE_COPY( range(minTotalQubits,NUM_QUBITS+1) );
                        
        // try every qubits subset and ordering 
        int* regs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), totalNumQubits) );
        
        // assign each sub-reg its minimum length
        int unallocQubits = totalNumQubits;
        int numQubitsPerReg[numRegs];
        for (int i=0; i<numRegs; i++) 
            if (encoding == UNSIGNED) {
                numQubitsPerReg[i] = 1;
                unallocQubits -= 1;
            }
            else if (encoding == TWOS_COMPLEMENT) {
                numQubitsPerReg[i] = 2;
                unallocQubits -= 2;
            }
        // and randomly allocate the remaining qubits between the registers
        while (unallocQubits > 0) {
            numQubitsPerReg[getRandomInt(0,numRegs)] += 1;
            unallocQubits--;
        }
        
        
        // choose a random number of overrides (even overriding every amplitude)
        int numOverrides = getRandomInt(0, (1<<totalNumQubits) + 1);
        
        // randomise each override index (uniqueness isn't checked)
        long long int overrideInds[numOverrides*numRegs];
        int flatInd = 0;
        for (int v=0; v<numOverrides; v++) {
            for (int r=0; r<numRegs; r++) {
                if (encoding == UNSIGNED)
                    overrideInds[flatInd] = getRandomInt(0, 1<<numQubitsPerReg[r]);
                else if (encoding == TWOS_COMPLEMENT) 
                    overrideInds[flatInd] = getRandomInt(-(1<<(numQubitsPerReg[r]-1)), (1<<(numQubitsPerReg[r]-1))-1);
                flatInd++;
            }
        }

        // override to a random phase
        qreal overridePhases[numOverrides];
        for (int v=0; v<numOverrides; v++)
            overridePhases[v] = getRandomReal(-4, 4); // periodic in [-pi, pi]


        // determine the values corresponding to each register for all basis states
        qreal regVals[1<<totalNumQubits][numRegs];
        for (long long int i=0; i<(1<<totalNumQubits); i++) {
            
            long long int bits = i;
            for (int r=0; r<numRegs; r++) {            
                regVals[i][r] = bits % (1 << numQubitsPerReg[r]);
                bits = bits >> numQubitsPerReg[r];
                
                if (encoding == TWOS_COMPLEMENT)
                    regVals[i][r] = getTwosComplement(regVals[i][r], numQubitsPerReg[r]);
            }
        }
        
        /* produce a reference diagonal matrix which assumes the qubits are
         * contiguous and strictly increasing between the registers, and hence 
         * only depends on the number of qubits in each register.
         */
        QMatrix diagMatr = getZeroMatrix(1 << totalNumQubits);
        
        
        SECTION( "INVERSE_NORM" ) {
            
            enum phaseFunc func = INVERSE_NORM;
            qreal divPhase = getRandomReal(-4, 4);
            qreal params[] = {divPhase};
            int numParams = 1;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = (phase == 0.)? divPhase : 1/sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            
            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_NORM" ) {
            
            enum phaseFunc func = SCALED_NORM;
            qreal coeff = getRandomReal(-10, 10);
            qreal params[] = {coeff};
            int numParams = 1;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = coeff * sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_INVERSE_NORM" ) {
            
            enum phaseFunc func = SCALED_INVERSE_NORM;
            qreal coeff = getRandomReal(-10, 10);
            qreal divPhase = getRandomReal(-4, 4);
            qreal params[] = {coeff, divPhase};
            int numParams = 2;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r], 2);
                phase = (phase == 0.)? divPhase : coeff/sqrt(phase);
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            
            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_INVERSE_SHIFTED_NORM" ) {
            
            enum phaseFunc func = SCALED_INVERSE_SHIFTED_NORM;
            int numParams = 2 + numRegs;
            qreal params[numParams];
            params[0] = getRandomReal(-10, 10); // scaling
            params[1] = getRandomReal(-4, 4); // divergence override
            for (int r=0; r<numRegs; r++)
                params[2+r] = getRandomReal(-8, 8); // shifts
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 0;
                for (int r=0; r<numRegs; r++)
                    phase += pow(regVals[i][r] - params[2+r], 2);
                phase = sqrt(phase);
                phase = (phase <= REAL_EPS)? params[1] : params[0]/phase;
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }        
        SECTION( "INVERSE_PRODUCT" ) {
            
            enum phaseFunc func = INVERSE_PRODUCT;
            qreal divPhase = getRandomReal(-4, 4);
            qreal params[] = {divPhase};
            int numParams = 1;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                phase = (phase == 0.)? divPhase : 1. / phase;
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);

            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_PRODUCT" ) {
            
            enum phaseFunc func = SCALED_PRODUCT;
            qreal coeff = getRandomReal(-10, 10);
            qreal params[] = {coeff};
            int numParams = 1;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                diagMatr[i][i] = expI(coeff * phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            
            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "SCALED_INVERSE_PRODUCT" ) {
            
            enum phaseFunc func = SCALED_INVERSE_PRODUCT;
            qreal coeff = getRandomReal(-10, 10);
            qreal divPhase = getRandomReal(-4, 4);
            qreal params[] = {coeff, divPhase};
            int numParams = 2;
            
            for (size_t i=0; i<diagMatr.size(); i++) {
                qreal phase = 1;
                for (int r=0; r<numRegs; r++)
                    phase *= regVals[i][r];
                phase = (phase == 0)? divPhase : coeff / phase;
                diagMatr[i][i] = expI(phase);
            }
            setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            
            SECTION( "state-vector" ) {
                
                applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
            }
            SECTION( "density-matrix" ) {
                
                applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
            }
        }
        SECTION( "INVERSE_DISTANCE" ) {
            
            enum phaseFunc func = INVERSE_DISTANCE;
            qreal divPhase = getRandomReal( -4, 4 );
            qreal params[] = {divPhase};
            int numParams = 1;
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = (phase == 0.)? divPhase : 1./sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
                setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
        SECTION( "SCALED_DISTANCE" ) {
            
            enum phaseFunc func = SCALED_DISTANCE;
            qreal coeff = getRandomReal( -10, 10 );                
            qreal params[] = {coeff};
            int numParams = 1;
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = coeff * sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
                setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
        SECTION( "SCALED_INVERSE_DISTANCE" ) {
            
            enum phaseFunc func = SCALED_INVERSE_DISTANCE;
            qreal coeff = getRandomReal( -10, 10 );
            qreal divPhase = getRandomReal( -4, 4 );
            qreal params[] = {coeff, divPhase};
            int numParams = 2;
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r+1]-regVals[i][r], 2);
                    phase = (phase == 0.)? divPhase : coeff/sqrt(phase);
                    diagMatr[i][i] = expI(phase);
                }
                setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
        SECTION( "SCALED_INVERSE_SHIFTED_DISTANCE" ) {
            
            enum phaseFunc func = SCALED_INVERSE_SHIFTED_DISTANCE;
            int numParams = 2 + numRegs/2;
            qreal params[numParams];
            
            // test only if there are an even number of registers
            if (numRegs%2 == 0) {

                params[0] = getRandomReal( -10, 10 ); // scaling
                params[1] = getRandomReal( -4, 4 ); // divergence override
                for (int r=0; r<numRegs/2; r++)
                    params[2+r] = getRandomReal( -8, 8 ); // shifts
                
                for (size_t i=0; i<diagMatr.size(); i++) {
                    qreal phase = 0;
                    for (int r=0; r<numRegs; r+=2)
                        phase += pow(regVals[i][r]-regVals[i][r+1]-params[2+r/2], 2);
                    phase = sqrt(phase);
                    phase = (phase <= REAL_EPS)? params[1] : params[0]/phase;
                    diagMatr[i][i] = expI(phase);
                }
                
                setDiagMatrixOverrides(diagMatr, numQubitsPerReg, numRegs, encoding, overrideInds, overridePhases, numOverrides);
            }
            
            SECTION( "state-vector" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refVec, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregVec, refVec, 1E2*REAL_EPS) );
                }
            }
            SECTION( "density-matrix" ) {
                
                if (numRegs%2 == 0) {
                    applyParamNamedPhaseFuncOverrides(quregMatr, regs, numQubitsPerReg, numRegs, encoding, func, params, numParams, overrideInds, overridePhases, numOverrides);
                    applyReferenceOp(refMatr, regs, totalNumQubits, diagMatr);
                    REQUIRE( areEqual(quregMatr, refMatr, 1E2*REAL_EPS) );
                }
            }
        }
    }
    SECTION( "input validation" ) {
        
        int numRegs = 2;
        int numQubitsPerReg[] = {2,3};
        int regs[] = {0,1,2,3,4};
        
        SECTION( "number of registers" ) {
            
            numRegs = GENERATE_COPY( -1, 0, 1+MAX_NUM_REGS_APPLY_ARBITRARY_PHASE );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0, NULL, NULL, 0), Contains("Invalid number of qubit subregisters") );
        }
        SECTION( "number of qubits" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = GENERATE( -1, 0, 1+NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0, NULL, NULL, 0), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            regs[GENERATE(2,3,4)] = regs[1];
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0, NULL, NULL, 0), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {

            regs[GENERATE(range(0,NUM_QUBITS))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0, NULL, NULL, 0), Contains("Invalid qubit index") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, 0, NULL, NULL, 0), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubitsPerReg[GENERATE_COPY(range(0,numRegs))] = 1;
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, TWOS_COMPLEMENT, NORM, NULL, 0, NULL, NULL, 0), Contains("A sub-register contained too few qubits to employ TWOS_COMPLEMENT encoding") );
        }
        SECTION( "phase function name" ) {
            
            enum phaseFunc func = (enum phaseFunc) GENERATE( -1, 14 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, NULL, 0, NULL, NULL, 0), Contains("Invalid named phase function") );
        }
        SECTION( "phase function parameters" ) {
            
            qreal params[numRegs + 3];
            for (int r=0; r<numRegs + 3; r++)
                params[r] = 0;
            
            SECTION( "no parameter functions" ) {
                
                enum phaseFunc func = GENERATE( NORM, PRODUCT, DISTANCE  );
                int numParams = GENERATE( -1, 1, 2 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams, NULL, NULL, 0), Contains("Invalid number of parameters") );
            }
            SECTION( "single parameter functions" ) {
                
                enum phaseFunc func = GENERATE( SCALED_NORM, INVERSE_NORM, SCALED_PRODUCT, INVERSE_PRODUCT, SCALED_DISTANCE, INVERSE_DISTANCE );
                int numParams = GENERATE( -1, 0, 2 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams, NULL, NULL, 0), Contains("Invalid number of parameters") );
            }
            SECTION( "two parameter functions" ) {    
                
                enum phaseFunc func = GENERATE( SCALED_INVERSE_NORM, SCALED_INVERSE_PRODUCT, SCALED_INVERSE_DISTANCE );
                int numParams = GENERATE( 0, 1, 3 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams, NULL, NULL, 0), Contains("Invalid number of parameters") );
            }
            SECTION( "shifted distance" ) {
                
                if (numRegs%2 == 0) {
                    enum phaseFunc func = SCALED_INVERSE_SHIFTED_DISTANCE;
                    int numParams = GENERATE_COPY( 0, 1, numRegs/2 - 1, numRegs/2, numRegs/2 + 1, numRegs/2 + 3 );
                    REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams, NULL, NULL, 0), Contains("Invalid number of parameters") );
                }
            }
            SECTION( "shifted norm" ) {
                
                enum phaseFunc func = SCALED_INVERSE_SHIFTED_NORM;
                int numParams = GENERATE_COPY( 0, 1, numRegs-1, numRegs, numRegs+1, numRegs+3 );
                REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, func, params, numParams, NULL, NULL, 0), Contains("Invalid number of parameters") );
            }
        }
        SECTION( "distance pair registers" ) {
            
            int numQb[] = {1,1,1,1,1};
            int qb[] = {0,1,2,3,4};
            
            numRegs = GENERATE( 1, 3, 5 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, qb, numQb, numRegs, UNSIGNED, DISTANCE, NULL, 0, NULL, NULL, 0), Contains("Phase functions DISTANCE") && Contains("even number of sub-registers") );
        }
        SECTION( "number of overrides" ) {
   
            int numOverrides = -1;
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, UNSIGNED, NORM, NULL, 0, NULL, NULL, numOverrides), Contains("Invalid number of phase function overrides specified") );
        }
        SECTION( "override indices" ) {
            
            // numQubitsPerReg = {2, 3}
            int numOverrides = 3;
            long long int overrideInds[] = {0,0, 0,0, 0,0}; // repetition not checked
            qreal overridePhases[]       = {.1,  .1,  .1};
            
            // first element of overrideInds coordinate is a 2 qubit register
            enum bitEncoding enc = UNSIGNED;
            int badInd = GENERATE(0, 2, 4);
            overrideInds[badInd] = GENERATE( -1, (1<<2) );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, 0, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            overrideInds[badInd] = 0;
            
            // second element of overrideInds coordinate is a 3 qubit register
            badInd += 1;
            overrideInds[badInd] = GENERATE( -1, (1<<3) );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, 0, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            overrideInds[badInd] = 0;
            badInd -= 1;
            
            enc = TWOS_COMPLEMENT;
            int minInd = -(1<<(numQubitsPerReg[0]-1));
            int maxInd = (1<<(numQubitsPerReg[0]-1)) - 1;
            overrideInds[badInd] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, 0, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
            overrideInds[badInd] = 0;
            
            badInd++;
            minInd = -(1<<(numQubitsPerReg[1]-1));
            maxInd = (1<<(numQubitsPerReg[1]-1)) -1;
            overrideInds[badInd] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyParamNamedPhaseFuncOverrides(quregVec, regs, numQubitsPerReg, numRegs, enc, NORM, NULL, 0, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
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



/** @sa applyPhaseFunc
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyPhaseFunc", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encodings
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every sub-register size
        int numQubits = GENERATE_COPY( range(1,NUM_QUBITS+1) );
        
        // force at least 2 qubits in two's compement though
        if (encoding == TWOS_COMPLEMENT && numQubits == 1)
            numQubits++;
        
        // try every possible sub-register
        int* qubits = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numQubits) );
        
        // choose a random number of terms in the phase function
        int numTerms = getRandomInt(1, 5);
        
        // populate the phase function with random but POSITIVE-power terms,
        // and in two's complement mode, strictly integer powers
        qreal coeffs[numTerms];
        qreal expons[numTerms];
        for (int t=0; t<numTerms; t++) {
            coeffs[t] = getRandomReal(-10,10);
            if (encoding == TWOS_COMPLEMENT)
                expons[t] = getRandomInt(0, 3+1);  // repetition of power is ok
            else 
                expons[t] = getRandomReal(0, 3);
        }
        
        // build a reference diagonal matrix, on the reduced Hilbert space
        QMatrix matr = getZeroMatrix( 1 << numQubits );
        for (size_t i=0; i<matr.size(); i++) {
            
            long long int ind = 0;
            if (encoding == UNSIGNED)
                ind = i;
            if (encoding == TWOS_COMPLEMENT)
                ind = getTwosComplement(i, numQubits);
            
            qreal phase = 0;
            for (int t=0; t<numTerms; t++)
                phase += coeffs[t] * pow(ind, expons[t]);
                
            matr[i][i] = expI(phase);
        }
        
        SECTION( "state-vector" ) {
            
            applyPhaseFunc(quregVec, qubits, numQubits, encoding, coeffs, expons, numTerms);
            applyReferenceOp(refVec, qubits, numQubits, matr);
            REQUIRE( areEqual(quregVec, refVec, 1E4*REAL_EPS) );
        }
        SECTION( "density-matrix" ) {
        
            applyPhaseFunc(quregMatr, qubits, numQubits, encoding, coeffs, expons, numTerms);
            applyReferenceOp(refMatr, qubits, numQubits, matr);
            REQUIRE( areEqual(quregMatr, refMatr, 1E6*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        int numQubits = 3;
        int qubits[] = {0,1,2};
        
        SECTION( "number of qubits" ) {
            
            numQubits = GENERATE_COPY( -1, 0, NUM_QUBITS+1 );
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, NULL, numQubits, UNSIGNED, NULL, NULL, 1), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition of qubits" ) {
            
            qubits[GENERATE(1,2)] = qubits[0];
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, UNSIGNED, NULL, NULL, 1), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {
            
            int inv = GENERATE( -1, NUM_QUBITS );
            qubits[ GENERATE_COPY( range(0,numQubits) )] = inv;
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, UNSIGNED, NULL, NULL, 1), Contains("Invalid qubit index") );
        }
        SECTION( "number of terms" ) {
            
            int numTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, UNSIGNED, NULL, NULL, numTerms), Contains("Invalid number of terms in the phase function") );
        }
        SECTION( "bit encoding name" ) {
            
            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1,2 );
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, enc, NULL, NULL, 1), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubits = 1;
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, TWOS_COMPLEMENT, NULL, NULL, 1), Contains("too few qubits to employ TWOS_COMPLEMENT") );
        }
        SECTION( "fractional exponent" ) {
            
            int numTerms = 3;
            qreal coeffs[] = {0,0,0};
            qreal expos[] = {1,2,3};
            expos[GENERATE_COPY( range(0,numTerms) )] = GENERATE( 0.5, 1.999, 5.0001 );
            
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, TWOS_COMPLEMENT, coeffs, expos, numTerms), Contains("fractional exponent") && Contains("TWOS_COMPLEMENT") && Contains("negative indices were not overriden") );
        }
        SECTION( "negative exponent" ) {

            int numTerms = 3;
            qreal coeffs[] = {0,0,0};
            qreal expos[] = {1,2,3};
            expos[GENERATE_COPY( range(0,numTerms) )] = GENERATE( -1, -2, -1.5 );
            
            enum bitEncoding encoding = GENERATE( UNSIGNED, TWOS_COMPLEMENT );
    
            REQUIRE_THROWS_WITH( applyPhaseFunc(quregVec, qubits, numQubits, encoding, coeffs, expos, numTerms), Contains("The phase function contained a negative exponent which would diverge at zero, but the zero index was not overriden") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyPhaseFuncOverrides
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyPhaseFuncOverrides", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );

    SECTION( "correctness" ) {
        
        // try every kind of binary encodings
        enum bitEncoding encoding = GENERATE( UNSIGNED,TWOS_COMPLEMENT );
        
        // try every sub-register size
        int numQubits = GENERATE_COPY( range(1,NUM_QUBITS+1) );
        
        // force at least 2 qubits in two's compement though
        if (encoding == TWOS_COMPLEMENT && numQubits == 1)
            numQubits++;
        
        // try every possible sub-register
        int* qubits = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numQubits) );
        
        // choose a random number of terms in the phase function
        int numTerms = getRandomInt(1, 21);
        
        // populate the phase function with random powers.
        // this includes negative powers, so we must always override 0.
        // in two's complement, we must use only integer powers
        qreal coeffs[numTerms];
        qreal expons[numTerms];
        for (int t=0; t<numTerms; t++) {
            coeffs[t] = getRandomReal(-10,10);
            if (encoding == TWOS_COMPLEMENT)
                expons[t] = getRandomInt(-3, 3+1); // note we COULD do getRandomReal(), and override all negatives
            else
                expons[t] = getRandomReal(-3, 3);
        }
        
        // choose a random number of overrides (even overriding every amplitude)
        int numOverrides = getRandomInt(1, (1<<numQubits) + 1);
        
        // randomise each override index (uniqueness isn't checked)
        long long int overrideInds[numOverrides];
        overrideInds[0] = 0LL;
        for (int i=1; i<numOverrides; i++)
            if (encoding == UNSIGNED)
                overrideInds[i] = getRandomInt(0, 1<<numQubits);
            else if (encoding == TWOS_COMPLEMENT)
                overrideInds[i] = getRandomInt(-(1<<(numQubits-1)), (1<<(numQubits-1))-1);
            
        // override to a random phase
        qreal overridePhases[numOverrides];
        for (int i=0; i<numOverrides; i++)
            overridePhases[i] = getRandomReal(-4, 4); // periodic in [-pi, pi]
            
        // build a reference diagonal matrix, on the reduced Hilbert space
        QMatrix matr = getZeroMatrix( 1 << numQubits );
        for (size_t i=0; i<matr.size(); i++) {
            
            long long int ind = 0;
            if (encoding == UNSIGNED)
                ind = i;
            if (encoding == TWOS_COMPLEMENT)
                ind = getTwosComplement(i, numQubits);
                
            // reference diagonal matrix incorporates overriden phases
            qreal phase;
            bool overriden = false;
            for (int v=0; v<numOverrides; v++) {
                if (ind == overrideInds[v]) {
                    phase = overridePhases[v];
                    overriden = true;
                    break;
                }
            }
            
            if (!overriden) {
                phase = 0;
                for (int t=0; t<numTerms; t++)
                    phase += coeffs[t] * pow(ind, expons[t]);
            }
                
            matr[i][i] = expI(phase);
        }         
        
        SECTION( "state-vector" ) {
            
            applyPhaseFuncOverrides(quregVec, qubits, numQubits, encoding, coeffs, expons, numTerms, overrideInds, overridePhases, numOverrides);
            applyReferenceOp(refVec, qubits, numQubits, matr);
            REQUIRE( areEqual(quregVec, refVec, 1E4*REAL_EPS) );
        }
        SECTION( "density-matrix" ) {
            
            applyPhaseFuncOverrides(quregMatr, qubits, numQubits, encoding, coeffs, expons, numTerms, overrideInds, overridePhases, numOverrides);
            applyReferenceOp(refMatr, qubits, numQubits, matr);
            REQUIRE( areEqual(quregMatr, refMatr, 1E6*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        int numQubits = 3;
        int qubits[] = {0,1,2};

        SECTION( "number of qubits" ) {
            
            int numQubits = GENERATE_COPY( -1, 0, NUM_QUBITS+1 );
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, NULL, numQubits, UNSIGNED, NULL, NULL, 1, NULL, NULL, 0), Contains("Invalid number of qubits") );
        }
        SECTION( "repetition qubits" ) {
            
            qubits[GENERATE(1,2)] = qubits[0];
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, UNSIGNED, NULL, NULL, 1, NULL, NULL, 0), Contains("The qubits must be unique") );
        }
        SECTION( "qubit indices" ) {
            
            int inv = GENERATE( -1, NUM_QUBITS );
            qubits[ GENERATE_COPY( range(0,numQubits) )] = inv;
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, UNSIGNED, NULL, NULL, 1, NULL, NULL, 0), Contains("Invalid qubit index") );
        }
        SECTION( "number of terms" ) {
            
            int numTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, UNSIGNED, NULL, NULL, numTerms, NULL, NULL, 0), Contains("Invalid number of terms in the phase function") );
        }
        SECTION( "bit encoding name" ) {

            enum bitEncoding enc = (enum bitEncoding) GENERATE( -1,2 );
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, enc, NULL, NULL, 1, NULL, NULL, 0), Contains("Invalid bit encoding") );
        }
        SECTION( "two's complement register" ) {
            
            numQubits = 1;
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, TWOS_COMPLEMENT, NULL, NULL, 1, NULL, NULL, 0), Contains("too few qubits to employ TWOS_COMPLEMENT") );
        }
        SECTION( "number of overrides" ) {

            qreal dummyTerms[] = {0};
            
            int numOverrides = GENERATE_COPY( -1, 1 + (1<<numQubits) );
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, UNSIGNED, dummyTerms, dummyTerms, 1, NULL, NULL, numOverrides), Contains("Invalid number of phase function overrides") );
        }
        SECTION( "override indices" ) {
            
            int numOverrides = 3;
            long long int overrideInds[] = {0,1,2};
            qreal overridePhases[] = {.1,.1,.1};
            qreal dummyTerms[] = {0};
            
            enum bitEncoding encoding = UNSIGNED;
            overrideInds[GENERATE(0,1,2)] = GENERATE_COPY( -1, (1<<numQubits) );
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, encoding, dummyTerms, dummyTerms, 1, overrideInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the UNSIGNED encoding") );
            
            encoding = TWOS_COMPLEMENT;
            long long int newInds[] = {0,1,2};
            int minInd = -(1<<(numQubits-1));
            int maxInd = (1<<(numQubits-1)) -1;
            newInds[GENERATE(0,1,2)] = GENERATE_COPY( minInd-1, maxInd+1 );
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, encoding, dummyTerms, dummyTerms, 1, newInds, overridePhases, numOverrides), Contains("Invalid phase function override index, in the TWOS_COMPLEMENT encoding") );
        }
        SECTION( "fractional exponent" ) {
            
            int numTerms = 3;
            qreal coeffs[] = {0,0,0};
            qreal expos[] = {1,2,3};
            
            // make one exponent fractional, thereby requiring negative overrides
            expos[GENERATE_COPY( range(0,numTerms) )] = GENERATE( 0.5, 1.999, 5.0001 );
            
            // catch when no negative indices are overridden
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, TWOS_COMPLEMENT, coeffs, expos, numTerms, NULL, NULL, 0), Contains("fractional exponent") && Contains("TWOS_COMPLEMENT") && Contains("negative indices were not overriden") );
            
            int numNegs = 1 << (numQubits-1);
            long long int overrideInds[numNegs];
            qreal overridePhases[numNegs];
            for (int i=0; i<numNegs; i++) {
                overrideInds[i] = -(i+1);
                overridePhases[i] = 0;
            }
            
            // ensure no throw when all are overriden
            REQUIRE_NOTHROW( applyPhaseFuncOverrides(quregVec, qubits, numQubits, TWOS_COMPLEMENT, coeffs, expos, numTerms, overrideInds, overridePhases, numNegs) );

            // catch when at least one isn't overriden
            overrideInds[GENERATE_COPY( range(0,numNegs) )] = 0; // override a non-negative
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, TWOS_COMPLEMENT, coeffs, expos, numTerms, overrideInds, overridePhases, numNegs), Contains("fractional exponent") && Contains("TWOS_COMPLEMENT") && Contains("negative indices were not overriden") );
        }
        SECTION( "negative exponent" ) {
            
            int numTerms = 3;
            qreal coeffs[] = {0,0,0};
            qreal expos[] = {1,2,3};
            expos[GENERATE_COPY( range(0,numTerms) )] = GENERATE( -1, -2 );
            
            enum bitEncoding encoding = GENERATE( UNSIGNED, TWOS_COMPLEMENT );
            
            // test both when giving no overrides, and giving all non-zero overrides
            int numOverrides = GENERATE( 0, 3 ); 
            long long int overrideInds[] = {1,2,3};
            qreal overridePhases[] = {0,0,0};
            REQUIRE_THROWS_WITH( applyPhaseFuncOverrides(quregVec, qubits, numQubits, encoding, coeffs, expos, numTerms, overrideInds, overridePhases, numOverrides), Contains("The phase function contained a negative exponent which would diverge at zero, but the zero index was not overriden") );
            
            // but ensure that when the zero IS overriden (anywhere), there's no error 
            numOverrides = 3;
            overrideInds[GENERATE_COPY(range(0,numOverrides))] = 0;
            REQUIRE_NOTHROW( applyPhaseFuncOverrides(quregVec, qubits, numQubits, encoding, coeffs, expos, numTerms, overrideInds, overridePhases, numOverrides) );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
}



/** @sa applyProjector
 * @ingroup unittest 
 * @author Tyson Jones 
 */
TEST_CASE( "applyProjector", "[operators]" ) {
    
    Qureg vec = createQureg(NUM_QUBITS, QUEST_ENV);
    Qureg mat = createDensityQureg(NUM_QUBITS, QUEST_ENV);
    
    SECTION( "correctness" ) {
        
        int qubit = GENERATE( range(0,NUM_QUBITS) );
        int outcome = GENERATE( 0, 1 );
        
        // repeat these random tests 10 times on every qubit, and for both outcomes
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            SECTION( "normalised" ) {
                
                // use a random L2 state for every qubit & outcome
                QVector vecRef = getRandomStateVector(NUM_QUBITS);
                toQureg(vec, vecRef);
                
                // zero non-outcome reference amps
                for (size_t ind=0; ind<vecRef.size(); ind++) {
                    int bit = (ind >> qubit) & 1; // target-th bit
                    if (bit != outcome)
                        vecRef[ind] = 0;
                }
                
                applyProjector(vec, qubit, outcome);
                REQUIRE( areEqual(vec, vecRef) );
            }
            SECTION( "unnormalised" ) {
                
                // use a random non-physical state for every qubit & outcome
                QVector vecRef = getRandomQVector(1 << NUM_QUBITS);
                toQureg(vec, vecRef);
                
                // zero non-outcome reference amps
                for (size_t ind=0; ind<vecRef.size(); ind++) {
                    int bit = (ind >> qubit) & 1; // target-th bit
                    if (bit != outcome)
                        vecRef[ind] = 0;
                }
                
                applyProjector(vec, qubit, outcome);
                REQUIRE( areEqual(vec, vecRef) );
            }
        }        
        SECTION( "density-matrix" ) {
            
            SECTION( "pure" ) {

                QVector vecRef = getRandomStateVector(NUM_QUBITS);
                QMatrix matRef = getPureDensityMatrix(vecRef);
                
                toQureg(mat, matRef);
                applyProjector(mat, qubit, outcome);
                
                // zero any amplitudes that aren't |outcome><outcome|
                for (size_t r=0; r<matRef.size(); r++) {
                    for (size_t c=0; c<matRef.size(); c++) {
                        int ketBit = (c >> qubit) & 1;
                        int braBit = (r >> qubit) & 1;
                        if (!(ketBit == outcome && braBit == outcome))
                            matRef[r][c] = 0;
                    }
                }
                
                REQUIRE( areEqual(mat, matRef) );
            }
            SECTION( "mixed" ) {
                
                QMatrix matRef = getRandomDensityMatrix(NUM_QUBITS);
                
                toQureg(mat, matRef);
                applyProjector(mat, qubit, outcome);
                
                // zero any amplitudes that aren't |outcome><outcome|
                for (size_t r=0; r<matRef.size(); r++) {
                    for (size_t c=0; c<matRef.size(); c++) {
                        int ketBit = (c >> qubit) & 1;
                        int braBit = (r >> qubit) & 1;
                        if (!(ketBit == outcome && braBit == outcome))
                            matRef[r][c] = 0;
                    }
                }
                
                REQUIRE( areEqual(mat, matRef) );
            }
            SECTION( "unnormalised" ) {
                
                QMatrix matRef = getRandomQMatrix(1 << NUM_QUBITS);
                
                toQureg(mat, matRef);
                applyProjector(mat, qubit, outcome);
                
                // zero any amplitudes that aren't |outcome><outcome|
                for (size_t r=0; r<matRef.size(); r++) {
                    for (size_t c=0; c<matRef.size(); c++) {
                        int ketBit = (c >> qubit) & 1;
                        int braBit = (r >> qubit) & 1;
                        if (!(ketBit == outcome && braBit == outcome))
                            matRef[r][c] = 0;
                    }
                }
                
                REQUIRE( areEqual(mat, matRef) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
            int qubit = GENERATE( -1, NUM_QUBITS );
            int outcome = 0;
            REQUIRE_THROWS_WITH( applyProjector(mat, qubit, outcome), Contains("Invalid target qubit") );
        }
        SECTION( "outcome value" ) {
            
            int qubit = 0;
            int outcome = GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( applyProjector(mat, qubit, outcome), Contains("Invalid measurement outcome") );
        }
    }
    destroyQureg(vec, QUEST_ENV);
    destroyQureg(mat, QUEST_ENV);
}



/** @sa applyQFT
 * @ingroup unittest
 * @author Tyson Jones 
 */
TEST_CASE( "applyQFT", "[operators]" ) {
    
    PREPARE_TEST( quregVec, quregMatr, refVec, refMatr );
    
    SECTION( "correctness" ) {
        
        // try every sub-register size
        int numQubits = GENERATE_COPY( range(1,NUM_QUBITS+1) );
        
        // try every possible sub-register
        int* qubits = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numQubits) );
        
        SECTION( "state-vector" ) {
            
            SECTION( "normalised" ) {
        
                QVector refVec = getRandomStateVector(NUM_QUBITS);
                toQureg(quregVec, refVec);
                
                applyQFT(quregVec, qubits, numQubits);
                refVec = getDFT(refVec, qubits, numQubits);
                
                REQUIRE( areEqual(quregVec, refVec) );
            }
            SECTION( "unnormalised" ) {
                
                QVector refVec = getRandomQVector(1 << NUM_QUBITS);
                toQureg(quregVec, refVec);
                
                applyQFT(quregVec, qubits, numQubits);
                refVec = getDFT(refVec, qubits, numQubits);
                
                REQUIRE( areEqual(quregVec, refVec) );
            }
        }
        SECTION( "density-matrix" ) {
            
            SECTION( "pure" ) {
                
                /* a pure density matrix should be mapped to a pure state 
                 * corresponding to the state-vector DFT
                 */
                
                refVec = getRandomStateVector(NUM_QUBITS);
                refMatr = getPureDensityMatrix(refVec);
                toQureg(quregMatr, refMatr);
                
                applyQFT(quregMatr, qubits, numQubits);
                refVec = getDFT(refVec, qubits, numQubits);
                refMatr = getPureDensityMatrix(refVec);
                
                REQUIRE( areEqual(quregMatr, refMatr) );
            }
            SECTION( "mixed" ) {
                
                /* a mixed density matrix, conceptualised as a mixture of orthogonal
                 * state-vectors, should be mapped to an equally weighted mixture 
                 * of DFTs of each state-vector (because QFT is unitary and hence 
                 * maintains state orthogonality)
                 */
                
                int numStates = (1 << NUM_QUBITS)/4; // a quarter of as many states as are possible
                std::vector<QVector> states = getRandomOrthonormalVectors(NUM_QUBITS, numStates);
                std::vector<qreal> probs = getRandomProbabilities(numStates);
                
                // set qureg to random mixture of state-vectors
                refMatr = getMixedDensityMatrix(probs, states);
                toQureg(quregMatr, refMatr);
                
                // apply QFT to mixture
                applyQFT(quregMatr, qubits, numQubits);
                
                // compute dft of mixture, via dft of each state
                refMatr = getZeroMatrix(1 << NUM_QUBITS);
                for (int i=0; i<numStates; i++) {
                    QVector dft = getDFT(states[i], qubits, numQubits);
                    refMatr += probs[i] * getPureDensityMatrix(dft);
                }
                
                REQUIRE( areEqual(quregMatr, refMatr) );
            }
            SECTION( "unnormalised" ) {
                
                /* repeat method above, except that we use unnormalised vectors, 
                 * and mix them with arbitrary complex numbers instead of probabilities,
                 * yielding an unnormalised density matrix 
                 */
                
                int numVecs = (1 << NUM_QUBITS)/4; // a quarter of as many states as are possible
                std::vector<QVector> vecs;
                std::vector<qcomp> coeffs;
                for (int i=0; i<numVecs; i++) {
                    vecs.push_back(getRandomQVector(1 << NUM_QUBITS));
                    coeffs.push_back(getRandomComplex());
                }
                
                // produce unnormalised matrix via random complex sum of random unnormalised vectors
                refMatr = getZeroMatrix(1 << NUM_QUBITS);
                for (int i=0; i<numVecs; i++)
                    refMatr += coeffs[i] * getPureDensityMatrix(vecs[i]);
                    
                toQureg(quregMatr, refMatr);
                applyQFT(quregMatr, qubits, numQubits);
                
                // compute target matrix via dft of each unnormalised vector 
                refMatr = getZeroMatrix(1 << NUM_QUBITS);
                for (int i=0; i<numVecs; i++) {
                    QVector dft = getDFT(vecs[i], qubits, numQubits);
                    refMatr += coeffs[i] * getPureDensityMatrix(dft);
                }
                    
                REQUIRE( areEqual(quregMatr, refMatr) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "number of targets" ) {
            
            // there cannot be more targets than qubits in register
            int numQubits = GENERATE( -1, 0, NUM_QUBITS+1 );
            int qubits[NUM_QUBITS+1];
            
            REQUIRE_THROWS_WITH( applyQFT(quregVec, qubits, numQubits), Contains("Invalid number of target"));
        }
        SECTION( "repetition in targets" ) {
            
            int numQubits = 3;
            int qubits[] = {1,2,2};
            
            REQUIRE_THROWS_WITH( applyQFT(quregVec, qubits, numQubits), Contains("target") && Contains("unique"));
        }
        SECTION( "qubit indices" ) {
            
            int numQubits = 3;
            int qubits[] = {1,2,3};
            
            int inv = GENERATE( -1, NUM_QUBITS );
            qubits[GENERATE_COPY( range(0,numQubits) )] = inv; // make invalid target
            REQUIRE_THROWS_WITH( applyQFT(quregVec, qubits, numQubits), Contains("Invalid target") );
        }
    }
    CLEANUP_TEST( quregVec, quregMatr );
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
                REQUIRE( areEqual(vec, vecRef, 10*REAL_EPS) );
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
            
            PauliHamil hamil = createPauliHamil(NUM_QUBITS, 3);
            qreal coeffs[] = {(qreal) (M_PI * sqrt(2.0)), M_PI, M_PI};
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

