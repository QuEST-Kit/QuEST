/** @file
 * Ported tests of the deprecated QuEST v3 interface,
 * unit testing the "calculations" module.
 * 
 * This file should be excluded from doxygen parsing so 
 * as not to conflict with the doc of the v4 unit tests.
 * 
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

// must define preprocessors to enable quest's
// deprecated v3 API, and disable the numerous
// warnings issued by its compilation
#define INCLUDE_DEPRECATED_FUNCTIONS 1
#define DISABLE_DEPRECATION_WARNINGS 1
#include "quest/include/quest.h"

#include "test_utilities.hpp"

/* allows concise use of ContainsSubstring in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::ContainsSubstring;
using Catch::Approx;


/** @sa calcDensityInnerProduct
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcDensityInnerProduct", "[calculations]" ) {

    Qureg mat1 = createForcedDensityQureg(NUM_QUBITS);
    Qureg mat2 = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        // repeat these random tests 10 times 
        GENERATE( range(0,10) );
        
        SECTION( "density-matrix" ) {
            
            SECTION( "pure" ) {
                
                // mat1 = |r1><r1|
                QVector r1 = getRandomStateVector(NUM_QUBITS);
                toQureg(mat1, getKetBra(r1,r1));
                
                // mat2 = |r2><r2|
                QVector r2 = getRandomStateVector(NUM_QUBITS);
                toQureg(mat2, getKetBra(r2,r2));
                
                // prod( |r1><r1|, |r2><r2| ) = |<r1|r2>|^2 
                qcomp prod = 0;
                for (size_t i=0; i<r1.size(); i++)
                    prod += conj(r1[i]) * r2[i];
                qreal densProd = pow(abs(prod),2);
                
                REQUIRE( real(calcDensityInnerProduct(mat1,mat2)) == Approx(densProd).margin(100 * REAL_EPS) );
            }
            SECTION( "mixed" ) {
                
                QMatrix ref1 = getRandomDensityMatrix(NUM_QUBITS);
                QMatrix ref2 = getRandomDensityMatrix(NUM_QUBITS);
                toQureg(mat1, ref1);
                toQureg(mat2, ref2);
                
                // prod(mat1, mat2) = sum_{ij} conj(mat1_{ij}) * mat2_{ij}
                qcomp refProd = 0;
                for (size_t i=0; i<ref1.size(); i++)
                    for (size_t j=0; j<ref1.size(); j++)
                        refProd += conj(ref1[i][j]) * ref2[i][j];
                REQUIRE( imag(refProd) == Approx(0).margin(REAL_EPS) );
                
                REQUIRE( real(calcDensityInnerProduct(mat1,mat2)) == Approx(real(refProd)).margin(100 * REAL_EPS) );
                
                // should be invariant under ordering
                REQUIRE( real(calcDensityInnerProduct(mat1,mat2)) == Approx(real(calcDensityInnerProduct(mat2,mat1))).margin(100 * REAL_EPS) );
            }
            SECTION( "unnormalised" ) {
                
                // set both to random (non-Hermitian) complex matrices
                QMatrix ref1 = getRandomQMatrix(1<<NUM_QUBITS);
                QMatrix ref2 = getRandomQMatrix(1<<NUM_QUBITS);
                toQureg(mat1, ref1);
                toQureg(mat2, ref2);
                
                // prod(mat1, mat2) = real(sum_{ij} conj(mat1_{ij}) * mat2_{ij})
                qcomp refProd = 0;
                for (size_t i=0; i<ref1.size(); i++)
                    for (size_t j=0; j<ref1.size(); j++)
                        refProd += conj(ref1[i][j]) * ref2[i][j];
                        
                REQUIRE( real(calcDensityInnerProduct(mat1,mat2)) == Approx(real(refProd)).margin(100 * REAL_EPS) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "dimensions" ) {
            
            Qureg mat3 = createDensityQureg(NUM_QUBITS + 1);
            REQUIRE_THROWS_WITH( calcDensityInnerProduct(mat1,mat3), ContainsSubstring("differing numbers of qubits") );
            destroyQureg(mat3);
        }

        // in v4, calcDensityInnerProduct() redirects to calcInnerProduct() which is
        // valid for both statevectors and density matrices (and combinations thereof!)

            // SECTION( "state-vectors" ) {
                
            //     Qureg vec = createForcedQureg(NUM_QUBITS);
                
            //     REQUIRE_THROWS_WITH( calcDensityInnerProduct(mat1,vec), ContainsSubstring("valid only for density matrices") );
            //     REQUIRE_THROWS_WITH( calcDensityInnerProduct(vec,mat1), ContainsSubstring("valid only for density matrices") );
            //     REQUIRE_THROWS_WITH( calcDensityInnerProduct(vec,vec),  ContainsSubstring("valid only for density matrices") );        
                
            //     destroyQureg(vec);
            // }
    }
    destroyQureg(mat1);
    destroyQureg(mat2);
}



/** @sa calcExpecDiagonalOp
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcExpecDiagonalOp", "[calculations]" ) {

    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    initDebugState(vec);
    initDebugState(mat);
    QVector vecRef = toQVector(vec);
    QMatrix matRef = toQMatrix(mat);
    
    SECTION( "correctness" ) {
        
        // try 10 random operators 
        GENERATE( range(0,10) );
        
        // make a totally random (non-Hermitian) diagonal oeprator
        DiagonalOp op = createDiagonalOp(NUM_QUBITS, getQuESTEnv());
        for (long long int i=0; i<op.numElemsPerNode; i++)
            op.cpuElems[i] = getRandomComplex();
        syncDiagonalOp(op);
        
        SECTION( "state-vector" ) {

            /* calcExpecDiagOp calculates <qureg|diag|qureg> */
            
            QVector sumRef = toQMatrix(op) * vecRef;
            qcomp prod = 0;
            for (size_t i=0; i<vecRef.size(); i++)
                prod += conj(vecRef[i]) * sumRef[i];
            
            qcomp res = calcExpecDiagonalOp(vec, op);
            REQUIRE( real(res) == Approx(real(prod)).margin(REAL_EPS) );
            REQUIRE( imag(res) == Approx(imag(prod)).margin(REAL_EPS) );
        } 
        SECTION( "density-matrix" ) {
            
            /* calcExpecDiagOp calculates Trace( diag * qureg ) */
            matRef = toQMatrix(op) * matRef;            
            qcomp tr = 0;
            for (size_t i=0; i<matRef.size(); i++)
                tr += matRef[i][i];

            qcomp res = calcExpecDiagonalOp(mat, op);
            REQUIRE( real(res) == Approx(real(tr)).margin(100*REAL_EPS) );
            REQUIRE( imag(res) == Approx(imag(tr)).margin(100*REAL_EPS) );
        }
        
        destroyDiagonalOp(op, getQuESTEnv());
    }
    SECTION( "input validation" ) {
        
        SECTION( "mismatching size" ) {
            
            DiagonalOp op = createDiagonalOp(NUM_QUBITS + 1, getQuESTEnv());
            
            REQUIRE_THROWS_WITH( calcExpecDiagonalOp(vec, op), ContainsSubstring("different number of qubits"));
            REQUIRE_THROWS_WITH( calcExpecDiagonalOp(mat, op), ContainsSubstring("different number of qubits"));
            
            destroyDiagonalOp(op, getQuESTEnv());
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



// calcExpecPauliHamil removed because PauliHamil is deprecated,
// and replacement PauliStrSum presently has no compatible constructor

    // /** @sa calcExpecPauliHamil
    //  * @ingroup deprecatedtests 
    //  * @author Tyson Jones 
    //  */
    // TEST_CASE( "calcExpecPauliHamil", "[calculations]" ) {
        
    //     Qureg vec = createForcedQureg(NUM_QUBITS);
    //     Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    //     initDebugState(vec);
    //     initDebugState(mat);
    //     QVector vecRef = toQVector(vec);
    //     QMatrix matRef = toQMatrix(mat);
        
    //     Qureg vecWork = createForcedQureg(NUM_QUBITS);
    //     Qureg matWork = createForcedDensityQureg(NUM_QUBITS);
        
    //     SECTION( "correctness" ) {
            
    //         /* it's too expensive to try every possible Pauli configuration, so
    //          * we'll try 10 random codes, and for each, random coefficients
    //          */
    //         GENERATE( range(0,10) );
                    
    //         int numTerms = GENERATE( 1, 2, 10, 15 );
    //         PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);
    //         setRandomPauliSum(hamil);
    //         QMatrix refHamil = toQMatrix(hamil);
            
    //         SECTION( "state-vector" ) {

    //             /* calcExpecPauliHamil calculates <qureg|pauliHum|qureg> */
                
    //             QVector sumRef = refHamil * vecRef;
    //             qcomp prod = 0;
    //             for (size_t i=0; i<vecRef.size(); i++)
    //                 prod += conj(vecRef[i]) * sumRef[i];
    //             REQUIRE( imag(prod) == Approx(0).margin(10*REAL_EPS) );
                
    //             qreal res = calcExpecPauliHamil(vec, hamil, vecWork);
    //             REQUIRE( res == Approx(real(prod)).margin(10*REAL_EPS) );
    //         } 
    //         SECTION( "density-matrix" ) {
                
    //             /* calcExpecPauliHamil calculates Trace( pauliHamil * qureg ) */
    //             matRef = refHamil * matRef;            
    //             qreal tr = 0;
    //             for (size_t i=0; i<matRef.size(); i++)
    //                 tr += real(matRef[i][i]);
    //             // (get real, since we start in a non-Hermitian state, hence diagonal isn't real)
                
    //             qreal res = calcExpecPauliHamil(mat, hamil, matWork);
    //             REQUIRE( res == Approx(tr).margin(1E2*REAL_EPS) );
    //         }
            
    //         destroyPauliHamil(hamil);
    //     }
    //     SECTION( "input validation" ) {
            
    //         SECTION( "pauli codes" ) {
                
    //             int numTerms = 3;
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);

    //             // make one pauli code wrong
    //             hamil.pauliCodes[GENERATE_COPY( range(0,numTerms*NUM_QUBITS) )] = (pauliOpType) GENERATE( -1, 4 );
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(vec, hamil, vecWork), ContainsSubstring("Invalid Pauli code") );
                
    //             destroyPauliHamil(hamil);
    //         }
    //         SECTION( "workspace type" ) {
                
    //             int numTerms = 1;
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);
                
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(vec, hamil, mat), ContainsSubstring("Registers must both be state-vectors or both be density matrices") );
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(mat, hamil, vec), ContainsSubstring("Registers must both be state-vectors or both be density matrices") );
                
    //             destroyPauliHamil(hamil);
    //         }
    //         SECTION( "workspace dimensions" ) {
                    
    //             int numTerms = 1;
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS, numTerms);
        
    //             Qureg vec2 = createQureg(NUM_QUBITS + 1);
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(vec, hamil, vec2), ContainsSubstring("Dimensions") && ContainsSubstring("don't match") );
    //             destroyQureg(vec2);
                
    //             Qureg mat2 = createDensityQureg(NUM_QUBITS + 1);
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(mat, hamil, mat2), ContainsSubstring("Dimensions") && ContainsSubstring("don't match") );
    //             destroyQureg(mat2);
                
    //             destroyPauliHamil(hamil);
    //         }
    //         SECTION( "matching hamiltonian qubits" ) {
                
    //             int numTerms = 1;
    //             PauliHamil hamil = createPauliHamil(NUM_QUBITS + 1, numTerms);
                
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(vec, hamil, vecWork), ContainsSubstring("same number of qubits") );
    //             REQUIRE_THROWS_WITH( calcExpecPauliHamil(mat, hamil, matWork), ContainsSubstring("same number of qubits") );
                
    //             destroyPauliHamil(hamil);
    //         }
    //     }
    //     destroyQureg(vec);
    //     destroyQureg(mat);
    //     destroyQureg(vecWork);
    //     destroyQureg(matWork);
    // }



/** @sa calcExpecPauliProd
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcExpecPauliProd", "[calculations]" ) {
    
    QuESTEnv env = getQuESTEnv();
    Qureg vec = createCustomQureg(NUM_QUBITS, 0, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded);
    Qureg mat = createCustomQureg(NUM_QUBITS, 1, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded);

    initDebugState(vec);
    initDebugState(mat);
    QVector vecRef = toQVector(vec);
    QMatrix matRef = toQMatrix(mat);
    
    Qureg vecWork = createForcedQureg(NUM_QUBITS);
    Qureg matWork = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        int numTargs = GENERATE( range(1,NUM_QUBITS+1) );
        int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        
        /* it's too expensive to try ALL Pauli sequences, via 
         *      pauliOpType* paulis = GENERATE_COPY( pauliseqs(numTargs) );.
         * Furthermore, take(10, pauliseqs(numTargs)) will try the same pauli codes.
         * Hence, we instead opt to repeatedlyrandomly generate pauliseqs
         */
        GENERATE( range(0,10) ); // gen 10 random pauli-codes for every targs
        vector<pauliOpType> paulis(numTargs);
        for (int i=0; i<numTargs; i++)
            paulis[i] = (pauliOpType) getRandomInt(0,4);

        // produce a numTargs-big matrix 'pauliProd' by pauli-matrix tensoring
        QMatrix iMatr{{1,0},{0,1}};
        QMatrix xMatr{{0,1},{1,0}};
        QMatrix yMatr{{0,-qcomp(0,1)},{qcomp(0,1),0}};
        QMatrix zMatr{{1,0},{0,-1}};
        QMatrix pauliProd{{1}};
        for (int i=0; i<numTargs; i++) {
            QMatrix fac;
            if (paulis[i] == PAULI_I) fac = iMatr;
            if (paulis[i] == PAULI_X) fac = xMatr;
            if (paulis[i] == PAULI_Y) fac = yMatr;
            if (paulis[i] == PAULI_Z) fac = zMatr;
            pauliProd = getKroneckerProduct(fac, pauliProd);
        }
        
        SECTION( "state-vector" ) {

            /* calcExpecPauliProd calculates <qureg|pauliProd|qureg> */
            
            QVector prodRef = vecRef; 
            applyReferenceOp(prodRef, targs, numTargs, pauliProd);
            qcomp prod = 0;
            for (size_t i=0; i<vecRef.size(); i++)
                prod += conj(vecRef[i]) * prodRef[i];
            REQUIRE( imag(prod) == Approx(0).margin(REAL_EPS) );
            
            qreal res = calcExpecPauliProd(vec, targs, paulis.data(), numTargs, vecWork);
            REQUIRE( res == Approx(real(prod)).margin(REAL_EPS) );
        }
        SECTION( "density-matrix" ) {
            
            /* calcExpecPauliProd calculates Trace( pauliProd * qureg ) */

            // produce (pauliProd * mat)
            QMatrix fullOp = getFullOperatorMatrix(NULL, 0, targs, numTargs, pauliProd, NUM_QUBITS);
            matRef = fullOp * matRef;
            
            // compute real(trace(pauliProd * mat))
            qreal tr = 0;
            for (size_t i=0; i<matRef.size(); i++)
                tr += real(matRef[i][i]);
            // (get real, since we start in a non-Hermitian state, hence diagonal isn't real)

            // disable validation during call, because result is non-real and will upset post-check
            setValidationOff();
            qreal res = calcExpecPauliProd(mat, targs, paulis.data(), numTargs, matWork);
            setValidationOn();

            REQUIRE( res == Approx(tr).margin(10*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {

        pauliOpType pauliOpsToAvoidDeprecSegFault[100];
        for (int i=0; i<100; i++)
            pauliOpsToAvoidDeprecSegFault[i] = PAULI_X;
        
        SECTION( "number of targets" ) {
            
            int targs[NUM_QUBITS+1];
            for (int i=0; i<NUM_QUBITS+1; i++)
                targs[i] = i;

            // too few
            int numTargs = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, pauliOpsToAvoidDeprecSegFault, numTargs, vecWork), ContainsSubstring("Invalid number of Paulis") );

            // too many
            numTargs = NUM_QUBITS + 1;
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, pauliOpsToAvoidDeprecSegFault, numTargs, vecWork), ContainsSubstring("highest-index non-identity Pauli operator") && ContainsSubstring("exceeds the maximum target") );
        }
        SECTION( "target indices" ) {
            
            int numTargs = 3;
            int targs[3] = {0, 1, 2};

            int badInd = GENERATE( range(0,3) );
            
            // make one index too small
            targs[badInd] = -1;
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, pauliOpsToAvoidDeprecSegFault, numTargs, vecWork), ContainsSubstring("Pauli indices must be non-negative") );

            // make one index too big
            targs[badInd] = NUM_QUBITS;
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, pauliOpsToAvoidDeprecSegFault, numTargs, vecWork), ContainsSubstring("highest-index non-identity Pauli operator") );

            // make one index WAY too big
            targs[badInd] = 65;
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, pauliOpsToAvoidDeprecSegFault, numTargs, vecWork), ContainsSubstring("exceed the maximum number of representable Pauli operators") );
        }
        SECTION( "repetition in targets" ) {
            
            int numTargs = 3;
            int targs[3] = {0, 1, 1};
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, pauliOpsToAvoidDeprecSegFault, numTargs, vecWork), ContainsSubstring("Indices must be unique") );
        }
        SECTION( "pauli codes" ) {
            
            int numTargs = 3;
            int targs[3] = {0, 1, 2};
            pauliOpType codes[3] = {PAULI_X, PAULI_Y, PAULI_Z};
            
            // make one pauli wrong
            codes[GENERATE( range(0,3) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, codes, numTargs, vecWork), ContainsSubstring("invalid Pauli code") );
        }

        // workspace is no longer needed; argument is ignored

            // SECTION( "workspace type" ) {
                
            //     int numTargs = 1;
            //     int targs[1] = {0};
            //     pauliOpType codes[1] = {PAULI_I};
                
            //     REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, codes, numTargs, matWork), ContainsSubstring("Registers must both be state-vectors or both be density matrices") );
            //     REQUIRE_THROWS_WITH( calcExpecPauliProd(mat, targs, codes, numTargs, vecWork), ContainsSubstring("Registers must both be state-vectors or both be density matrices") );
            // }
            // SECTION( "workspace dimensions" ) {
                    
            //     int numTargs = 1;
            //     int targs[1] = {0};
            //     pauliOpType codes[1] = {PAULI_I};
        
            //     Qureg vec2 = createQureg(NUM_QUBITS + 1);
            //     REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, codes, numTargs, vec2), ContainsSubstring("Dimensions") && ContainsSubstring("don't match") );
            //     destroyQureg(vec2);
                
            //     Qureg mat2 = createDensityQureg(NUM_QUBITS + 1);
            //     REQUIRE_THROWS_WITH( calcExpecPauliProd(mat, targs, codes, numTargs, mat2), ContainsSubstring("Dimensions") && ContainsSubstring("don't match") );
            //     destroyQureg(mat2);
            // }
    }
    destroyQureg(vec);
    destroyQureg(mat);
    destroyQureg(vecWork);
    destroyQureg(matWork);
}



/** @sa calcExpecPauliSum
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcExpecPauliSum", "[calculations]" ) {
   
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);

    QVector vecRef = getRandomStateVector(NUM_QUBITS);
    QMatrix matRef = getRandomDensityMatrix(NUM_QUBITS);
    toQureg(vec, vecRef);
    toQureg(mat, matRef);
    
    // accepted by v3 deprecated API but discarded by v4
    Qureg vecWork = createForcedQureg(NUM_QUBITS);
    Qureg matWork = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        int numSumTerms = GENERATE( 1, 2, 10, 15 );
        
        /* it's too expensive to try every possible Pauli configuration, so
         * we'll try 10 random codes, and for each, random coefficients
         */
        GENERATE( range(0,10) );
        int totNumCodes = numSumTerms * NUM_QUBITS;
        vector<pauliOpType> paulis(totNumCodes);
        vector<qreal> coeffs(numSumTerms);
        setRandomPauliSum(coeffs.data(), paulis.data(), NUM_QUBITS, numSumTerms);
        
        // produce a numTargs-big matrix 'pauliSum' by pauli-matrix tensoring and summing
        QMatrix pauliSum = toQMatrix(coeffs.data(), paulis.data(), NUM_QUBITS, numSumTerms);
        
        SECTION( "state-vector" ) {

            /* calcExpecPauliSum calculates <qureg|pauliSum|qureg> */
            
            QVector sumRef = pauliSum * vecRef;
            qcomp prod = 0;
            for (size_t i=0; i<vecRef.size(); i++)
                prod += conj(vecRef[i]) * sumRef[i];
            REQUIRE( imag(prod) == Approx(0).margin(10*REAL_EPS) );
            
            qreal res = calcExpecPauliSum(vec, paulis.data(), coeffs.data(), numSumTerms, vecWork);
            REQUIRE( res == Approx(real(prod)).margin(10*REAL_EPS) );
        } 
        SECTION( "density-matrix" ) {
            
            /* calcExpecPauliSum calculates Trace( pauliSum * qureg ) */
            matRef = pauliSum * matRef;            
            qreal tr = 0;
            for (size_t i=0; i<matRef.size(); i++)
                tr += real(matRef[i][i]);
            // (get real, since we start in a non-Hermitian state, hence diagonal isn't real)
            
            qreal res = calcExpecPauliSum(mat, paulis.data(), coeffs.data(), numSumTerms, matWork);
            REQUIRE( res == Approx(tr).margin(1E2*REAL_EPS) );
        }
    }
    SECTION( "input validation" ) {
        
        // cannot be validated; deprecated API copies before validating numSumTerms, causing segfault

            // SECTION( "number of sum terms" ) {
                
            //     int numSumTerms = GENERATE( -1, 0 );
            //     REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, NULL, NULL, numSumTerms, vecWork), ContainsSubstring("The number of terms must be a positive integer") );
            // }

        SECTION( "pauli codes" ) {
            
            // make valid params
            int numSumTerms = 3;
            vector<qreal> coeffs(numSumTerms);
            vector<pauliOpType> codes(numSumTerms*NUM_QUBITS);
            for (int i=0; i<numSumTerms*NUM_QUBITS; i++)
                codes[i] = PAULI_I;

            // make one pauli wrong
            codes[GENERATE_COPY( range(0,numSumTerms*NUM_QUBITS) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, codes.data(), coeffs.data(), numSumTerms, vecWork), ContainsSubstring("invalid Pauli code") );
        }

        // the v4 API does not use a workspace, so it is discarded by the v3 deprecation layer

            // SECTION( "workspace type" ) {
                
            //     // make valid params
            //     int numSumTerms = 1;
            //     qreal coeffs[1] = {0};
            //     pauliOpType codes[NUM_QUBITS];
            //     for (int i=0; i<NUM_QUBITS; i++)
            //         codes[i] = PAULI_I;
                
            //     REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, codes, coeffs, numSumTerms, mat), ContainsSubstring("Registers must both be state-vectors or both be density matrices") );
            //     REQUIRE_THROWS_WITH( calcExpecPauliSum(mat, codes, coeffs, numSumTerms, vec), ContainsSubstring("Registers must both be state-vectors or both be density matrices") );
            // }
            
            // SECTION( "workspace dimensions" ) {
                    
            //     // make valid params
            //     int numSumTerms = 1;
            //     qreal coeffs[1] = {0};
            //     pauliOpType codes[NUM_QUBITS];
            //     for (int i=0; i<NUM_QUBITS; i++)
            //         codes[i] = PAULI_I;
        
            //     Qureg vec2 = createQureg(NUM_QUBITS + 1);
            //     REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, codes, coeffs, numSumTerms, vec2), ContainsSubstring("Dimensions") && ContainsSubstring("don't match") );
            //     destroyQureg(vec2);
                
            //     Qureg mat2 = createDensityQureg(NUM_QUBITS + 1);
            //     REQUIRE_THROWS_WITH( calcExpecPauliSum(mat, codes, coeffs, numSumTerms, mat2), ContainsSubstring("Dimensions") && ContainsSubstring("don't match") );
            //     destroyQureg(mat2);
            // }
    }
    destroyQureg(vec);
    destroyQureg(mat);
    destroyQureg(vecWork);
    destroyQureg(matWork);
}



/** @sa calcFidelity
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcFidelity", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    Qureg pure = createForcedQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        // repeat the below random tests 10 times 
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            /* calcFidelity computes |<vec|pure>|^2 */
             
            SECTION( "normalised" ) {
                 
                // random L2 vectors
                QVector vecRef = getRandomStateVector(NUM_QUBITS);
                QVector pureRef = getRandomStateVector(NUM_QUBITS);
                toQureg(vec, vecRef);
                toQureg(pure, pureRef);
                
                // |<vec|vec>|^2 = |1|^2 = 1
                REQUIRE( calcFidelity(vec,vec) == Approx(1) );
                
                // |<vec|pure>|^2 = |sum_j conj(vec_j) * pure_j|^2
                qcomp dotProd = 0;
                for (size_t i=0; i<vecRef.size(); i++)
                    dotProd += conj(vecRef[i]) * pureRef[i];
                qreal refFid = pow(abs(dotProd), 2);
                
                REQUIRE( calcFidelity(vec,pure) == Approx(refFid) );
            }

            // unnormalised test is no longer supported, since v4 calcFidelity
            // validates thet the fidelity is correctly approximately real

                // SECTION( "unnormalised" ) {
                
                //     // random unnormalised vectors
                //     QVector vecRef = getRandomQVector(1<<NUM_QUBITS);
                //     QVector pureRef = getRandomQVector(1<<NUM_QUBITS);
                //     toQureg(vec, vecRef);
                //     toQureg(pure, pureRef);
                    
                //     // Let nv be magnitude of vec, hence |unit-vec> = 1/sqrt(nv)|vec>
                //     qreal nv = 0;
                //     for (size_t i=0; i<vecRef.size(); i++)
                //         nv += pow(abs(vecRef[i]), 2);
                //     // then <vec|vec> = sqrt(nv)*sqrt(nv) <unit-vec|unit-vec> = nv,
                //     // hence |<vec|vec>|^2 = nv*nv
                //     REQUIRE( calcFidelity(vec,vec) == Approx( nv*nv ) );
                    
                //     qcomp dotProd = 0;
                //     for (size_t i=0; i<vecRef.size(); i++)
                //         dotProd += conj(vecRef[i]) * pureRef[i];
                //     qreal refFid = pow(abs(dotProd), 2);
                    
                //     REQUIRE( calcFidelity(vec,pure) == Approx(refFid) ); 
                // }
        }
        SECTION( "density-matrix" ) {
            
            /* calcFidelity computes <pure|mat|pure> */
            
            SECTION( "pure" ) {
                
                QVector pureRef = getRandomStateVector(NUM_QUBITS);
                toQureg(pure, pureRef);
                
                // test when density matrix is the same pure state 
                QMatrix matRef = getKetBra(pureRef, pureRef);
                toQureg(mat, matRef);
                REQUIRE( calcFidelity(mat,pure) == Approx(1) ); 
                
                // test when density matrix is a random pure state
                QVector r1 = getRandomStateVector(NUM_QUBITS);
                matRef = getKetBra(r1, r1); // actually pure |r1><r1|
                toQureg(mat, matRef);
                
                // <pure|r1><r1|pure> = |<r1|pure>|^2 = |sum_j conj(r1_j) * pure_j|^2
                qcomp dotProd = 0;
                for (size_t i=0; i<r1.size(); i++)
                    dotProd += conj(r1[i]) * pureRef[i];
                qreal refFid = pow(abs(dotProd), 2);
                
                REQUIRE( calcFidelity(mat,pure) == Approx(refFid).margin(100 * REAL_EPS) ); 
            }
            SECTION( "mixed" ) {
                
                QVector pureRef = getRandomStateVector(NUM_QUBITS);
                toQureg(pure, pureRef);
            
                // test when density matrix is mixed 
                QMatrix matRef = getRandomDensityMatrix(NUM_QUBITS);
                toQureg(mat, matRef);
                
                // <pure|mat|pure> = <pure| (Mat|pure>)
                QVector rhs = matRef * pureRef;
                qcomp dotProd = 0;
                for (size_t i=0; i<rhs.size(); i++)
                    dotProd += conj(pureRef[i]) * rhs[i];

                REQUIRE( imag(dotProd) == Approx(0).margin(REAL_EPS) );
                REQUIRE( calcFidelity(mat,pure) == Approx(real(dotProd)).margin(100 * REAL_EPS) );
            }

            // unnormalised test is no longer supported, since v4 calcFidelity
            // validates thet the fidelity is correctly approximately real

                // SECTION( "unnormalised" ) {
                    
                //     // test when both density matrix and pure state are unnormalised
                //     QVector pureRef = getRandomQVector(1<<NUM_QUBITS);
                //     QMatrix matRef = getRandomQMatrix(1<<NUM_QUBITS);
                //     toQureg(pure, pureRef);
                //     toQureg(mat, matRef);
                    
                //     // real[ <pure|mat|pure> ] = real[ <pure| (Mat|pure>) ]
                //     QVector rhs = matRef * pureRef;
                //     qcomp dotProd = 0;
                //     for (size_t i=0; i<rhs.size(); i++)
                //         dotProd += conj(pureRef[i]) * rhs[i];
                    
                //     REQUIRE( calcFidelity(mat,pure) == Approx(real(dotProd)) );
                // }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "dimensions" ) {
            
            // two state-vectors
            Qureg vec2 = createQureg(vec.numQubits + 1);
            REQUIRE_THROWS_WITH( calcFidelity(vec2,vec), ContainsSubstring("differing numbers of qubits") );
            destroyQureg(vec2);
        
            // density-matrix and state-vector
            Qureg mat2 = createDensityQureg(vec.numQubits + 1);
            REQUIRE_THROWS_WITH( calcFidelity(mat2,vec), ContainsSubstring("differing numbers of qubits") );
            destroyQureg(mat2);
        }
        SECTION( "density-matrices" ) {
            
            // two mixed statess
            REQUIRE_THROWS_WITH( calcFidelity(mat,mat), ContainsSubstring("Quregs cannot both be density matrices") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
    destroyQureg(pure);
}



/** @sa calcHilbertSchmidtDistance
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcHilbertSchmidtDistance", "[calculations]" ) {
    
    Qureg mat1 = createForcedDensityQureg(NUM_QUBITS);
    Qureg mat2 = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        // perform these random tests 10 times 
        GENERATE( range(0,10) );
        
        SECTION( "density-matrix" ) {
            
            SECTION( "pure" ) {
                
                // create random |r1><r1| and |r2><r2| states
                QVector r1 = getRandomStateVector(NUM_QUBITS);
                QMatrix m1 = getKetBra(r1,r1);
                toQureg(mat1, m1);
                QVector r2 = getRandomStateVector(NUM_QUBITS);
                QMatrix m2 = getKetBra(r2,r2);
                toQureg(mat2, m2);
                
                // Tr{ (a-b)(a-b)^dagger } =  sum_{ij} |a_{ij} - b_{ij}|^2
                qreal tr = 0;
                for (size_t i=0; i<m1.size(); i++)
                    for (size_t j=0; j<m1.size(); j++)
                        tr += pow(abs(m1[i][j] - m2[i][j]), 2);
                
                qreal res = calcHilbertSchmidtDistance(mat1, mat2);
                REQUIRE( res == Approx(sqrt(tr)) );
            
            }
            SECTION( "normalised" ) {
                
                QMatrix ref1 = getRandomDensityMatrix(NUM_QUBITS);
                QMatrix ref2 = getRandomDensityMatrix(NUM_QUBITS);
                toQureg(mat1, ref1);
                toQureg(mat2, ref2);
                
                // Tr{ (a-b)(a-b)^dagger } =  sum_{ij} |a_{ij} - b_{ij}|^2
                qreal tr = 0;
                for (size_t i=0; i<ref1.size(); i++)
                    for (size_t j=0; j<ref1.size(); j++)
                        tr += pow(abs(ref1[i][j] - ref2[i][j]), 2);
                
                qreal res = calcHilbertSchmidtDistance(mat1, mat2);
                REQUIRE( res == Approx(sqrt(tr)) );
            }
            SECTION( "unnormalised" ) {
                
                // mat1 and mat2 are both random matrices
                QMatrix ref1 = getRandomQMatrix(1<<NUM_QUBITS);
                QMatrix ref2 = getRandomQMatrix(1<<NUM_QUBITS);
                toQureg(mat1, ref1);
                toQureg(mat2, ref2);
                
                // Tr{ (a-b)(a-b)^dagger } =  sum_{ij} |a_{ij} - b_{ij}|^2
                qreal tr = 0;
                for (size_t i=0; i<ref1.size(); i++)
                    for (size_t j=0; j<ref1.size(); j++)
                        tr += pow(abs(ref1[i][j] - ref2[i][j]), 2);
                
                qreal res = calcHilbertSchmidtDistance(mat1, mat2);
                REQUIRE( res == Approx(sqrt(tr)) );
            }
        }
    }
    SECTION( "input validation") {
        
        SECTION( "dimensions" ) {
            
            Qureg mat3 = createDensityQureg(NUM_QUBITS + 1);
            REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(mat1,mat3), ContainsSubstring("differing numbers of qubits") );
            destroyQureg(mat3);
        }

        // in v4, this function calls calcDistance() which has state-vector overloads

            // SECTION( "state-vector" ) {
                
            //     Qureg vec = createForcedQureg(NUM_QUBITS);
                
            //     REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(vec,mat1), ContainsSubstring("valid only for density matrices") );
            //     REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(mat1,vec), ContainsSubstring("valid only for density matrices") );
            //     REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(vec,vec), ContainsSubstring("valid only for density matrices") );
                
            //     destroyQureg(vec);
            // }
    }
    destroyQureg(mat1);
    destroyQureg(mat2);
}



/** @sa calcInnerProduct
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcInnerProduct", "[calculations]" ) {
    
    Qureg vec1 = createForcedQureg(NUM_QUBITS);
    Qureg vec2 = createForcedQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        // perform these random tests 10 times 
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            SECTION( "normalised" ) {
                
                // <r1|r2> = sum_j conj(r1_j) * r2_j
                QVector r1 = getRandomStateVector(NUM_QUBITS);
                QVector r2 = getRandomStateVector(NUM_QUBITS);
                qcomp prod = 0;
                for (size_t i=0; i<r1.size(); i++)
                    prod += conj(r1[i]) * r2[i];
                    
                toQureg(vec1, r1);
                toQureg(vec2, r2);
                qcomp res = calcInnerProduct(vec1,vec2);
                
                REQUIRE( real(res) == Approx(real(prod)) );
                REQUIRE( imag(res) == Approx(imag(prod)) );
            }
            SECTION( "unnormalised" ) {
                
                // <r1|r2> = sum_j conj(r1_j) * r2_j
                QVector r1 = getRandomQVector(1<<NUM_QUBITS);
                QVector r2 = getRandomQVector(1<<NUM_QUBITS);
                qcomp prod = 0;
                for (size_t i=0; i<r1.size(); i++)
                    prod += conj(r1[i]) * r2[i];
                    
                toQureg(vec1, r1);
                toQureg(vec2, r2);
                qcomp res = calcInnerProduct(vec1,vec2);
                
                REQUIRE( real(res) == Approx(real(prod)) );
                REQUIRE( imag(res) == Approx(imag(prod)) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "dimensions" ) {
            
            Qureg vec3 = createQureg(NUM_QUBITS + 1);
            REQUIRE_THROWS_WITH( calcInnerProduct(vec1,vec3), ContainsSubstring("differing numbers of qubits") );
            destroyQureg(vec3);
        }

        // density-matrix arguments are permitted in v4

            // SECTION( "density-matrix" ) {
                
            //     Qureg mat = createForcedDensityQureg(NUM_QUBITS);
                
            //     REQUIRE_THROWS_WITH( calcInnerProduct(vec1,mat), ContainsSubstring("valid only for state-vectors") );
            //     REQUIRE_THROWS_WITH( calcInnerProduct(mat,vec1), ContainsSubstring("valid only for state-vectors") );
            //     REQUIRE_THROWS_WITH( calcInnerProduct(mat,mat), ContainsSubstring("valid only for state-vectors") );
                
            //     destroyQureg(mat);
            // }
    }
    destroyQureg(vec1);
    destroyQureg(vec2);
}



/** @sa calcProbOfAllOutcomes
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcProbOfAllOutcomes", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
    
        // generate all possible qubit arrangements
        int numQubits = GENERATE_COPY( range(1,NUM_QUBITS+1) );
        int* qubits = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numQubits) );
        
        int numOutcomes = 1<<numQubits;
        vector<qreal> probs(numOutcomes);
        QVector refProbs = QVector(numOutcomes);
            
        SECTION( "state-vector" ) {
            
            SECTION( "normalised" ) {
                
                QVector ref = getRandomStateVector(NUM_QUBITS);
                toQureg(vec, ref);

                // prob is sum of |amp|^2 of basis states which encode outcome
                for (size_t i=0; i<ref.size(); i++) {
                    int outcome = 0;
                    for (int q=0; q<numQubits; q++) {
                        int bit = (i >> qubits[q]) & 1;
                        outcome += bit * (1 << q);
                    }
                    refProbs[outcome] += pow(abs(ref[i]), 2);
                }

                calcProbOfAllOutcomes(probs.data(), vec, qubits, numQubits);
                REQUIRE( areEqual(refProbs, probs.data()) );
            }
            SECTION( "unnormalised" ) {
                
                QVector ref = getRandomQVector(1<<NUM_QUBITS);
                toQureg(vec, ref);
                
                // prob is sum of |amp|^2 of basis states which encode outcome
                for (size_t i=0; i<ref.size(); i++) {
                    int outcome = 0;
                    for (int q=0; q<numQubits; q++) {
                        int bit = (i >> qubits[q]) & 1;
                        outcome += bit * (1 << q);
                    }
                    refProbs[outcome] += pow(abs(ref[i]), 2);
                }

                calcProbOfAllOutcomes(probs.data(), vec, qubits, numQubits);
                REQUIRE( areEqual(refProbs, probs.data()) );
            }
        }
        SECTION( "density-matrix" ) {
            
            SECTION( "normalised" ) {
            
                QMatrix ref = getRandomDensityMatrix(NUM_QUBITS);
                toQureg(mat, ref);
                
                // prob is sum of diagonals which encode outcome 
                for (size_t i=0; i<ref.size(); i++) {
                    int outcome = 0;
                    for (int q=0; q<numQubits; q++) {
                        int bit = (i >> qubits[q]) & 1;
                        outcome += bit * (1 << q);
                    }
                    refProbs[outcome] += real(ref[i][i]);
                }
                
                calcProbOfAllOutcomes(probs.data(), mat, qubits, numQubits);
                REQUIRE( areEqual(refProbs, probs.data()) );
            }
            SECTION( "unnormalised" ) {
            
                QMatrix ref = getRandomQMatrix(1<<NUM_QUBITS);
                toQureg(mat, ref);
                
                // prob is sum of diagonals which encode outcome 
                for (size_t i=0; i<ref.size(); i++) {
                    int outcome = 0;
                    for (int q=0; q<numQubits; q++) {
                        int bit = (i >> qubits[q]) & 1;
                        outcome += bit * (1 << q);
                    }
                    refProbs[outcome] += real(ref[i][i]);
                }
                
                calcProbOfAllOutcomes(probs.data(), mat, qubits, numQubits);
                REQUIRE( areEqual(refProbs, probs.data()) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        int numQubits = 3;
        int qubits[] = {0, 1, 2};
        qreal probs[8];
        
        SECTION( "number of qubits" ) {
            
            // too small
            numQubits = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( calcProbOfAllOutcomes(probs, mat, qubits, numQubits), ContainsSubstring("specified number of target qubits") && ContainsSubstring("invalid.") );

            // too big
            numQubits = NUM_QUBITS + 1;
            REQUIRE_THROWS_WITH( calcProbOfAllOutcomes(probs, mat, qubits, numQubits), ContainsSubstring("target qubits") && ContainsSubstring("exceeds the number of qubits in the Qureg") );
        }
        SECTION( "qubit indices" ) {
            
            qubits[GENERATE_COPY(range(0,numQubits))] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( calcProbOfAllOutcomes(probs, mat, qubits, numQubits), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "repetition of qubits" ) {
            
            qubits[GENERATE_COPY(1,2)] = qubits[0];
            REQUIRE_THROWS_WITH( calcProbOfAllOutcomes(probs, mat, qubits, numQubits), ContainsSubstring("qubits must be unique") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat); 
}
        



/** @sa calcProbOfOutcome
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcProbOfOutcome", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        int target = GENERATE( range(0,NUM_QUBITS) );
        int outcome = GENERATE( 0, 1 );
        
        SECTION( "state-vector" ) {
            
            SECTION( "normalised" ) {
                
                QVector ref = getRandomStateVector(NUM_QUBITS);
                toQureg(vec, ref);
                
                // prob is sum of |amp|^2 of amplitudes where target bit is outcome
                qreal prob = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == outcome)
                        prob += pow(abs(ref[ind]), 2);
                }
                
                REQUIRE( calcProbOfOutcome(vec, target, outcome) == Approx(prob) );
            }
            SECTION( "unnormalised" ) {
                
                QVector ref = getRandomQVector(1<<NUM_QUBITS);
                toQureg(vec, ref);
                
                // prob is sum of |amp|^2 of amplitudes where target bit is outcome
                qreal prob = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == outcome)
                        prob += pow(abs(ref[ind]), 2);
                }
                
                REQUIRE( calcProbOfOutcome(vec, target, outcome) == Approx(prob) );
            }
        }
        SECTION( "density-matrix" ) {

            SECTION( "pure" ) {
                
                // set mat to a random |r><r|
                QVector ref = getRandomStateVector(NUM_QUBITS);
                toQureg(mat, getKetBra(ref, ref));
                
                // calc prob of the state-vector
                qreal prob = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == outcome)
                        prob += pow(abs(ref[ind]), 2);
                }
                
                REQUIRE( calcProbOfOutcome(mat, target, outcome) == Approx(prob) );
            }
            SECTION( "mixed" ) {
    
                QMatrix ref = getRandomDensityMatrix(NUM_QUBITS);
                toQureg(mat, ref);
                
                // prob is sum of diagonal amps (should be real) where target bit is outcome
                qcomp tr = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == outcome)
                        tr += ref[ind][ind];
                }
                
                REQUIRE( imag(tr) == Approx(0).margin(REAL_EPS) );
                
                REQUIRE( calcProbOfOutcome(mat, target, outcome) == Approx(real(tr)) );
            }
            SECTION( "unnormalised" ) {
                
                QMatrix ref = getRandomQMatrix(1<<NUM_QUBITS);
                toQureg(mat, ref);
                
                // prob is (sum of real of diagonal amps where target bit is outcome)
                qreal tr = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == outcome)
                        tr += real(ref[ind][ind]);
                }
                
                REQUIRE( calcProbOfOutcome(mat, target, outcome) == Approx(tr) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( calcProbOfOutcome(vec, target, 0), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "outcome value" ) {
            
            int outcome = GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( calcProbOfOutcome(vec, 0, outcome), ContainsSubstring("measurement outcome") && ContainsSubstring("invalid") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa calcPurity
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcPurity", "[calculations]" ) {
    
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        // perform the following random tests 10 times 
        GENERATE( range(1,10) );
        
        SECTION( "density-matrix" ) {
            
            SECTION( "pure" ) {
            
                // pure states have unity purity 
                initZeroState(mat);
                REQUIRE( calcPurity(mat) == 1 );
                
                // (try also a pure random L2-vector)
                QVector r1 = getRandomStateVector(NUM_QUBITS); // |r>
                QMatrix m1 = getKetBra(r1, r1); // |r><r|
                toQureg(mat, m1);
                REQUIRE( calcPurity(mat) == Approx(1) );
            
            }
            SECTION( "mixed" ) {
            
                // mixed states have 1/2^N < purity < 1
                QMatrix ref = getRandomDensityMatrix(NUM_QUBITS);
                toQureg(mat, ref);
                qreal purity = calcPurity(mat);
                REQUIRE( purity < 1 );
                REQUIRE( purity >= 1/pow(2.,NUM_QUBITS) );
                
                // compare to Tr(rho^2)
                QMatrix prod = ref*ref;
                qreal tr = 0;
                for (size_t i=0; i<prod.size(); i++)
                    tr += real(prod[i][i]);
                REQUIRE( purity == Approx(tr) );
            }
            SECTION( "unnormalised" ) {
            
                // unphysical states give sum_{ij} |rho_ij|^2
                QMatrix ref = getRandomQMatrix(1<<NUM_QUBITS);
                qreal tot = 0;
                for (size_t i=0; i<ref.size(); i++)
                    for (size_t j=0; j<ref.size(); j++)
                        tot += pow(abs(ref[i][j]), 2);
                    
                toQureg(mat, ref);
                REQUIRE( calcPurity(mat) == Approx(tot) );
            }
        }
    }
    SECTION( "input validation" ) {

        // in v4, this accepts state-vectors
        
            // SECTION( "state-vector" ) {
                
            //     Qureg vec = createForcedQureg(NUM_QUBITS);
            //     REQUIRE_THROWS_WITH( calcPurity(vec), ContainsSubstring("valid only for density matrices") );
            //     destroyQureg(vec);
            // }
    }
    destroyQureg(mat);
}



/** @sa calcTotalProb
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "calcTotalProb", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
        
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            // normalised: prob(vec) = 1
            initPlusState(vec);
            REQUIRE( calcTotalProb(vec) == Approx(1) );
            
            // zero norm: prob(vec) = 0
            initBlankState(vec);
            REQUIRE( calcTotalProb(vec) == 0 );
            
            // random L2 state: prob(vec) = 1
            toQureg(vec, getRandomStateVector(NUM_QUBITS));
            REQUIRE( calcTotalProb(vec) == Approx(1) );
            
            // unnormalised: prob(vec) = sum_i |vec_i|^2
            initDebugState(vec);
            QVector ref = toQVector(vec);
            qreal refProb = 0;
            for (size_t i=0; i<ref.size(); i++)
                refProb += pow(abs(ref[i]), 2);
            REQUIRE( calcTotalProb(vec) == Approx(refProb) );
        }
        SECTION( "density-matrix" ) {
            
            // normalised: prob(mat) = 1
            initPlusState(mat);
            REQUIRE( calcTotalProb(mat) == Approx(1) );
            
            // zero norm: prob(mat) = 0
            initBlankState(mat);
            REQUIRE( calcTotalProb(mat) == 0 );
            
            // random density matrix: prob(mat) = 1
            toQureg(mat, getRandomDensityMatrix(NUM_QUBITS));
            REQUIRE( calcTotalProb(mat) == Approx(1) );
            
            // unnormalised: prob(mat) = sum_i real(mat_{ii})
            initDebugState(mat);
            QMatrix ref = toQMatrix(mat);
            qreal refProb = 0;
            for (size_t i=0; i<ref.size(); i++)
                refProb += real(ref[i][i]);
            REQUIRE( calcTotalProb(mat) == Approx(refProb) );
        }
    }
    SECTION( "input validation" ) {
        
        // no validation 
        SUCCEED();
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa getAmp
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getAmp", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initDebugState(vec);
            QVector ref = toQVector(vec);

            int ind = GENERATE( range(0,1<<NUM_QUBITS) );
            qcomp amp = getAmp(vec,ind);
            REQUIRE( amp == ref[ind] );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getAmp(vec,ind), ContainsSubstring("Basis state index") && ContainsSubstring("invalid") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createForcedDensityQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( getAmp(mat,0), ContainsSubstring("Expected a statevector Qureg but received a density matrix") );
            destroyQureg(mat);
        }
    }
    destroyQureg(vec);
}



/** @sa getDensityAmp
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getDensityAmp", "[calculations]" ) {
    
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "density-matrix" ) {
            
            initDebugState(mat);
            QMatrix ref = toQMatrix(mat);

            int row = GENERATE( range(0,1<<NUM_QUBITS) );
            int col = GENERATE( range(0,1<<NUM_QUBITS) );
            
            qcomp amp = getDensityAmp(mat,row,col);
            REQUIRE( amp == ref[row][col] );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getDensityAmp(mat,ind,0), ContainsSubstring("The row and column indices") && ContainsSubstring("invalid") );
            REQUIRE_THROWS_WITH( getDensityAmp(mat,0,ind), ContainsSubstring("The row and column indices") && ContainsSubstring("invalid") );

        }
        SECTION( "state-vector" ) {
            
            Qureg vec = createForcedQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( getDensityAmp(vec,0,0), ContainsSubstring("Expected a density matrix Qureg but received a statevector") );
            destroyQureg(vec);
        }
    }
    destroyQureg(mat);
}



/** @sa getImagAmp
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getImagAmp", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initDebugState(vec);
            QVector ref = toQVector(vec);

            int ind = GENERATE( range(0,1<<NUM_QUBITS) );
            REQUIRE( getImagAmp(vec,ind) == imag(ref[ind]) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getImagAmp(vec,ind), ContainsSubstring("Basis state index") && ContainsSubstring("invalid") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createForcedDensityQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( getImagAmp(mat,0), ContainsSubstring("Expected a statevector Qureg but received a density matrix") );
            destroyQureg(mat);
        }
    }
    destroyQureg(vec);
}



/** @sa getNumAmps
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getNumAmps", "[calculations]" ) {
        
    SECTION( "correctness" ) {
        
        // test >= NUM_QUBITS so as not to limit distribution size
        int numQb = GENERATE( range(NUM_QUBITS, NUM_QUBITS+10) );
        
        SECTION( "state-vector" ) {
            
            Qureg vec = createForcedQureg(numQb);
            REQUIRE( getNumAmps(vec) == (1<<numQb) );
            destroyQureg(vec);
        }
    }

    // in v4, argument can be a density matrix (return total num amps)

        // SECTION( "input validation" ) {
            
        //     SECTION( "density-matrix" ) {
        //         Qureg mat = createForcedDensityQureg(NUM_QUBITS);
        //         REQUIRE_THROWS_WITH( getNumAmps(mat), ContainsSubstring("Expected a statevector Qureg but received a density matrix") );
        //         destroyQureg(mat);
        //     }
        // }
}



/** @sa getNumQubits
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getNumQubits", "[calculations]" ) {
        
    SECTION( "correctness" ) {
        
        // test >= NUM_QUBITS so as not to limit distribution size
        int numQb = GENERATE( range(NUM_QUBITS, NUM_QUBITS+10) );
        
        SECTION( "state-vector" ) {
            
            Qureg vec = createForcedQureg(numQb);
            REQUIRE( getNumQubits(vec) == numQb );
            destroyQureg(vec);
        }
        SECTION( "density-matrix" ) {
            
            // density matrices use square as much memory; we must be careful not to seg-fault!
            if (2*numQb > 25)
                numQb = 13; // max size

            Qureg mat = createForcedDensityQureg(numQb);
            REQUIRE( getNumQubits(mat) == numQb );
            destroyQureg(mat);
        }
    }
    SECTION( "input validation" ) {
        
        // no validation
        SUCCEED();
    }
}



/** @sa getProbAmp
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getProbAmp", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initDebugState(vec);
            QVector ref = toQVector(vec);

            int ind = GENERATE( range(0,1<<NUM_QUBITS) );
            qreal refCalc = pow(abs(ref[ind]), 2);
            REQUIRE( getProbAmp(vec,ind) == Approx(refCalc) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getProbAmp(vec,ind), ContainsSubstring("Basis state index") && ContainsSubstring("invalid") );
        }

        // in v4, this redirects to calcProbOfBasisState() which accepts
        // both statevectors and density matrices

            // SECTION( "density-matrix" ) {
                
            //     Qureg mat = createForcedDensityQureg(NUM_QUBITS);
            //     REQUIRE_THROWS_WITH( getProbAmp(mat,0), ContainsSubstring("Expected a statevector Qureg but received a density matrix") );
            //     destroyQureg(mat);
            // }
    }
    destroyQureg(vec);
}



/** @sa getRealAmp
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "getRealAmp", "[calculations]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initDebugState(vec);
            QVector ref = toQVector(vec);

            int ind = GENERATE( range(0,1<<NUM_QUBITS) );
            REQUIRE( getRealAmp(vec,ind) == real(ref[ind]) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getRealAmp(vec,ind), ContainsSubstring("Basis state index") && ContainsSubstring("invalid") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createForcedDensityQureg(NUM_QUBITS);
            REQUIRE_THROWS_WITH( getRealAmp(mat,0), ContainsSubstring("Expected a statevector Qureg but received a density matrix") );
            destroyQureg(mat);
        }
    }
    destroyQureg(vec);
}


