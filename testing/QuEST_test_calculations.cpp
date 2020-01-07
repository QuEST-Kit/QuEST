
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

/** Prepares the needed data structures for unit testing. This creates 
 * the QuEST environment, and a state-vector and density matrix of the size 'numQb',
 * both of which are (automatically) initialised to the zero state.
 * numQb should be NUM_QUBITS unless motivated otherwise.
 */
#define PREPARE_TEST(env, vec, mat, numQb) \
    QuESTEnv env = createQuESTEnv(); \
    Qureg vec = createQureg(numQb, env); \
    Qureg mat = createDensityQureg(numQb, env);

/** Destroys the data structures made by PREPARE_TEST */
#define CLEANUP_TEST(env, vec, mat) \
    destroyQureg(vec, env); \
    destroyQureg(mat, env); \
    destroyQuESTEnv(env);

/* allows concise use of Contains in catch's REQUIRE_THROWS_WITH */
using Catch::Matchers::Contains;


TEST_CASE( "calcDensityInnerProduct", "[calculations]" ) {

    QuESTEnv env = createQuESTEnv();
    Qureg mat1 = createDensityQureg(NUM_QUBITS, env);
    Qureg mat2 = createDensityQureg(NUM_QUBITS, env);
    
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
                
                REQUIRE( calcDensityInnerProduct(mat1,mat2) == Approx(densProd) );
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
                
                REQUIRE( calcDensityInnerProduct(mat1,mat2) == Approx(real(refProd)) );
                
                // should be invariant under ordering
                REQUIRE( calcDensityInnerProduct(mat1,mat2) == Approx(calcDensityInnerProduct(mat2,mat1)) );
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
                        
                REQUIRE( calcDensityInnerProduct(mat1,mat2) == Approx(real(refProd)) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "dimensions" ) {
            
            Qureg mat3 = createDensityQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcDensityInnerProduct(mat1,mat3), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(mat3, env);
        }
        SECTION( "state-vectors" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            
            REQUIRE_THROWS_WITH( calcDensityInnerProduct(mat1,vec), Contains("valid only for density matrices") );
            REQUIRE_THROWS_WITH( calcDensityInnerProduct(vec,mat1), Contains("valid only for density matrices") );
            REQUIRE_THROWS_WITH( calcDensityInnerProduct(vec,vec),  Contains("valid only for density matrices") );        
            
            destroyQureg(vec, env);
        }
    }
    destroyQureg(mat1, env);
    destroyQureg(mat2, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcExpecPauliProd", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    initDebugState(vec);
    initDebugState(mat);
    QVector vecRef = toQVector(vec);
    QMatrix matRef = toQMatrix(mat);
    
    Qureg vecWork = createQureg(NUM_QUBITS, env);
    Qureg matWork = createDensityQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        int numTargs = GENERATE( range(1,NUM_QUBITS+1) );
        int* targs = GENERATE_COPY( sublists(range(0,NUM_QUBITS), numTargs) );
        
        /* it's too expensive to try ALL Pauli sequences, via 
         *      pauliOpType* paulis = GENERATE_COPY( pauliseqs(numTargs) );.
         * Furthermore, take(10, pauliseqs(numTargs)) will try the same pauli codes.
         * Hence, we instead opt to repeatedlyrandomly generate pauliseqs
         */
        GENERATE( range(0,10) ); // gen 10 random pauli-codes for every targs
        pauliOpType paulis[numTargs];
        for (int i=0; i<numTargs; i++)
            paulis[i] = (pauliOpType) getRandomInt(0,4);
            
        // produce a numTargs-big matrix 'pauliProd' by pauli-matrix tensoring
        QMatrix iMatr{{1,0},{0,1}};
        QMatrix xMatr{{0,1},{1,0}};
        QMatrix yMatr{{0,-1i},{1i,0}};
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
            
            qreal res = calcExpecPauliProd(vec, targs, paulis, numTargs, vecWork);
            REQUIRE( res == Approx(real(prod)).margin(REAL_EPS) );
        }
        SECTION( "density matrix" ) {
            
            /* calcExpecPauliProd calculates Trace( pauliProd * qureg ) */

            // produce (pauliProd * mat)
            QMatrix fullOp = getFullOperatorMatrix(NULL, 0, targs, numTargs, pauliProd, NUM_QUBITS);
            matRef = fullOp * matRef;
            
            // compute real(trace(pauliProd * mat))
            qreal tr = 0;
            for (size_t i=0; i<matRef.size(); i++)
                tr += real(matRef[i][i]);
            // (get real, since we start in a non-Hermitian state, hence diagonal isn't real)
            
            qreal res = calcExpecPauliProd(mat, targs, paulis, numTargs, matWork);
            REQUIRE( res == Approx(tr).margin(10*REAL_EPS) );
        }
    }
    SECTION( "validation" ) {
        
        SECTION( "number of targets" ) {
            
            int numTargs = GENERATE( -1, 0, NUM_QUBITS+1 );
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, NULL, NULL, numTargs, vecWork), Contains("Invalid number of target") );
        }
        SECTION( "target indices" ) {
            
            int numTargs = 3;
            int targs[3] = {0, 1, 2};
            
            // make one index wrong
            targs[GENERATE( range(0,3) )] = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, NULL, numTargs, vecWork), Contains("Invalid target qubit") );
        }
        SECTION( "repetition in targets" ) {
            
            int numTargs = 3;
            int targs[3] = {0, 1, 1};
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, NULL, numTargs, vecWork), Contains("target qubits must be unique") );
        }
        SECTION( "pauli codes" ) {
            
            int numTargs = 3;
            int targs[3] = {0, 1, 2};
            pauliOpType codes[3] = {PAULI_X, PAULI_Y, PAULI_Z};
            
            // make one pauli wrong
            codes[GENERATE( range(0,3) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, codes, numTargs, vecWork), Contains("Invalid Pauli code") );
        }
        SECTION( "workspace type" ) {
            
            int numTargs = 1;
            int targs[1] = {0};
            pauliOpType codes[1] = {PAULI_I};
            
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, codes, numTargs, matWork), Contains("Registers must both be state-vectors or both be density matrices") );
            REQUIRE_THROWS_WITH( calcExpecPauliProd(mat, targs, codes, numTargs, vecWork), Contains("Registers must both be state-vectors or both be density matrices") );
        }
        SECTION( "workspace dimensions" ) {
                
            int numTargs = 1;
            int targs[1] = {0};
            pauliOpType codes[1] = {PAULI_I};
    
            Qureg vec2 = createQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcExpecPauliProd(vec, targs, codes, numTargs, vec2), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(vec2, env);
            
            Qureg mat2 = createDensityQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcExpecPauliProd(mat, targs, codes, numTargs, mat2), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(mat2, env);
        }
    }
    destroyQureg(vec, env);
    destroyQureg(mat, env);
    destroyQureg(vecWork, env);
    destroyQureg(matWork, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcExpecPauliSum", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    initDebugState(vec);
    initDebugState(mat);
    QVector vecRef = toQVector(vec);
    QMatrix matRef = toQMatrix(mat);
    
    Qureg vecWork = createQureg(NUM_QUBITS, env);
    Qureg matWork = createDensityQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        // try 1 to 5 sums of pauli products
        int numSumTerms = GENERATE( range(1,5+1) );
        
        /* it's too expensive to try every possible Pauli configuration, so
         * we'll try 10 random codes, and for each, random coefficients
         */
        GENERATE( range(0,10) );
        int totNumCodes = numSumTerms*NUM_QUBITS;
        pauliOpType paulis[totNumCodes];
        for (int i=0; i<totNumCodes; i++)
            paulis[i] = (pauliOpType) getRandomInt(0,4);
            
        // for every above param configuration, try random coefficients
        qreal coeffs[numSumTerms];
        for (int i=0; i<numSumTerms; i++)
            coeffs[i] = getRandomReal(-5, 5);
            
        // produce a numTargs-big matrix 'pauliSum' by pauli-matrix tensoring and summing
        QMatrix iMatr{{1,0},{0,1}};
        QMatrix xMatr{{0,1},{1,0}};
        QMatrix yMatr{{0,-1i},{1i,0}};
        QMatrix zMatr{{1,0},{0,-1}};
        QMatrix pauliSum = getZeroMatrix(1<<NUM_QUBITS);
        
        for (int t=0; t<numSumTerms; t++) {
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

            /* calcExpecPauliSum calculates <qureg|pauliSum|qureg> */
            
            QVector sumRef = pauliSum * vecRef;
            qcomp prod = 0;
            for (size_t i=0; i<vecRef.size(); i++)
                prod += conj(vecRef[i]) * sumRef[i];
            REQUIRE( imag(prod) == Approx(0).margin(10*REAL_EPS) );
            
            qreal res = calcExpecPauliSum(vec, paulis, coeffs, numSumTerms, vecWork);
            REQUIRE( res == Approx(real(prod)).margin(10*REAL_EPS) );
        } 
        SECTION( "density matrix" ) {
            
            /* calcExpecPauliSum calculates Trace( pauliSum * qureg ) */
            matRef = pauliSum * matRef;            
            qreal tr = 0;
            for (size_t i=0; i<matRef.size(); i++)
                tr += real(matRef[i][i]);
            // (get real, since we start in a non-Hermitian state, hence diagonal isn't real)
            
            qreal res = calcExpecPauliSum(mat, paulis, coeffs, numSumTerms, matWork);
            REQUIRE( res == Approx(tr).margin(1E2*REAL_EPS) );
        }
    }
    SECTION( "validation" ) {
        
        SECTION( "number of sum terms" ) {
            
            int numSumTerms = GENERATE( -1, 0 );
            REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, NULL, NULL, numSumTerms, vecWork), Contains("Invalid number of terms in the Pauli sum") );
        }
        SECTION( "pauli codes" ) {
            
            // make valid params
            int numSumTerms = 3;
            qreal coeffs[numSumTerms];
            pauliOpType codes[numSumTerms*NUM_QUBITS];
            for (int i=0; i<numSumTerms*NUM_QUBITS; i++)
                codes[i] = PAULI_I;

            // make one pauli wrong
            codes[GENERATE_COPY( range(0,numSumTerms*NUM_QUBITS) )] = (pauliOpType) GENERATE( -1, 4 );
            REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, codes, coeffs, numSumTerms, vecWork), Contains("Invalid Pauli code") );
        }
        SECTION( "workspace type" ) {
            
            // make valid params
            int numSumTerms = 1;
            qreal coeffs[1] = {0};
            pauliOpType codes[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                codes[i] = PAULI_I;
            
            REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, codes, coeffs, numSumTerms, mat), Contains("Registers must both be state-vectors or both be density matrices") );
            REQUIRE_THROWS_WITH( calcExpecPauliSum(mat, codes, coeffs, numSumTerms, vec), Contains("Registers must both be state-vectors or both be density matrices") );
        }
        SECTION( "workspace dimensions" ) {
                
            // make valid params
            int numSumTerms = 1;
            qreal coeffs[1] = {0};
            pauliOpType codes[NUM_QUBITS];
            for (int i=0; i<NUM_QUBITS; i++)
                codes[i] = PAULI_I;
    
            Qureg vec2 = createQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcExpecPauliSum(vec, codes, coeffs, numSumTerms, vec2), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(vec2, env);
            
            Qureg mat2 = createDensityQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcExpecPauliSum(mat, codes, coeffs, numSumTerms, mat2), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(mat2, env);
        }
    }
    destroyQureg(vec, env);
    destroyQureg(mat, env);
    destroyQureg(vecWork, env);
    destroyQureg(matWork, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcFidelity", "[calculations]" ) {
    
    PREPARE_TEST(env, vec, mat, NUM_QUBITS);
    Qureg pure = createQureg(NUM_QUBITS, env);
    
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
            SECTION( "unnormalised" ) {
            
                // random unnormalised vectors
                QVector vecRef = getRandomQVector(1<<NUM_QUBITS);
                QVector pureRef = getRandomQVector(1<<NUM_QUBITS);
                toQureg(vec, vecRef);
                toQureg(pure, pureRef);
                
                // Let nv be magnitude of vec, hence |unit-vec> = 1/sqrt(nv)|vec>
                qreal nv = 0;
                for (size_t i=0; i<vecRef.size(); i++)
                    nv += pow(abs(vecRef[i]), 2);
                // then <vec|vec> = sqrt(nv)*sqrt(nv) <unit-vec|unit-vec> = nv,
                // hence |<vec|vec>|^2 = nv*nv
                REQUIRE( calcFidelity(vec,vec) == Approx( nv*nv ) );
                
                qcomp dotProd = 0;
                for (size_t i=0; i<vecRef.size(); i++)
                    dotProd += conj(vecRef[i]) * pureRef[i];
                qreal refFid = pow(abs(dotProd), 2);
                
                REQUIRE( calcFidelity(vec,pure) == Approx(refFid) ); 
            }
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
                
                REQUIRE( calcFidelity(mat,pure) == Approx(refFid) ); 
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
                REQUIRE( calcFidelity(mat,pure) == Approx(real(dotProd)) );
            }
            SECTION( "unnormalised" ) {
                
                // test when both density matrix and pure state are unnormalised
                QVector pureRef = getRandomQVector(1<<NUM_QUBITS);
                QMatrix matRef = getRandomQMatrix(1<<NUM_QUBITS);
                toQureg(pure, pureRef);
                toQureg(mat, matRef);
                
                // real[ <pure|mat|pure> ] = real[ <pure| (Mat|pure>) ]
                QVector rhs = matRef * pureRef;
                qcomp dotProd = 0;
                for (size_t i=0; i<rhs.size(); i++)
                    dotProd += conj(pureRef[i]) * rhs[i];
                
                REQUIRE( calcFidelity(mat,pure) == Approx(real(dotProd)) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "dimensions" ) {
            
            // two state-vectors
            Qureg vec2 = createQureg(vec.numQubitsRepresented + 1, env);
            REQUIRE_THROWS_WITH( calcFidelity(vec2,vec), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(vec2, env);
        
            // density-matrix and state-vector
            Qureg mat2 = createDensityQureg(vec.numQubitsRepresented + 1, env);
            REQUIRE_THROWS_WITH( calcFidelity(mat2,vec), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(mat2, env);
        }
        SECTION( "density-matrices" ) {
            
            // two mixed statess
            REQUIRE_THROWS_WITH( calcFidelity(mat,mat), Contains("Second argument must be a state-vector") );
        }
    }
    destroyQureg(pure, env);
    CLEANUP_TEST(env, vec, mat);
}



TEST_CASE( "calcHilbertSchmidtDistance", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg mat1 = createDensityQureg(NUM_QUBITS, env);
    Qureg mat2 = createDensityQureg(NUM_QUBITS, env);
    
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
            
            Qureg mat3 = createDensityQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(mat1,mat3), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(mat3, env);
        }
        SECTION( "state-vector" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            
            REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(vec,mat1), Contains("valid only for density matrices") );
            REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(mat1,vec), Contains("valid only for density matrices") );
            REQUIRE_THROWS_WITH( calcHilbertSchmidtDistance(vec,vec), Contains("valid only for density matrices") );
            
            destroyQureg(vec, env);
        }
    }
    destroyQureg(mat1, env);
    destroyQureg(mat2, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcInnerProduct", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec1 = createQureg(NUM_QUBITS, env);
    Qureg vec2 = createQureg(NUM_QUBITS, env);
    
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
                Complex res = calcInnerProduct(vec1,vec2);
                
                REQUIRE( res.real == Approx(real(prod)) );
                REQUIRE( res.imag == Approx(imag(prod)) );
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
                Complex res = calcInnerProduct(vec1,vec2);
                
                REQUIRE( res.real == Approx(real(prod)) );
                REQUIRE( res.imag == Approx(imag(prod)) );
            }
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "dimensions" ) {
            
            Qureg vec3 = createQureg(NUM_QUBITS + 1, env);
            REQUIRE_THROWS_WITH( calcInnerProduct(vec1,vec3), Contains("Dimensions") && Contains("don't match") );
            destroyQureg(vec3, env);
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(NUM_QUBITS, env);
            
            REQUIRE_THROWS_WITH( calcInnerProduct(vec1,mat), Contains("valid only for state-vectors") );
            REQUIRE_THROWS_WITH( calcInnerProduct(mat,vec1), Contains("valid only for state-vectors") );
            REQUIRE_THROWS_WITH( calcInnerProduct(mat,mat), Contains("valid only for state-vectors") );
            
            destroyQureg(mat, env);
        }
    }
    destroyQureg(vec1, env);
    destroyQureg(vec2, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcProbOfOutcome", "[calculations]" ) {
    
    PREPARE_TEST(env, vec, mat, NUM_QUBITS);

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
                
                // prob is (sum of |amp|^2 of amplitudes where target bit is zero)
                // or 1 - (this) if outcome == 1
                qreal prob = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == 0)
                        prob += pow(abs(ref[ind]), 2);
                }
                if (outcome == 1)
                    prob = 1 - prob;
                
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
                // or 1 - (this) if outcome == 1
                qreal tr = 0;
                for (size_t ind=0; ind<ref.size(); ind++) {
                    int bit = (ind >> target) & 1; // target-th bit
                    if (bit == 0)
                        tr += real(ref[ind][ind]);
                }
                if (outcome == 1)
                    tr = 1 - tr;
                
                REQUIRE( calcProbOfOutcome(mat, target, outcome) == Approx(tr) );
            }
        }
    }
    SECTION( "validation" ) {
        
        SECTION( "qubit indices" ) {
            
            int target = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( calcProbOfOutcome(vec, target, 0), Contains("Invalid target qubit") );
        }
        SECTION( "outcome value" ) {
            
            int outcome = GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( calcProbOfOutcome(vec, 0, outcome), Contains("Invalid measurement outcome") );
        }
    }
    CLEANUP_TEST(env, vec, mat);
}



TEST_CASE( "calcPurity", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    
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
        
        SECTION( "state-vector" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( calcPurity(vec), Contains("valid only for density matrices") );
            destroyQureg(vec, env);
        }
    }
    destroyQureg(mat, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcTotalProb", "[calculations]" ) {
    
    PREPARE_TEST(env, vec, mat, NUM_QUBITS);
    
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
    CLEANUP_TEST(env, vec, mat);
}



TEST_CASE( "getAmp", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            initDebugState(vec);
            QVector ref = toQVector(vec);

            int ind = GENERATE( range(0,1<<NUM_QUBITS) );
            Complex amp = getAmp(vec,ind);
            REQUIRE( fromComplex(amp) == ref[ind] );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getAmp(vec,ind), Contains("Invalid amplitude index") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( getAmp(mat,0), Contains("valid only for state-vectors") );
            destroyQureg(mat, env);
        }
    }
    destroyQureg(vec, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "getDensityAmp", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        SECTION( "density-matrix" ) {
            
            initDebugState(mat);
            QMatrix ref = toQMatrix(mat);

            int row = GENERATE( range(0,1<<NUM_QUBITS) );
            int col = GENERATE( range(0,1<<NUM_QUBITS) );
            
            Complex amp = getDensityAmp(mat,row,col);
            REQUIRE( fromComplex(amp) == ref[row][col] );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "state index" ) {
            
            int ind = GENERATE( -1, 1<<NUM_QUBITS );
            REQUIRE_THROWS_WITH( getDensityAmp(mat,ind,0), Contains("Invalid amplitude index") );
            REQUIRE_THROWS_WITH( getDensityAmp(mat,0,ind), Contains("Invalid amplitude index") );

        }
        SECTION( "state-vector" ) {
            
            Qureg vec = createQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( getDensityAmp(vec,0,0), Contains("valid only for density matrices") );
            destroyQureg(vec, env);
        }
    }
    destroyQureg(mat, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "getImagAmp", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    
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
            REQUIRE_THROWS_WITH( getImagAmp(vec,ind), Contains("Invalid amplitude index") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( getImagAmp(mat,0), Contains("valid only for state-vectors") );
            destroyQureg(mat, env);
        }
    }
    destroyQureg(vec, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "getNumAmps", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    
    SECTION( "correctness" ) {
        
        // test >= NUM_QUBITS so as not to limit distribution size
        int numQb = GENERATE( range(NUM_QUBITS, NUM_QUBITS+10) );
        
        SECTION( "state-vector" ) {
            
            Qureg vec = createQureg(numQb, env);
            REQUIRE( getNumAmps(vec) == (1<<numQb) );
            destroyQureg(vec, env);
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "density-matrix" ) {
            Qureg mat = createDensityQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( getNumAmps(mat), Contains("valid only for state-vectors") );
            destroyQureg(mat, env);
        }
    }
    destroyQuESTEnv(env);
}



TEST_CASE( "getNumQubits", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    
    SECTION( "correctness" ) {
        
        // test >= NUM_QUBITS so as not to limit distribution size
        int numQb = GENERATE( range(NUM_QUBITS, NUM_QUBITS+10) );
        
        SECTION( "state-vector" ) {
            
            Qureg vec = createQureg(numQb, env);
            REQUIRE( getNumQubits(vec) == numQb );
            destroyQureg(vec, env);
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(numQb, env);
            REQUIRE( getNumQubits(mat) == numQb );
            destroyQureg(mat, env);
        }
    }
    SECTION( "input validation" ) {
        
        // no validation
        SUCCEED();
    }
    destroyQuESTEnv(env);
}



TEST_CASE( "getProbAmp", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    
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
            REQUIRE_THROWS_WITH( getProbAmp(vec,ind), Contains("Invalid amplitude index") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( getProbAmp(mat,0), Contains("valid only for state-vectors") );
            destroyQureg(mat, env);
        }
    }
    destroyQureg(vec, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "getRealAmp", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    
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
            REQUIRE_THROWS_WITH( getRealAmp(vec,ind), Contains("Invalid amplitude index") );
        }
        SECTION( "density-matrix" ) {
            
            Qureg mat = createDensityQureg(NUM_QUBITS, env);
            REQUIRE_THROWS_WITH( getRealAmp(mat,0), Contains("valid only for state-vectors") );
            destroyQureg(mat, env);
        }
    }
    destroyQureg(vec, env);
    destroyQuESTEnv(env);
}


