
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


TEST_CASE( "calcDensityInnerProduct", "[calculations]" ) {

    FAIL();
    
    QuESTEnv env = createQuESTEnv();
    
    SECTION( "correctness" ) {
        
        SECTION( "density-matrix" ) {
            
        }
    }
    SECTION( "input validation" ) {
        
    }
    destroyQuESTEnv(env);
}



TEST_CASE( "calcExpecPauliProd", "[calculations]" ) {
    
    FAIL();
    
}



TEST_CASE( "calcExpecPauliSum", "[calculations]" ) {
    
    FAIL();
    
}



TEST_CASE( "calcFidelity", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    Qureg pure = createQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        // repeat the below random tests 10 times 
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            /* calcFidelity computes |<vec|pure>|^2 */
             
            SECTION( "normalised" ) {
                 
                // random L2 vectors
                QVector vecRef = getNormalised(getRandomQVector(1<<NUM_QUBITS));
                QVector pureRef = getNormalised(getRandomQVector(1<<NUM_QUBITS));
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
                
                QVector pureRef = getNormalised(getRandomQVector(1<<NUM_QUBITS));
                toQureg(pure, pureRef);
                
                // test when density matrix is the same pure state 
                QMatrix matRef = getKetBra(pureRef, pureRef);
                toQureg(mat, matRef);
                REQUIRE( calcFidelity(mat,pure) == Approx(1) ); 
                
                // test when density matrix is a random pure state
                QVector r1 = getNormalised(getRandomQVector(1<<NUM_QUBITS));
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
                
                QVector pureRef = getNormalised(getRandomQVector(1<<NUM_QUBITS));
                toQureg(pure, pureRef);
            
                // test when density matrix is mixed 
                QVector r1 = getNormalised(getRandomQVector(1<<NUM_QUBITS));
                QVector r2 = getNormalised(getRandomQVector(1<<NUM_QUBITS));
                QMatrix matRef = .3*getKetBra(r1, r1) + .7*getKetBra(r2, r2); // .3|r1><r1| + .7|r2><r2|
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
    SECTION( "validation" ) {
        
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
    destroyQureg(vec, env);
    destroyQureg(mat, env);
    destroyQureg(pure, env);
    destroyQuESTEnv(env);
}



TEST_CASE( "calcHilbertSchmidtDistance", "[calculations]" ) {
    
    FAIL();
    
}



TEST_CASE( "calcInnerProduct", "[calculations]" ) {
    
    FAIL();
    
}



TEST_CASE( "calcProbOfOutcome", "[calculations]" ) {
    
    FAIL();
    
}



TEST_CASE( "calcPurity", "[calculations]" ) {
    
    QuESTEnv env = createQuESTEnv();
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        // perform the following random tests 10 times 
        GENERATE( range(1,10) );
        
        SECTION( "density-matrix" ) {
            
            // pure states have unity purity 
            initZeroState(mat);
            REQUIRE( calcPurity(mat) == 1 );
            
            // (try also a pure random L2-vector)
            QVector r1 = getNormalised(getRandomQVector(1<<NUM_QUBITS)); // |r>
            QMatrix m1 = getKetBra(r1, r1); // |r><r|
            toQureg(mat, m1);
            REQUIRE( calcPurity(mat) == Approx(1) );
            
            // mixed states have 0 < purity < 1, given by Tr(rho^2)
            QVector r2 = getNormalised(getRandomQVector(1<<NUM_QUBITS));
            QMatrix m2 = getKetBra(r2, r2);
            toQureg(mat, .3*m1 + .7*m2); // mix: .3|r1><r1| + .7|r2><r2|
            qreal purity = calcPurity(mat);
            REQUIRE( purity < 1 );
            REQUIRE( purity >= 1/pow(2.,NUM_QUBITS) );
            
            // (compare to Tr(mat^2))
            QMatrix ref = toQMatrix(mat);
            QMatrix prod = ref*ref;
            qreal tr = 0;
            for (size_t i=0; i<prod.size(); i++)
                tr += real(prod[i][i]);
            REQUIRE( purity == Approx(tr) );
            
            // unphysical states give sum_{ij} |rho_ij|^2
            ref = getRandomQMatrix(ref.size());
            qreal tot = 0;
            for (size_t i=0; i<ref.size(); i++)
                for (size_t j=0; j<ref.size(); j++)
                    tot += pow(abs(ref[i][j]), 2);
                
            toQureg(mat, ref);
            REQUIRE( calcPurity(mat) == Approx(tot) );
        }
    }
    SECTION( "validation" ) {
        
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
    
    QuESTEnv env = createQuESTEnv();
    Qureg vec = createQureg(NUM_QUBITS, env);
    Qureg mat = createDensityQureg(NUM_QUBITS, env);
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
            // normalised: prob(vec) = 1
            initPlusState(vec);
            REQUIRE( calcTotalProb(vec) == Approx(1) );
            
            // zero norm: prob(vec) = 0
            initBlankState(vec);
            REQUIRE( calcTotalProb(vec) == 0 );
            
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
    destroyQureg(vec, env);
    destroyQureg(mat, env);
    destroyQuESTEnv(env);
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


