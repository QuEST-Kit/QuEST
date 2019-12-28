
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
    
    FAIL();
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
        }
        SECTION( "density-matrix" ) {
            
        }
    }
    
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
            toQureg(mat, .5*m1 + .5*m2); // mix: .5|r1><r1| + .5|r2><r2|
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


