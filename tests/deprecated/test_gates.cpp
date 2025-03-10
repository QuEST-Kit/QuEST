/** @file
 * Ported tests of the deprecated QuEST v3 interface,
 * unit testing the "gates" module.
 * 
 * This file should be excluded from doxygen parsing so 
 * as not to conflict with the doc of the v4 unit tests.
 * 
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */

// must include catch.hpp first, because quest's
// deprecation library will include a definition
// of Vector{}
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


/** @sa collapseToOutcome
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "collapseToOutcome", "[gates]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
    
    SECTION( "correctness" ) {
        
        int qubit = GENERATE( range(0,NUM_QUBITS) );
        int outcome = GENERATE( 0, 1 );
        
        // repeat these random tests 10 times on every qubit, and for both outcomes
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            // use a random L2 state for every qubit & outcome
            QVector vecRef = getRandomStateVector(NUM_QUBITS);
            toQureg(vec, vecRef);
            
            // calculate prob of outcome
            qreal prob = 0;
            for (size_t ind=0; ind<vecRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // target-th bit
                if (bit == outcome)
                    prob += pow(abs(vecRef[ind]), 2);
            }
                
            // renormalise by the outcome prob
            for (size_t ind=0; ind<vecRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // target-th bit
                if (bit == outcome)
                    vecRef[ind] /= sqrt(prob);
                else 
                    vecRef[ind] = 0;
            }
            
            qreal res = collapseToOutcome(vec, qubit, outcome);
            REQUIRE( res == Approx(prob) );
            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            // use a random density matrix for every qubit & outcome 
            QMatrix matRef = getRandomDensityMatrix(NUM_QUBITS);
            toQureg(mat, matRef);
            
            // prob is sum of diagonal amps (should be real) where target bit is outcome
            qcomp tr = 0;
            for (size_t ind=0; ind<matRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // qubit-th bit
                if (bit == outcome)
                    tr += matRef[ind][ind];
            }
            REQUIRE( imag(tr) == Approx(0).margin(REAL_EPS) );
            qreal prob = real(tr);
            
            // renorm (/prob) every |*outcome*><*outcome*| state, zeroing all others 
            for (size_t r=0; r<matRef.size(); r++) {
                for (size_t c=0; c<matRef.size(); c++) {
                    int ketBit = (c >> qubit) & 1;
                    int braBit = (r >> qubit) & 1;
                    
                    if (ketBit == outcome && braBit == outcome)
                        matRef[r][c] /= prob;
                    else
                        matRef[r][c] = 0;
                }
            }
            
            qreal res = collapseToOutcome(mat, qubit, outcome);
            REQUIRE( res == Approx(prob) );
            REQUIRE( areEqual(mat, matRef) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
            int qubit = GENERATE( -1, NUM_QUBITS );
            int outcome = 0;
            REQUIRE_THROWS_WITH( collapseToOutcome(mat, qubit, outcome), ContainsSubstring("Invalid target qubit") );
        }
        SECTION( "outcome value" ) {
            
            int qubit = 0;
            int outcome = GENERATE( -1, 2 );
            REQUIRE_THROWS_WITH( collapseToOutcome(mat, qubit, outcome), ContainsSubstring("qubit measurement outcome") && ContainsSubstring("invalid") );
        }
        SECTION( "outcome probability" ) {
            
            initZeroState(vec);
            REQUIRE_THROWS_WITH( collapseToOutcome(vec, 0, 1), ContainsSubstring("impossibly unlikely") );
            initClassicalState(vec, 1);
            REQUIRE_THROWS_WITH( collapseToOutcome(vec, 0, 0), ContainsSubstring("impossibly unlikely") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa measure
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "measure", "[gates]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
        
    SECTION( "correctness" ) {
        
        int qubit = GENERATE( range(0,NUM_QUBITS) );
        
        // repeat these random tests 10 times on every qubit 
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            QVector vecRef = getRandomStateVector(NUM_QUBITS);
            toQureg(vec, vecRef);
            
            int outcome = measure(vec, qubit);
            REQUIRE( (outcome == 0 || outcome == 1) );
            
            // calculate prob of this outcome
            qreal prob = 0;
            for (size_t ind=0; ind<vecRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // target-th bit
                if (bit == outcome)
                    prob += pow(abs(vecRef[ind]), 2);
            }
                
            REQUIRE( prob > REAL_EPS );
                
            // renormalise by the outcome prob
            for (size_t ind=0; ind<vecRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // target-th bit
                if (bit == outcome)
                    vecRef[ind] /= sqrt(prob);
                else 
                    vecRef[ind] = 0;
            }
            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            QMatrix matRef = getRandomDensityMatrix(NUM_QUBITS);
            toQureg(mat, matRef);
            
            int outcome = measure(mat, qubit);
            REQUIRE( (outcome == 0 || outcome == 1) );
            
            // compute prob of this outcome
            qreal prob = 0;
            for (size_t ind=0; ind<matRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // qubit-th bit
                if (bit == outcome)
                    prob += real(matRef[ind][ind]);
            }
            
            REQUIRE( prob > REAL_EPS );
                        
            // renorm (/prob) every |*outcome*><*outcome*| state, zeroing all others 
            for (size_t r=0; r<matRef.size(); r++) {
                for (size_t c=0; c<matRef.size(); c++) {
                    int ketBit = (c >> qubit) & 1;
                    int braBit = (r >> qubit) & 1;
                    
                    if (ketBit == outcome && braBit == outcome)
                        matRef[r][c] /= prob;
                    else
                        matRef[r][c] = 0;
                }
            }
            
            REQUIRE( areEqual(mat, matRef) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
            int qubit = GENERATE( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( measure(vec, qubit), ContainsSubstring("Invalid target qubit") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}



/** @sa measureWithStats
 * @ingroup deprecatedtests 
 * @author Tyson Jones 
 */
TEST_CASE( "measureWithStats", "[gates]" ) {
    
    Qureg vec = createForcedQureg(NUM_QUBITS);
    Qureg mat = createForcedDensityQureg(NUM_QUBITS);
        
    SECTION( "correctness" ) {
        
        int qubit = GENERATE( range(0,NUM_QUBITS) );
        
        // repeat these random tests 10 times on every qubit 
        GENERATE( range(0,10) );
        
        SECTION( "state-vector" ) {
            
            QVector vecRef = getRandomStateVector(NUM_QUBITS);
            toQureg(vec, vecRef);
            
            qreal res;
            int outcome = measureWithStats(vec, qubit, &res);
            REQUIRE( (outcome == 0 || outcome == 1) );
            
            // calculate prob of this outcome
            qreal prob = 0;
            for (size_t ind=0; ind<vecRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // target-th bit
                if (bit == outcome)
                    prob += pow(abs(vecRef[ind]), 2);
            }
            
            REQUIRE( prob == Approx(res) );
            
            // renormalise by the outcome prob
            for (size_t ind=0; ind<vecRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // target-th bit
                if (bit == outcome)
                    vecRef[ind] /= sqrt(prob);
                else 
                    vecRef[ind] = 0;
            }
            REQUIRE( areEqual(vec, vecRef) );
        }
        SECTION( "density-matrix" ) {
            
            QMatrix matRef = getRandomDensityMatrix(NUM_QUBITS);
            toQureg(mat, matRef);
            
            qreal res;
            int outcome = measureWithStats(mat, qubit, &res);
            REQUIRE( (outcome == 0 || outcome == 1) );
            
            // compute prob of this outcome
            qreal prob = 0;
            for (size_t ind=0; ind<matRef.size(); ind++) {
                int bit = (ind >> qubit) & 1; // qubit-th bit
                if (bit == outcome)
                    prob += real(matRef[ind][ind]);
            }
            
            REQUIRE( prob == Approx(res) );
                        
            // renorm (/prob) every |*outcome*><*outcome*| state, zeroing all others 
            for (size_t r=0; r<matRef.size(); r++) {
                for (size_t c=0; c<matRef.size(); c++) {
                    int ketBit = (c >> qubit) & 1;
                    int braBit = (r >> qubit) & 1;
                    
                    if (ketBit == outcome && braBit == outcome)
                        matRef[r][c] /= prob;
                    else
                        matRef[r][c] = 0;
                }
            }
            
            REQUIRE( areEqual(mat, matRef) );
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
            int qubit = GENERATE( -1, NUM_QUBITS );
            qreal res;
            REQUIRE_THROWS_WITH( measureWithStats(vec, qubit, &res), ContainsSubstring("Invalid target qubit") );
        }
    }
    destroyQureg(vec);
    destroyQureg(mat);
}