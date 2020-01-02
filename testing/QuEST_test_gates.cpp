
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



TEST_CASE( "collapseToOutcome", "[gates]" ) {
    
    FAIL();
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
        }
        SECTION( "density-matrix" ) {
            
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
        }
        SECTION( "outcome value" ) {
            
        }
        SECTION( "outcome probability" ) {
            
        }
    }
}



TEST_CASE( "measure", "[gates]" ) {
    
    FAIL();
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
        }
        SECTION( "density-matrix" ) {
            
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
        }
    }
}



TEST_CASE( "measureWithStats", "[gates]" ) {
    
    FAIL();
    
    SECTION( "correctness" ) {
        
        SECTION( "state-vector" ) {
            
        }
        SECTION( "density-matrix" ) {
            
        }
    }
    SECTION( "input validation" ) {
        
        SECTION( "qubit index" ) {
            
        }
    }
}