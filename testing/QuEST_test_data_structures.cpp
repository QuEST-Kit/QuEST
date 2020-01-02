
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



TEST_CASE( "fromComplex", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "getStaticComplexMatrixN", "[data_structures]" ) {
    
    /* use of this function may be illegal in C++ */
    FAIL();
}



TEST_CASE( "toComplex", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "createCloneQureg", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "createComplexMatrixN", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "createDensityQureg", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "createQuESTEnv", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "createQureg", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "destroyComplexMatrixN", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "destroyQuESTEnv", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "destroyQureg", "[data_structures]" ) {
    FAIL();
}



TEST_CASE( "initComplexMatrixN", "[data_structures]" ) {
    
    /* use of this function may be illegal in C++ */
    FAIL();
}