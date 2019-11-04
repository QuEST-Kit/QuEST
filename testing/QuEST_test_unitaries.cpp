/** @file
 * Unit testing for QuEST's 'unitaries' API. 
 * The tests are in alphabetical order of the API doc. 
 *
 * These tests work by constructing, from the unitary specification (e.g. 
 * control and target qubits), a full-Hilbert space complex matrix. This is 
 * then multiplied onto statevectors, or multiplied and it's conjugate-transpose
 * right-multiplied onto density matrices.
 *
 * QuEST's user validation handling is unit tested by redefining exitWithError 
 * (a weak C symbol) to throw a C++ exception, caught by the Catch2 library.
 *
 * @author Tyson Jones
 */
 
#include "catch.hpp"
#include "QuEST.h"
#include "QuEST_test_utils.hpp"

using Catch::Matchers::Contains;

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

/** Redefinition of QuEST_validation's invalidQuESTInputError function, called when a 
 * user passes an incorrect parameter (e.g. an negative qubit index). This is 
 * redefined here to, in lieu of printing and exiting, throw a C++ exception
 * which can be caught (and hence unit tested for) by Catch2
 */
 extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
     throw errMsg;
 } 

/** Prepares the needed data structures for unit testing. This creates 
 * the QuEST environment, a statevector and density matrix of the size 'numQb',
 * and corresponding QVector and QMatrix instances for analytic comparison.
 * numQb should be NUM_QUBITS unless motivated otherwise.
 */
#define PREPARE_TEST(env, quregVec, quregMatr, refVec, refMatr, numQb) \
    QuESTEnv env = createQuESTEnv(); \
    Qureg quregVec = createQureg(numQb, env); \
    Qureg quregMatr = createDensityQureg(numQb, env); \
    initDebugState(quregVec); \
    initDebugState(quregMatr); \
    QVector refVec = toQVector(quregVec); \
    QMatrix refMatr = toQMatrix(quregMatr);

/** Destroys the data structures made by PREPARE_TEST */
#define CLEANUP_TEST(env, quregVec, quregMatr) \
    destroyQureg(quregVec, env); \
    destroyQureg(quregMatr, env); \
    destroyQuESTEnv(env);
    
    
    
TEST_CASE( "pauliX", "[unitaries]" ) {
    
    PREPARE_TEST(env, quregVec, quregMatr, refVec, refMatr, NUM_QUBITS);

    QMatrix op{{0,1},{1,0}};
    
    SECTION( "state-vector correctness" ) {

        int target = GENERATE_COPY( range(0,NUM_QUBITS) );
        pauliX(quregVec, target);
        applyUnitaryOp(refVec, target, op);
        REQUIRE( areEqual(quregVec, refVec) );
    }
    
    SECTION( "density-matrix correctness" ) {
        
        int target = GENERATE_COPY( range(0,NUM_QUBITS) );
        pauliX(quregMatr, target);
        applyUnitaryOp(refMatr, target, op);
        REQUIRE( areEqual(quregMatr, refMatr) );
    }
    
    SECTION( "input validation" ) {
                
        SECTION( "qubit indices" ) {
            
            REQUIRE_THROWS_WITH( pauliX(quregVec, -1), Contains("Invalid target") );
            REQUIRE_THROWS_WITH( pauliX(quregVec, NUM_QUBITS), Contains("Invalid target") );
        }
    }
    
    CLEANUP_TEST(env, quregVec, quregMatr);
}



TEST_CASE( "compactUnitary", "[unitaries]" ) {
    
    PREPARE_TEST(env, quregVec, quregMatr, refVec, refMatr, NUM_QUBITS);
    
    qcomp a = .3 * exp(2i);
    qcomp b = sqrt(1-abs(a)*abs(a)) * exp(-3i);
    Complex alpha = toComplex( a );
    Complex beta = toComplex( b );
    QMatrix op = toQMatrix(alpha, beta);
    
    SECTION( "state-vector correctness" ) {
        
        int target = GENERATE_COPY( range(0,NUM_QUBITS) );
        compactUnitary(quregVec, target, alpha, beta);
        applyUnitaryOp(refVec, target, op);    
        REQUIRE( areEqual(quregVec, refVec) );
    }
    
    SECTION( "density-matrix correctness" ) {
        
        int target = GENERATE_COPY( range(0,NUM_QUBITS) );
        compactUnitary(quregMatr, target, alpha, beta);
        applyUnitaryOp(refMatr, target, op);    
        REQUIRE( areEqual(quregMatr, refMatr) );
    }
    
    SECTION( "input validation" ) {
        
        SECTION( "qubit indices" ) {
            
            REQUIRE_THROWS_WITH( compactUnitary(quregVec, -1, alpha, beta), Contains("Invalid target") );
            REQUIRE_THROWS_WITH( compactUnitary(quregVec, NUM_QUBITS, alpha, beta), Contains("Invalid target") );
        }
        
        SECTION( "unitarity" ) {
            
            // unitary when |a|^2 + |b^2 = 1
            alpha = {.real=1, .imag=2};
            beta = {.real=3, .imag=4};
            REQUIRE_THROWS_WITH( compactUnitary(quregVec, 0, alpha, beta), Contains("unitary") );
        }
    }
        
    CLEANUP_TEST(env, quregVec, quregMatr);
}



TEST_CASE( "controlledCompactUnitary", "[unitaries]" ) {
    
    PREPARE_TEST(env, quregVec, quregMatr, refVec, refMatr, NUM_QUBITS);
    
    qcomp a = .8 * exp(-1.5i);
    qcomp b = sqrt(1-abs(a)*abs(a)) * exp(.3i);
    Complex alpha = toComplex( a );
    Complex beta = toComplex( b );
    QMatrix op = toQMatrix(alpha, beta);
    
    SECTION( "state-vector correctness" ) {
        
        int target = GENERATE_COPY( range(0,NUM_QUBITS) );
        int control = GENERATE_COPY( filter([=](int c){ return c!=target; }, range(0,NUM_QUBITS)) );
        controlledCompactUnitary(quregVec, control, target, alpha, beta);
        applyUnitaryOp(refVec, control, target, op);    
        REQUIRE( areEqual(quregVec, refVec) );
    }
    
    SECTION( "density-matrix correctness" ) {
        
        int target = GENERATE_COPY( range(0,NUM_QUBITS) );
        int control = GENERATE_COPY( filter([=](int c){ return c!=target; }, range(0,NUM_QUBITS)) );
        controlledCompactUnitary(quregMatr, control, target, alpha, beta);
        applyUnitaryOp(refMatr, control, target, op);    
        REQUIRE( areEqual(quregMatr, refMatr) );
    }
    
    SECTION( "input validation" ) {
        
        SECTION( "control and target collision" ) {
            
            REQUIRE_THROWS_WITH( controlledCompactUnitary(quregVec, 0, 0, alpha, beta), Contains("Control") && Contains("target") );
        }
        
        SECTION( "qubit indices" ) {
            
            int qb = GENERATE_COPY( -1, NUM_QUBITS );
            REQUIRE_THROWS_WITH( controlledCompactUnitary(quregVec, qb, 0, alpha, beta), Contains("Invalid control") );
            REQUIRE_THROWS_WITH( controlledCompactUnitary(quregVec, 0, qb, alpha, beta), Contains("Invalid target") );
        }
        
        SECTION( "unitarity" ) {

            // unitary when |a|^2 + |b^2 = 1
            alpha = {.real=1, .imag=2};
            beta = {.real=3, .imag=4};
            REQUIRE_THROWS_WITH( controlledCompactUnitary(quregVec, 0, 1, alpha, beta), Contains("unitary") );
        }
    }
        
    CLEANUP_TEST(env, quregVec, quregMatr);
}