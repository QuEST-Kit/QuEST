// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Provides validation of user input
 */

#ifdef __cplusplus
extern "C" {
#endif
	
# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"

# include <stdio.h>
# include <stdlib.h>
		
void myTestFunc(const char* caller) {
	printf("Called from: %s\n", caller);
}

static const char* errorMessages[] = {
	[E_INVALID_NUM_QUBITS] = "Invalid number of qubits. Must create >0.",
	[E_INVALID_TARGET_QUBIT] = "Invalid target qubit. Note qubits are zero indexed.",
	[E_INVALID_CONTROL_QUBIT] = "Invalid control qubit. Note qubits are zero indexed.",
	[E_INVALID_STATE_INDEX] = "Invalid state index. Must be >=0 and <2^numQubits.",
	[E_TARGET_IS_CONTROL] = "Control qubit cannot equal target qubit.",
	[E_TARGET_IN_CONTROLS] = "Control qubits cannot include target qubit.",
	[E_INVALID_NUM_CONTROLS] = "Invalid number of control qubits. Must be >0 and <numQubits.",
	[E_NON_UNITARY_MATRIX] = "Matrix is not unitary.",
	[E_NON_UNITARY_COMPACT] = "Compact matrix is not unitary.",
	[E_SYS_TOO_BIG_TO_PRINT] = "Invalid system size. Cannot print output for systems greater than 5 qubits.",
	[E_COLLAPSE_STATE_ZERO_PROB] = "Can't collapse to state with zero probability.",
	[E_INVALID_QUBIT_OUTCOME] = "Invalid measurement outcome -- must be either 0 or 1.",
	[E_CANNOT_OPEN_FILE] = "Could not open file",
	[E_SECOND_ARG_MUST_BE_STATEVEC] = "Second argument must be a state-vector.",
	[E_MISMATCHING_REGISTER_DIMENSIONS] = "Dimensions of the qubit registers don't match.",
	[E_DEFINED_ONLY_FOR_STATEVECS] = "Valid only for state-vectors.",
	[E_DEFINED_ONLY_FOR_DENSMATRS] = "Valid only for density matrices."
};

/*
 OLD CODES: MIGHT BE NEEDED IN REFACTORING OLD VALIDATION

const char* errorCodes[] = {
    "Success",                                              // 0
    "Invalid target qubit. Note qubits are zero indexed.",  // 1
    "Invalid control qubit. Note qubits are zero indexed.", // 2 
    "Control qubit cannot equal target qubit.",             // 3
    "Invalid number of control qubits",                     // 4
    "Invalid unitary matrix.",                              // 5
    "as" // wot, this is actually compact non-initary       // 6
    "Invalid system size. Cannot print output for systems greater than 5 qubits.", // 7
    "Can't collapse to state with zero probability.", 		// 8
    "Invalid number of qubits.", 							// 9
    "Invalid measurement outcome -- must be either 0 or 1.",// 10
    "Could not open file.",									// 11
	"Second argument must be a pure state, not a density matrix.", // 12
	"Dimensions of the qubit registers do not match.", 		// 13
	"This operation is only defined for density matrices.",	// 14
	"This operation is only defined for two pure states.",	// 15
	"An non-unitary internal operation (phaseShift) occured.", //16
};
*/

void exitWithError(ErrorCode code, const char* func){
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, errorMessages[code]);
    printf("!!!\n");
    printf("exiting..\n");
    exit(code);
}

void QuESTAssert(int isValid, ErrorCode code, const char* func){
    if (!isValid) exitWithError(code, func);
}

int validateAlphaBeta(Complex alpha, Complex beta){
    if ( absReal(alpha.real*alpha.real 
                + alpha.imag*alpha.imag
                + beta.real*beta.real 
                + beta.imag*beta.imag - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateUnitVector(REAL ux, REAL uy, REAL uz){
    if ( absReal(sqrt(ux*ux + uy*uy + uz*uz) - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateMatrixIsUnitary(ComplexMatrix2 u){

    if ( absReal(u.r0c0.real*u.r0c0.real 
                + u.r0c0.imag*u.r0c0.imag
                + u.r1c0.real*u.r1c0.real
                + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    if ( absReal(u.r0c1.real*u.r0c1.real 
                + u.r0c1.imag*u.r0c1.imag
                + u.r1c1.real*u.r1c1.real
                + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;
    if ( absReal(u.r0c0.real*u.r0c1.real 
                + u.r0c0.imag*u.r0c1.imag
                + u.r1c0.real*u.r1c1.real
                + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;
    if ( absReal(u.r0c1.real*u.r0c0.imag
                - u.r0c0.real*u.r0c1.imag
                + u.r1c1.real*u.r1c0.imag
                - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;
    return 1;
}


#ifdef __cplusplus
}
#endif