// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Provides validation of user input
 */

#ifdef __cplusplus
extern "C" {
#endif
	
# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"
# include "QuEST_validation.h"

# include <stdio.h>
# include <stdlib.h>
		
void auto_validateCreateNumQubits(int numQubits, const char* caller) {
	QuESTAssert(numQubits>0, E_INVALID_NUM_QUBITS, caller);
}

void auto_validateStateIndex(QubitRegister qureg, long long int stateInd, const char* caller) {
	QuESTAssert(stateInd>=0 && stateInd<qureg.numAmpsTotal, E_INVALID_STATE_INDEX, caller);
}

void auto_validateTarget(QubitRegister qureg, int targetQubit, const char* caller) {
	QuESTAssert(targetQubit>=0 && targetQubit<qureg.numQubitsRepresented, E_INVALID_TARGET_QUBIT, caller);
}

void auto_validateControl(QubitRegister qureg, int controlQubit, const char* caller) {
	QuESTAssert(controlQubit>=0 && controlQubit<qureg.numQubitsRepresented, E_INVALID_CONTROL_QUBIT, caller);
}

void auto_validateControlTarget(QubitRegister qureg, int controlQubit, int targetQubit, const char* caller) {
	auto_validateTarget(qureg, targetQubit, caller);
	auto_validateControl(qureg, controlQubit, caller);
	QuESTAssert(controlQubit != targetQubit, E_TARGET_IS_CONTROL, caller);
}

// this could reject repeated qubits but cite "too many"
void auto_validateNumControls(QubitRegister qureg, const int numControlQubits, const char* caller) {
	QuESTAssert(numControlQubits>0 && numControlQubits<=qureg.numQubitsRepresented, E_INVALID_NUM_CONTROLS, caller);
}

void auto_validateMultiControls(QubitRegister qureg, int* controlQubits, const int numControlQubits, const char* caller) {
	auto_validateNumControls(qureg, numControlQubits, caller);
	for (int i=0; i < numControlQubits; i++) {
		auto_validateControl(qureg, controlQubits[i], caller);
	}
}

void auto_validateMultiControlsTarget(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, const char* caller) {
	auto_validateTarget(qureg, targetQubit, caller);
	auto_validateMultiControls(qureg, controlQubits, numControlQubits, caller);
	for (int i=0; i < numControlQubits; i++)
		QuESTAssert(controlQubits[i] != targetQubit, E_TARGET_IN_CONTROLS, caller);
}

void auto_validateUnitaryMatrix(ComplexMatrix2 u, const char* caller) {
	QuESTAssert(isMatrixUnitary(u), E_NON_UNITARY_MATRIX, caller);
}

void auto_validateUnitaryComplexPair(Complex alpha, Complex beta, const char* caller) {
	QuESTAssert(isComplexPairUnitary(alpha, beta), E_NON_UNITARY_COMPLEX_PAIR, caller);
}

void auto_validateVector(Vector vec, const char* caller) {
	QuESTAssert(getVectorMagnitude(vec) > REAL_EPS, E_ZERO_VECTOR, caller);
}

void auto_validateStateVecQureg(QubitRegister qureg, const char* caller) {
	QuESTAssert( ! qureg.isDensityMatrix, E_DEFINED_ONLY_FOR_STATEVECS, caller);
}

void auto_validatDensityMatrQureg(QubitRegister qureg, const char* caller) {
	QuESTAssert(qureg.isDensityMatrix, E_DEFINED_ONLY_FOR_DENSMATRS, caller);
}

void auto_validateOutcome(int outcome, const char* caller) {
	QuESTAssert(outcome==0 || outcome==1, E_INVALID_QUBIT_OUTCOME, caller);
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
	[E_NON_UNITARY_COMPLEX_PAIR] = "Compact matrix formed by given complex numbers is not unitary.",
	[E_ZERO_VECTOR] = "Invalid axis vector. Must be non-zero.",
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

int isComplexUnit(Complex alpha) {
	return (absReal(1 - sqrt(alpha.real*alpha.real + alpha.imag*alpha.imag)) < REAL_EPS); 
}

/* this is only being used internally, which is a mistrust of our own code */
/* we normalise the vectors given by the user, which aren't requiredly unit */
int isVectorUnit(REAL ux, REAL uy, REAL uz) {
    return (absReal(1 - sqrt(ux*ux + uy*uy + uz*uz)) < REAL_EPS );
}

int isComplexPairUnitary(Complex alpha, Complex beta) {
    return ( absReal( -1
				+ alpha.real*alpha.real 
                + alpha.imag*alpha.imag
                + beta.real*beta.real 
                + beta.imag*beta.imag) < REAL_EPS );
}

int isMatrixUnitary(ComplexMatrix2 u) {
    if ( absReal( u.r0c0.real*u.r0c0.real 
                + u.r0c0.imag*u.r0c0.imag
                + u.r1c0.real*u.r1c0.real
                + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    if ( absReal( u.r0c1.real*u.r0c1.real 
                + u.r0c1.imag*u.r0c1.imag
                + u.r1c1.real*u.r1c1.real
                + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;
    if ( absReal( u.r0c0.real*u.r0c1.real 
                + u.r0c0.imag*u.r0c1.imag
                + u.r1c0.real*u.r1c1.real
                + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;
    if ( absReal( u.r0c1.real*u.r0c0.imag
                - u.r0c0.real*u.r0c1.imag
                + u.r1c1.real*u.r1c0.imag
                - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;
    return 1;
}



#ifdef __cplusplus
}
#endif