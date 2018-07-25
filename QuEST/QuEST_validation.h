// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Provides defined in QuEST_validation.c which are used by QuEST.c
 */
 
# ifndef QuEST_VALIDATION
# define QuEST_VALIDATION

# include "QuEST.h"

# ifdef __cplusplus
extern "C" {
# endif
	
/* Wraps validation to automatically pass along the caller's signature.
 * N = number, Q = qureg, S = state, T = target, C = control, CS = controls, U = unitary
 * A = alpha, B = beta, V = vector, M = measurement outcome
 */
# define validateCreateNumQubits(N) auto_validateCreateNumQubits(N, __func__)
# define validateStateIndex(Q, S) auto_validateStateIndex(Q, S, __func__)
# define validateTarget(Q, T) auto_validateTarget(Q, T, __func__)
# define validateControlTarget(Q, C, T) auto_validateControlTarget(Q, C, T, __func__)
# define validateMultiControls(Q, CS, N) auto_validateMultiControls(Q, CS, N, __func__)
# define validateMultiControlsTarget(Q, CS, N, T) auto_validateMultiControlsTarget(Q, CS, N, T, __func__)
# define validateUnitaryMatrix(U) auto_validateUnitaryMatrix(U, __func__)
# define validateUnitaryComplexPair(A, B) auto_validateUnitaryComplexPair(A, B, __func__)
# define validateVector(V) auto_validateVector(V, __func__)
# define validateStateVecQureg(Q) auto_validateStateVecQureg(Q, __func__)
# define validatDensityMatrQureg(Q) auto_validatDensityMatrQureg(Q, __func__)
# define validateOutcome(M) auto_validateOutcome(M, __func__)

void auto_validateCreateNumQubits(int numQubits, const char* caller);
void auto_validateStateIndex(QubitRegister qureg, long long int stateInd, const char* caller);
void auto_validateTarget(QubitRegister qureg, int targetQubit, const char* caller);
void auto_validateControlTarget(QubitRegister qureg, int controlQubit, int targetQubit, const char* caller);
void auto_validateMultiControls(QubitRegister qureg, int* controlQubits, const int numControlQubits, const char* caller);
void auto_validateMultiControlsTarget(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, const char* caller);
void auto_validateUnitaryMatrix(ComplexMatrix2 u, const char* caller);
void auto_validateUnitaryComplexPair(Complex alpha, Complex beta, const char* caller);
void auto_validateVector(Vector vector, const char* caller);
void auto_validateStateVecQureg(QubitRegister qureg, const char* caller);
void auto_validatDensityMatrQureg(QubitRegister qureg, const char* caller);
void auto_validateOutcome(int outcome, const char* caller);

/* BELOW WON'T NEED TO BE EXPOSED AFTER ALL INTEGRATED INTO QuEST_validation.h! :D */
typedef enum {
	E_SUCCESS=0,
	E_INVALID_NUM_QUBITS,
	E_INVALID_TARGET_QUBIT,
	E_INVALID_CONTROL_QUBIT,
	E_INVALID_STATE_INDEX,
	E_TARGET_IS_CONTROL,
	E_TARGET_IN_CONTROLS,
	E_INVALID_NUM_CONTROLS,
	E_NON_UNITARY_MATRIX,
	E_NON_UNITARY_COMPLEX_PAIR,
	E_ZERO_VECTOR,
	E_SYS_TOO_BIG_TO_PRINT,
	E_COLLAPSE_STATE_ZERO_PROB,
	E_INVALID_QUBIT_OUTCOME,
	E_CANNOT_OPEN_FILE,
	E_SECOND_ARG_MUST_BE_STATEVEC,
	E_MISMATCHING_REGISTER_DIMENSIONS,
	E_DEFINED_ONLY_FOR_STATEVECS,
	E_DEFINED_ONLY_FOR_DENSMATRS,
} ErrorCode;

/* BELOW WON'T NEED TO BE EXPOSED AFTER ALL INTEGRATED INTO QuEST_validation.h! :D */

void exitWithError(ErrorCode code, const char *func);
void QuESTAssert(int isValid, ErrorCode code, const char *func);

/* BELOW WON'T NEED TO BE EXPOSED AFTER ALL INTEGRATED INTO QuEST_validation.h! :D */
int isComplexUnit(Complex alpha);
int isMatrixUnitary(ComplexMatrix2 u);
int isComplexPairUnitary(Complex alpha, Complex beta);
int isVectorUnit(REAL ux, REAL uy, REAL uz);


# ifdef __cplusplus
}
# endif

# endif // QuEST_VALIDATION