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
 * Q = qureg, S = state, T = target, C = control
 */
# define validateStateIndex(Q, S) auto_validateStateIndex(Q, S, __func__)
# define validateTarget(Q, T) auto_validateTarget(Q, T, __func__)
# define validateControlTarget(Q, C, T) auto_validateControlTarget(Q, C, T, __func__)

void auto_validateStateIndex(QubitRegister qureg, long long int stateInd, const char* caller);
void auto_validateTarget(QubitRegister qureg, int targetQubit, const char* caller);
void auto_validateControlTarget(QubitRegister qureg, int controlQubit, int targetQubit, const char* caller);

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
	E_NON_UNITARY_COMPACT,
	E_SYS_TOO_BIG_TO_PRINT,
	E_COLLAPSE_STATE_ZERO_PROB,
	E_INVALID_QUBIT_OUTCOME,
	E_CANNOT_OPEN_FILE,
	E_SECOND_ARG_MUST_BE_STATEVEC,
	E_MISMATCHING_REGISTER_DIMENSIONS,
	E_DEFINED_ONLY_FOR_STATEVECS,
	E_DEFINED_ONLY_FOR_DENSMATRS,
} ErrorCode;

void exitWithError(ErrorCode code, const char *func);

void QuESTAssert(int isValid, ErrorCode code, const char *func);

// need renaming
int validateUnitComplex(Complex alpha);

int validateMatrixIsUnitary(ComplexMatrix2 u);

int validateAlphaBeta(Complex alpha, Complex beta);

int validateUnitVector(REAL ux, REAL uy, REAL uz);


# ifdef __cplusplus
}
# endif

# endif // QuEST_VALIDATION