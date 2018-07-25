// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Functions defined in QuEST_common which are only for internal use by hardware-specific backends
 */

# ifndef QuEST_INTERNAL
# define QuEST_INTERNAL

# include "QuEST.h"

# ifdef __cplusplus
extern "C" {
# endif

typedef enum {
	E_SUCCESS=0,
	E_INVALID_NUM_QUBITS,
	E_INVALID_TARGET_QUBIT,
	E_INVALID_CONTROL_QUBIT,
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
	E_DEFINED_ONLY_FOR_DENSMATRS
} ErrorCode;

void exitWithError(ErrorCode code, const char *func);

void QuESTAssert(int isValid, ErrorCode code, const char *func);

unsigned long int hashString(char *str);

int validateUnitComplex(Complex alpha);

int validateMatrixIsUnitary(ComplexMatrix2 u);

int validateAlphaBeta(Complex alpha, Complex beta);

int validateUnitVector(REAL ux, REAL uy, REAL uz);

Complex getConjugateScalar(Complex scalar);

ComplexMatrix2 getConjugateMatrix(ComplexMatrix2 matr);

void shiftIndices(int* indices, int numIndices, int shift);

# ifdef __cplusplus
}
# endif

# endif // QuEST_INTERNAL