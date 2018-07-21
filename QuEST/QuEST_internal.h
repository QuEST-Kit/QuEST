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

extern const char* errorCodes[];

void exitWithError(int errorCode, const char *func);

void QuESTAssert(int isValid, int errorCode, const char *func);

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