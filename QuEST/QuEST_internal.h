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

unsigned long int hashString(char *str);

Complex getConjugateScalar(Complex scalar);

ComplexMatrix2 getConjugateMatrix(ComplexMatrix2 matr);

void shiftIndices(int* indices, int numIndices, int shift);

# ifdef __cplusplus
}
# endif

# endif // QuEST_INTERNAL