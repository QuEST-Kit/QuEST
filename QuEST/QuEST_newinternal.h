
# include "QuEST.h"

#ifdef __cplusplus
extern "C" {
#endif

// errors

extern const char* errorCodes[];

void exitWithError(int errorCode, const char *func);

void QuESTAssert(int isValid, int errorCode, const char *func);

unsigned long int hashString(char *str);

// validation

int validateMatrixIsUnitary(ComplexMatrix2 u);

int validateAlphaBeta(Complex alpha, Complex beta);

int validateUnitVector(REAL ux, REAL uy, REAL uz);

#ifdef __cplusplus
}
#endif