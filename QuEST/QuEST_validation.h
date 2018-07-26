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

void validateCreateNumQubits(int numQubits, const char* caller);

void validateStateIndex(QubitRegister qureg, long long int stateInd, const char* caller);

void validateTarget(QubitRegister qureg, int targetQubit, const char* caller);

void validateControlTarget(QubitRegister qureg, int controlQubit, int targetQubit, const char* caller);

void validateMultiControls(QubitRegister qureg, int* controlQubits, const int numControlQubits, const char* caller);

void validateMultiControlsTarget(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, const char* caller);

void validateUnitaryMatrix(ComplexMatrix2 u, const char* caller);

void validateUnitaryComplexPair(Complex alpha, Complex beta, const char* caller);

void validateVector(Vector vector, const char* caller);

void validateStateVecQureg(QubitRegister qureg, const char* caller);

void validatDensityMatrQureg(QubitRegister qureg, const char* caller);

void validateOutcome(int outcome, const char* caller);

void validateMeasurementProb(REAL prob, const char* caller);

void validateMatchingQuregDims(QubitRegister qureg1, QubitRegister qureg2, const char *caller);

void validateSecondQuregStateVec(QubitRegister qureg2, const char *caller);

void validateFileOpened(int opened, const char* caller);

# ifdef __cplusplus
}
# endif

# endif // QuEST_VALIDATION