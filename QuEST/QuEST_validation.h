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
 * A = alpha, B = beta, V = vector, M = measurement outcome, P = probability
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
# define validateDensityMatrQureg(Q) auto_validatDensityMatrQureg(Q, __func__)
# define validateOutcome(M) auto_validateOutcome(M, __func__)
# define validateMeasurementProb(P) auto_validateMeasurementProb(P, __func__)
# define validateMatchingQuregDims(Q1, Q2) auto_validateMatchingQuregDims(Q1, Q2, __func__)
# define validateSecondQuregStateVec(Q2) auto_validateSecondQuregStateVec(Q2, __func__)

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
void auto_validateMeasurementProb(REAL prob, const char* caller);
void auto_validateMatchingQuregDims(QubitRegister qureg1, QubitRegister qureg2, const char *caller);
void auto_validateSecondQuregStateVec(QubitRegister qureg2, const char *caller);

# ifdef __cplusplus
}
# endif

# endif // QuEST_VALIDATION