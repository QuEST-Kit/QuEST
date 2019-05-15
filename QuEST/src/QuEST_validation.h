// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Provides validation defined in QuEST_validation.c which is used exclusively by QuEST.c
 */
 
# ifndef QUEST_VALIDATION_H
# define QUEST_VALIDATION_H

# include "QuEST.h"

# ifdef __cplusplus
extern "C" {
# endif

void validateCreateNumQubits(int numQubits, const char* caller);

void validateStateIndex(Qureg qureg, long long int stateInd, const char* caller);

void validateTarget(Qureg qureg, int targetQubit, const char* caller);

void validateControlTarget(Qureg qureg, int controlQubit, int targetQubit, const char* caller);

void validateUniqueTargets(Qureg qureg, int qubit1, int qubit2, const char* caller);

void validateMultiControls(Qureg qureg, int* controlQubits, const int numControlQubits, const char* caller);

void validateMultiControlsTarget(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, const char* caller);

void validateUnitaryMatrix(ComplexMatrix2 u, const char* caller);

void validateUnitaryComplexPair(Complex alpha, Complex beta, const char* caller);

void validateVector(Vector vector, const char* caller);

void validateStateVecQureg(Qureg qureg, const char* caller);

void validateDensityMatrQureg(Qureg qureg, const char* caller);

void validateOutcome(int outcome, const char* caller);

void validateMeasurementProb(qreal prob, const char* caller);

void validateMatchingQuregDims(Qureg qureg1, Qureg qureg2, const char *caller);

void validateMatchingQuregTypes(Qureg qureg1, Qureg qureg2, const char *caller);

void validateSecondQuregStateVec(Qureg qureg2, const char *caller);

void validateNumAmps(Qureg qureg, long long int startInd, long long int numAmps, const char* caller);

void validateFileOpened(int opened, const char* caller);

void validateProb(qreal prob, const char* caller);

void validateNormProbs(qreal prob1, qreal prob2, const char* caller);

void validateOneQubitDephaseProb(qreal prob, const char* caller);

void validateTwoQubitDephaseProb(qreal prob, const char* caller);

void validateOneQubitDepolProb(qreal prob, const char* caller);

void validateTwoQubitDepolProb(qreal prob, const char* caller);

void validateOneQubitDampingProb(qreal prob, const char* caller);


# ifdef __cplusplus
}
# endif

# endif // QUEST_VALIDATION_H
