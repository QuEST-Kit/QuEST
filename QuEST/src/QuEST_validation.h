// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Provides validation defined in QuEST_validation.c which is used exclusively by QuEST.c
 *
 * @author Tyson Jones
 * @author Ania Brown (original testing of qubit indices, unitarity, valid collapse probability)
 * @author Balint Koczor (Kraus maps)
 */
 
# ifndef QUEST_VALIDATION_H
# define QUEST_VALIDATION_H

# include "QuEST.h"
# include <stdio.h>

# ifdef __cplusplus
extern "C" {
# endif

void validateNumRanks(int numRanks, const char* caller);

void validateNumQubitsInQureg(int numQubits, int numRanks, const char* caller);

void validateNumQubitsInMatrix(int numQubits, const char* caller);

void validateNumQubitsInDiagOp(int numQubits, int numRanks, const char* caller);

void validateAmpIndex(Qureg qureg, long long int ampInd, const char* caller);

void validateStateIndex(Qureg qureg, long long int stateInd, const char* caller);

void validateTarget(Qureg qureg, int targetQubit, const char* caller);

void validateControlTarget(Qureg qureg, int controlQubit, int targetQubit, const char* caller);

void validateUniqueTargets(Qureg qureg, int qubit1, int qubit2, const char* caller);

void validateMultiQubits(Qureg qureg, int* qubits, int numQubits, const char* caller);

void validateMultiTargets(Qureg qurge, int* targetQubits, int numTargetQubits, const char* caller);

void validateMultiControls(Qureg qureg, int* controlQubits, int numControlQubits, const char* caller);

void validateMultiControlsTarget(Qureg qureg, int* controlQubits, int numControlQubits, int targetQubit, const char* caller);

void validateMultiControlsMultiTargets(Qureg qureg, int* controlQubits, int numControlQubits, int* targetQubits, int numTargetQubits, const char* caller);

void validateControlState(int* controlState, int numControlQubits, const char* caller);

void validateOneQubitUnitaryMatrix(ComplexMatrix2 u, const char* caller);

void validateTwoQubitUnitaryMatrix(Qureg qureg, ComplexMatrix4 u, const char* caller);

void validateMultiQubitMatrix(Qureg qureg, ComplexMatrixN u, int numTargs, const char* caller);

void validateMultiQubitUnitaryMatrix(Qureg qureg, ComplexMatrixN u, int numTargs, const char* caller);

void validateMultiQubitMatrixFitsInNode(Qureg qureg, int numTargets, const char* caller);

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

void validateFileOpened(int opened, char* fn, const char* caller);

void validateProb(qreal prob, const char* caller);

void validateNormProbs(qreal prob1, qreal prob2, const char* caller);

void validateOneQubitDephaseProb(qreal prob, const char* caller);

void validateTwoQubitDephaseProb(qreal prob, const char* caller);

void validateOneQubitDepolProb(qreal prob, const char* caller);

void validateTwoQubitDepolProb(qreal prob, const char* caller);

void validateOneQubitPauliProbs(qreal probX, qreal probY, qreal probZ, const char* caller);

void validatePauliCodes(enum pauliOpType* pauliCodes, int numPauliCodes, const char* caller);

void validateNumPauliSumTerms(int numTerms, const char* caller);

void validateMatrixInit(ComplexMatrixN matr, const char* caller);

void validateOneQubitKrausMap(Qureg qureg, ComplexMatrix2* ops, int numOps, const char* caller);

void validateTwoQubitKrausMap(Qureg qureg, ComplexMatrix4* ops, int numOps, const char* caller);

void validateMultiQubitKrausMap(Qureg qureg, int numTargs, ComplexMatrixN* ops, int numOps, const char* caller);

void validateOneQubitDampingProb(qreal prob, const char* caller);

void validateHamilParams(int numQubits, int numTerms, const char* caller);

void validatePauliHamil(PauliHamil hamil, const char* caller);

void validateMatchingQuregPauliHamilDims(Qureg qureg, PauliHamil hamil, const char* caller);

void validateHamilFileParams(int numQubits, int numTerms, FILE* file, char* fn, const char* caller);

void validateHamilFileCoeffParsed(int parsed, PauliHamil h, FILE* file, char* fn, const char* caller);

void validateHamilFilePauliParsed(int parsed, PauliHamil h, FILE* file, char* fn, const char* caller);

void validateHamilFilePauliCode(enum pauliOpType code, PauliHamil h, FILE* file, char* fn, const char* caller);

void validateTrotterParams(int order, int reps, const char* caller);

void validateDiagOpInit(DiagonalOp, const char* caller);

void validateDiagonalOp(Qureg qureg, DiagonalOp op, const char* caller);

void validateDiagPauliHamil(DiagonalOp op, PauliHamil hamil, const char *caller);

void validateDiagPauliHamilFromFile(PauliHamil h, int numRanks, const char* caller);

void validateNumElems(DiagonalOp op, long long int startInd, long long int numElems, const char* caller);

void validateQubitSubregs(Qureg qureg, int* qubits, int* numQubitsPerReg, const int numReg, const char* caller);

void validatePhaseFuncOverrides(const int numQubits, enum bitEncoding encoding, long long int* overrideInds, int numOverrides, const char* caller);

void validateMultiVarPhaseFuncOverrides(int* numQubitsPerReg, const int numRegs, enum bitEncoding encoding, long long int* overrideInds, int numOverrides, const char* caller);

void validatePhaseFuncTerms(int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, int numOverrides, const char* caller);

void validateMultiVarPhaseFuncTerms(int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* exponents, int* numTermsPerReg, const char* caller);

void validatePhaseFuncName(enum phaseFunc funcCode, int numRegs, int numParams, const char* caller);

void validateBitEncoding(int numQubits, enum bitEncoding encoding, const char* caller);

void validateMultiRegBitEncoding(int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, const char* caller);

# ifdef __cplusplus
}
# endif

# endif // QUEST_VALIDATION_H
