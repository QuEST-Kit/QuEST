// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Provides validation of user input
 */

#ifdef __cplusplus
extern "C" {
#endif
    
# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"
# include "QuEST_validation.h"

# include <stdio.h>
# include <stdlib.h>

typedef enum {
    E_SUCCESS=0,
    E_INVALID_NUM_QUBITS,
    E_INVALID_TARGET_QUBIT,
    E_INVALID_CONTROL_QUBIT,
    E_INVALID_STATE_INDEX,
    E_INVALID_NUM_AMPS,
    E_INVALID_OFFSET_NUM_AMPS,
    E_TARGET_IS_CONTROL,
    E_TARGET_IN_CONTROLS,
    E_TARGETS_NOT_UNIQUE,
    E_INVALID_NUM_CONTROLS,
    E_NON_UNITARY_MATRIX,
    E_NON_UNITARY_COMPLEX_PAIR,
    E_ZERO_VECTOR,
    E_SYS_TOO_BIG_TO_PRINT,
    E_COLLAPSE_STATE_ZERO_PROB,
    E_INVALID_QUBIT_OUTCOME,
    E_CANNOT_OPEN_FILE,
    E_SECOND_ARG_MUST_BE_STATEVEC,
    E_MISMATCHING_QUREG_DIMENSIONS,
    E_MISMATCHING_QUREG_TYPES,
    E_DEFINED_ONLY_FOR_STATEVECS,
    E_DEFINED_ONLY_FOR_DENSMATRS,
    E_INVALID_PROB,
    E_UNNORM_PROBS,
    E_INVALID_ONE_QUBIT_DEPHASE_PROB,
    E_INVALID_TWO_QUBIT_DEPHASE_PROB,
    E_INVALID_ONE_QUBIT_DEPOL_PROB,
    E_INVALID_TWO_QUBIT_DEPOL_PROB
} ErrorCode;

static const char* errorMessages[] = {
    [E_INVALID_NUM_QUBITS] = "Invalid number of qubits. Must create >0.",
    [E_INVALID_TARGET_QUBIT] = "Invalid target qubit. Note qubits are zero indexed.",
    [E_INVALID_CONTROL_QUBIT] = "Invalid control qubit. Note qubits are zero indexed.",
    [E_INVALID_STATE_INDEX] = "Invalid state index. Must be >=0 and <2^numQubits.",
    [E_INVALID_NUM_AMPS] = "Invalid number of amplitudes. Must be >=0 and <=2^numQubits.",
    [E_INVALID_OFFSET_NUM_AMPS] = "More amplitudes given than exist in the statevector from the given starting index.",
    [E_TARGET_IS_CONTROL] = "Control qubit cannot equal target qubit.",
    [E_TARGET_IN_CONTROLS] = "Control qubits cannot include target qubit.",
    [E_TARGETS_NOT_UNIQUE] = "The two target qubits must be unique.",
    [E_INVALID_NUM_CONTROLS] = "Invalid number of control qubits. Must be >0 and <numQubits.",
    [E_NON_UNITARY_MATRIX] = "Matrix is not unitary.",
    [E_NON_UNITARY_COMPLEX_PAIR] = "Compact matrix formed by given complex numbers is not unitary.",
    [E_ZERO_VECTOR] = "Invalid axis vector. Must be non-zero.",
    [E_SYS_TOO_BIG_TO_PRINT] = "Invalid system size. Cannot print output for systems greater than 5 qubits.",
    [E_COLLAPSE_STATE_ZERO_PROB] = "Can't collapse to state with zero probability.",
    [E_INVALID_QUBIT_OUTCOME] = "Invalid measurement outcome -- must be either 0 or 1.",
    [E_CANNOT_OPEN_FILE] = "Could not open file",
    [E_SECOND_ARG_MUST_BE_STATEVEC] = "Second argument must be a state-vector.",
    [E_MISMATCHING_QUREG_DIMENSIONS] = "Dimensions of the qubit registers don't match.",
    [E_MISMATCHING_QUREG_TYPES] = "Registers must both be state-vectors or both be density matrices.",
    [E_DEFINED_ONLY_FOR_STATEVECS] = "Operation valid only for state-vectors.",
    [E_DEFINED_ONLY_FOR_DENSMATRS] = "Operation valid only for density matrices.",
    [E_INVALID_PROB] = "Probabilities must be in [0, 1].",
    [E_UNNORM_PROBS] = "Probabilities must sum to ~1.",
    [E_INVALID_ONE_QUBIT_DEPHASE_PROB] = "The probability of a single qubit dephase error cannot exceed 1/2, which maximally mixes.",
    [E_INVALID_TWO_QUBIT_DEPHASE_PROB] = "The probability of a two-qubit qubit dephase error cannot exceed 3/4, which maximally mixes.",
    [E_INVALID_ONE_QUBIT_DEPOL_PROB] = "The probability of a single qubit depolarising error cannot exceed 3/4, which maximally mixes.",
    [E_INVALID_TWO_QUBIT_DEPOL_PROB] = "The probability of a two-qubit depolarising error cannot exceed 15/16, which maximally mixes."
};

void exitWithError(ErrorCode code, const char* func){
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, errorMessages[code]);
    printf("!!!\n");
    printf("exiting..\n");
    exit(code);
}

void QuESTAssert(int isValid, ErrorCode code, const char* func){
    if (!isValid) exitWithError(code, func);
}

int isComplexUnit(Complex alpha) {
    return (absReal(1 - sqrt(alpha.real*alpha.real + alpha.imag*alpha.imag)) < REAL_EPS); 
}

int isVectorUnit(qreal ux, qreal uy, qreal uz) {
    return (absReal(1 - sqrt(ux*ux + uy*uy + uz*uz)) < REAL_EPS );
}

int isComplexPairUnitary(Complex alpha, Complex beta) {
    return ( absReal( -1
                + alpha.real*alpha.real 
                + alpha.imag*alpha.imag
                + beta.real*beta.real 
                + beta.imag*beta.imag) < REAL_EPS );
}

int isMatrixUnitary(ComplexMatrix2 u) {
    if ( absReal( u.r0c0.real*u.r0c0.real 
                + u.r0c0.imag*u.r0c0.imag
                + u.r1c0.real*u.r1c0.real
                + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    if ( absReal( u.r0c1.real*u.r0c1.real 
                + u.r0c1.imag*u.r0c1.imag
                + u.r1c1.real*u.r1c1.real
                + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;
    if ( absReal( u.r0c0.real*u.r0c1.real 
                + u.r0c0.imag*u.r0c1.imag
                + u.r1c0.real*u.r1c1.real
                + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;
    if ( absReal( u.r0c1.real*u.r0c0.imag
                - u.r0c0.real*u.r0c1.imag
                + u.r1c1.real*u.r1c0.imag
                - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;
    return 1;
}

void validateCreateNumQubits(int numQubits, const char* caller) {
    QuESTAssert(numQubits>0, E_INVALID_NUM_QUBITS, caller);
}

void validateStateIndex(Qureg qureg, long long int stateInd, const char* caller) {
    long long int stateMax = 1LL << qureg.numQubitsRepresented;
    QuESTAssert(stateInd>=0 && stateInd<stateMax, E_INVALID_STATE_INDEX, caller);
}

void validateNumAmps(Qureg qureg, long long int startInd, long long int numAmps, const char* caller) {
    validateStateIndex(qureg, startInd, caller);
    QuESTAssert(numAmps >= 0 && numAmps <= qureg.numAmpsTotal, E_INVALID_NUM_AMPS, caller);
    QuESTAssert(numAmps + startInd <= qureg.numAmpsTotal, E_INVALID_OFFSET_NUM_AMPS, caller);
}

void validateTarget(Qureg qureg, int targetQubit, const char* caller) {
    QuESTAssert(targetQubit>=0 && targetQubit<qureg.numQubitsRepresented, E_INVALID_TARGET_QUBIT, caller);
}

void validateControl(Qureg qureg, int controlQubit, const char* caller) {
    QuESTAssert(controlQubit>=0 && controlQubit<qureg.numQubitsRepresented, E_INVALID_CONTROL_QUBIT, caller);
}

void validateControlTarget(Qureg qureg, int controlQubit, int targetQubit, const char* caller) {
    validateTarget(qureg, targetQubit, caller);
    validateControl(qureg, controlQubit, caller);
    QuESTAssert(controlQubit != targetQubit, E_TARGET_IS_CONTROL, caller);
}

void validateUniqueTargets(Qureg qureg, int qubit1, int qubit2, const char* caller) {
    validateTarget(qureg, qubit1, caller);
    validateTarget(qureg, qubit2, caller);
    QuESTAssert(qubit1 != qubit2, E_TARGETS_NOT_UNIQUE, caller);
}

void validateNumControls(Qureg qureg, const int numControlQubits, const char* caller) {
    // this could reject repeated qubits and cite "too many" but oh well
    QuESTAssert(numControlQubits>0 && numControlQubits<=qureg.numQubitsRepresented, E_INVALID_NUM_CONTROLS, caller);
}

void validateMultiControls(Qureg qureg, int* controlQubits, const int numControlQubits, const char* caller) {
    validateNumControls(qureg, numControlQubits, caller);
    for (int i=0; i < numControlQubits; i++) {
        validateControl(qureg, controlQubits[i], caller);
    }
}

void validateMultiControlsTarget(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, const char* caller) {
    validateTarget(qureg, targetQubit, caller);
    validateMultiControls(qureg, controlQubits, numControlQubits, caller);
    for (int i=0; i < numControlQubits; i++)
        QuESTAssert(controlQubits[i] != targetQubit, E_TARGET_IN_CONTROLS, caller);
}

void validateUnitaryMatrix(ComplexMatrix2 u, const char* caller) {
    QuESTAssert(isMatrixUnitary(u), E_NON_UNITARY_MATRIX, caller);
}

void validateUnitaryComplexPair(Complex alpha, Complex beta, const char* caller) {
    QuESTAssert(isComplexPairUnitary(alpha, beta), E_NON_UNITARY_COMPLEX_PAIR, caller);
}

void validateVector(Vector vec, const char* caller) {
    QuESTAssert(getVectorMagnitude(vec) > REAL_EPS, E_ZERO_VECTOR, caller);
}

void validateStateVecQureg(Qureg qureg, const char* caller) {
    QuESTAssert( ! qureg.isDensityMatrix, E_DEFINED_ONLY_FOR_STATEVECS, caller);
}

void validateDensityMatrQureg(Qureg qureg, const char* caller) {
    QuESTAssert(qureg.isDensityMatrix, E_DEFINED_ONLY_FOR_DENSMATRS, caller);
}

void validateOutcome(int outcome, const char* caller) {
    QuESTAssert(outcome==0 || outcome==1, E_INVALID_QUBIT_OUTCOME, caller);
}

void validateMeasurementProb(qreal prob, const char* caller) {
    QuESTAssert(prob>REAL_EPS, E_COLLAPSE_STATE_ZERO_PROB, caller);
}

void validateMatchingQuregDims(Qureg qureg1, Qureg qureg2, const char *caller) {
    QuESTAssert(qureg1.numQubitsRepresented==qureg2.numQubitsRepresented, E_MISMATCHING_QUREG_DIMENSIONS, caller);
}

void validateMatchingQuregTypes(Qureg qureg1, Qureg qureg2, const char *caller) {
    QuESTAssert(qureg1.isDensityMatrix==qureg2.isDensityMatrix, E_MISMATCHING_QUREG_TYPES, caller);
}

void validateSecondQuregStateVec(Qureg qureg2, const char *caller) {
    QuESTAssert( ! qureg2.isDensityMatrix, E_SECOND_ARG_MUST_BE_STATEVEC, caller);
}

void validateFileOpened(int found, const char* caller) {
    QuESTAssert(found, E_CANNOT_OPEN_FILE, caller);
}

void validateProb(qreal prob, const char* caller) {
    QuESTAssert(prob >= 0 && prob <= 1, E_INVALID_PROB, caller);
}

void validateNormProbs(qreal prob1, qreal prob2, const char* caller) {
    validateProb(prob1, caller);
    validateProb(prob2, caller);
    
    qreal sum = prob1 + prob2;
    QuESTAssert(absReal(1 - sum) < REAL_EPS, E_UNNORM_PROBS, caller);
}

void validateOneQubitDephaseProb(qreal prob, const char* caller) {
    validateProb(prob, caller);
    QuESTAssert(prob <= 1/2.0, E_INVALID_ONE_QUBIT_DEPHASE_PROB, caller);
}

void validateTwoQubitDephaseProb(qreal prob, const char* caller) {
    validateProb(prob, caller);
    QuESTAssert(prob <= 3/4.0, E_INVALID_TWO_QUBIT_DEPHASE_PROB, caller);
}

void validateOneQubitDepolProb(qreal prob, const char* caller) {
    validateProb(prob, caller);
    QuESTAssert(prob <= 3/4.0, E_INVALID_ONE_QUBIT_DEPOL_PROB, caller);
}

void validateOneQubitDampingProb(qreal prob, const char* caller) {
    validateProb(prob, caller);
    QuESTAssert(prob <= 1.0, E_INVALID_ONE_QUBIT_DEPOL_PROB, caller);
}

void validateTwoQubitDepolProb(qreal prob, const char* caller) {
    validateProb(prob, caller);
    QuESTAssert(prob <= 15/16.0, E_INVALID_TWO_QUBIT_DEPOL_PROB, caller);
}




#ifdef __cplusplus
}
#endif
