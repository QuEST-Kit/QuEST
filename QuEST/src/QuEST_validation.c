// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Provides validation of user input
 *
 * @author Tyson Jones
 * @author Ania Brown (exitWithError(), QuESTAssert(), original testing of:
 *        qubit indices, unitarity, valid collapse probability)
 * @author Balint Koczor (Kraus maps)
 * @author Nicolas Vogt of HQS (one-qubit damping)
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
# include <stdint.h>

/* buffer for an error message which contains formatters. This must be global, 
 * since if a function writes to a local buffer then throws an error, the 
 * local buffer may be cleared and dangle before the error catcher can process it.
 */
char errMsgBuffer[1024];

typedef enum {
    E_SUCCESS=0,
    E_INVALID_NUM_RANKS,
    E_INVALID_NUM_CREATE_QUBITS,
    E_INVALID_QUBIT_INDEX,
    E_INVALID_TARGET_QUBIT,
    E_INVALID_CONTROL_QUBIT,
    E_INVALID_STATE_INDEX,
    E_INVALID_AMP_INDEX,
    E_INVALID_ELEM_INDEX,
    E_INVALID_NUM_AMPS,
    E_INVALID_NUM_ELEMS,
    E_INVALID_OFFSET_NUM_AMPS_QUREG,
    E_INVALID_OFFSET_NUM_ELEMS_DIAG,
    E_TARGET_IS_CONTROL,
    E_TARGET_IN_CONTROLS,
    E_CONTROL_TARGET_COLLISION,
    E_QUBITS_NOT_UNIQUE,
    E_TARGETS_NOT_UNIQUE,
    E_CONTROLS_NOT_UNIQUE,
    E_INVALID_NUM_QUBITS,
    E_INVALID_NUM_TARGETS,
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
    E_INVALID_TWO_QUBIT_DEPOL_PROB,
    E_INVALID_ONE_QUBIT_PAULI_PROBS,
    E_INVALID_CONTROLS_BIT_STATE,
    E_INVALID_PAULI_CODE,
    E_INVALID_NUM_SUM_TERMS,
    E_CANNOT_FIT_MULTI_QUBIT_MATRIX,
    E_INVALID_UNITARY_SIZE,
    E_COMPLEX_MATRIX_NOT_INIT,
    E_INVALID_NUM_ONE_QUBIT_KRAUS_OPS,
    E_INVALID_NUM_TWO_QUBIT_KRAUS_OPS,
    E_INVALID_NUM_N_QUBIT_KRAUS_OPS,
    E_INVALID_KRAUS_OPS,
    E_MISMATCHING_NUM_TARGS_KRAUS_SIZE,
    E_DISTRIB_QUREG_TOO_SMALL,
    E_DISTRIB_DIAG_OP_TOO_SMALL,
    E_NUM_AMPS_EXCEED_TYPE,
    E_INVALID_PAULI_HAMIL_PARAMS,
    E_INVALID_PAULI_HAMIL_FILE_PARAMS,
    E_CANNOT_PARSE_PAULI_HAMIL_FILE_COEFF,
    E_CANNOT_PARSE_PAULI_HAMIL_FILE_PAULI,
    E_INVALID_PAULI_HAMIL_FILE_PAULI_CODE,
    E_MISMATCHING_PAULI_HAMIL_QUREG_NUM_QUBITS,
    E_INVALID_TROTTER_ORDER,
    E_INVALID_TROTTER_REPS,
    E_MISMATCHING_QUREG_DIAGONAL_OP_SIZE,
    E_DIAGONAL_OP_NOT_INITIALISED
} ErrorCode;

static const char* errorMessages[] = {
    [E_INVALID_NUM_RANKS] = "Invalid number of nodes. Distributed simulation can only make use of a power-of-2 number of node.",
    [E_INVALID_NUM_CREATE_QUBITS] = "Invalid number of qubits. Must create >0.",
    [E_INVALID_QUBIT_INDEX] = "Invalid qubit index. Must be >=0 and <numQubits.",
    [E_INVALID_TARGET_QUBIT] = "Invalid target qubit. Must be >=0 and <numQubits.",
    [E_INVALID_CONTROL_QUBIT] = "Invalid control qubit. Must be >=0 and <numQubits.",
    [E_INVALID_STATE_INDEX] = "Invalid state index. Must be >=0 and <2^numQubits.",
    [E_INVALID_AMP_INDEX] = "Invalid amplitude index. Must be >=0 and <2^numQubits.",
    [E_INVALID_ELEM_INDEX] = "Invalid element index. Must be >=0 and <2^numQubits.",
    [E_INVALID_NUM_AMPS] = "Invalid number of amplitudes. Must be >=0 and <=2^numQubits.",
    [E_INVALID_NUM_ELEMS] = "Invalid number of elements. Must be >=0 and <=2^numQubits.",
    [E_INVALID_OFFSET_NUM_AMPS_QUREG] = "More amplitudes given than exist in the statevector from the given starting index.",
    [E_INVALID_OFFSET_NUM_ELEMS_DIAG] = "More elements given than exist in the diagonal operator from the given starting index.",
    [E_TARGET_IS_CONTROL] = "Control qubit cannot equal target qubit.",
    [E_TARGET_IN_CONTROLS] = "Control qubits cannot include target qubit.",
    [E_CONTROL_TARGET_COLLISION] = "Control and target qubits must be disjoint.",
    [E_QUBITS_NOT_UNIQUE] = "The qubits must be unique.",
    [E_TARGETS_NOT_UNIQUE] = "The target qubits must be unique.",
    [E_CONTROLS_NOT_UNIQUE] = "The control qubits should be unique.",
    [E_INVALID_NUM_QUBITS] = "Invalid number of qubits. Must be >0 and <=numQubits.",
    [E_INVALID_NUM_TARGETS] = "Invalid number of target qubits. Must be >0 and <=numQubits.",
    [E_INVALID_NUM_CONTROLS] = "Invalid number of control qubits. Must be >0 and <numQubits.",
    [E_NON_UNITARY_MATRIX] = "Matrix is not unitary.",
    [E_NON_UNITARY_COMPLEX_PAIR] = "Compact matrix formed by given complex numbers is not unitary.",
    [E_ZERO_VECTOR] = "Invalid axis vector. Must be non-zero.",
    [E_SYS_TOO_BIG_TO_PRINT] = "Invalid system size. Cannot print output for systems greater than 5 qubits.",
    [E_COLLAPSE_STATE_ZERO_PROB] = "Can't collapse to state with zero probability.",
    [E_INVALID_QUBIT_OUTCOME] = "Invalid measurement outcome -- must be either 0 or 1.",
    [E_CANNOT_OPEN_FILE] = "Could not open file (%s).",
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
    [E_INVALID_TWO_QUBIT_DEPOL_PROB] = "The probability of a two-qubit depolarising error cannot exceed 15/16, which maximally mixes.",
    [E_INVALID_ONE_QUBIT_PAULI_PROBS] = "The probability of any X, Y or Z error cannot exceed the probability of no error.",
    [E_INVALID_CONTROLS_BIT_STATE] = "The state of the control qubits must be a bit sequence (0s and 1s).",
    [E_INVALID_PAULI_CODE] = "Invalid Pauli code. Codes must be 0 (or PAULI_I), 1 (PAULI_X), 2 (PAULI_Y) or 3 (PAULI_Z) to indicate the identity, X, Y and Z operators respectively.",
    [E_INVALID_NUM_SUM_TERMS] = "Invalid number of terms in the Pauli sum. The number of terms must be >0.",
    [E_CANNOT_FIT_MULTI_QUBIT_MATRIX] = "The specified matrix targets too many qubits; the batches of amplitudes to modify cannot all fit in a single distributed node's memory allocation.",
    [E_INVALID_UNITARY_SIZE] = "The matrix size does not match the number of target qubits.",
    [E_COMPLEX_MATRIX_NOT_INIT] = "The ComplexMatrixN was not successfully created (possibly insufficient memory available).",
    [E_INVALID_NUM_ONE_QUBIT_KRAUS_OPS] = "At least 1 and at most 4 single qubit Kraus operators may be specified.",
    [E_INVALID_NUM_TWO_QUBIT_KRAUS_OPS] = "At least 1 and at most 16 two-qubit Kraus operators may be specified.",
    [E_INVALID_NUM_N_QUBIT_KRAUS_OPS] = "At least 1 and at most 4*N^2 of N-qubit Kraus operators may be specified.",
    [E_INVALID_KRAUS_OPS] = "The specified Kraus map is not a completely positive, trace preserving map.",
    [E_MISMATCHING_NUM_TARGS_KRAUS_SIZE] = "Every Kraus operator must be of the same number of qubits as the number of targets.",
    [E_DISTRIB_QUREG_TOO_SMALL] = "Too few qubits. The created qureg must have at least one amplitude per node used in distributed simulation.",
    [E_DISTRIB_DIAG_OP_TOO_SMALL] = "Too few qubits. The created DiagonalOp must contain at least one element per node used in distributed simulation.",
    [E_NUM_AMPS_EXCEED_TYPE] = "Too many qubits (max of log2(SIZE_MAX)). Cannot store the number of amplitudes per-node in the size_t type.",
    [E_INVALID_PAULI_HAMIL_PARAMS] = "The number of qubits and terms in the PauliHamil must be strictly positive.",
    [E_INVALID_PAULI_HAMIL_FILE_PARAMS] = "The number of qubits and terms in the PauliHamil file (%s) must be strictly positive.",
    [E_CANNOT_PARSE_PAULI_HAMIL_FILE_COEFF] = "Failed to parse the next expected term coefficient in PauliHamil file (%s).",
    [E_CANNOT_PARSE_PAULI_HAMIL_FILE_PAULI] = "Failed to parse the next expected Pauli code in PauliHamil file (%s).",
    [E_INVALID_PAULI_HAMIL_FILE_PAULI_CODE] = "The PauliHamil file (%s) contained an invalid pauli code (%d). Codes must be 0 (or PAULI_I), 1 (PAULI_X), 2 (PAULI_Y) or 3 (PAULI_Z) to indicate the identity, X, Y and Z operators respectively.",
    [E_MISMATCHING_PAULI_HAMIL_QUREG_NUM_QUBITS] = "The PauliHamil must act on the same number of qubits as exist in the Qureg.",
    [E_INVALID_TROTTER_ORDER] = "The Trotterisation order must be 1, or an even number (for higher-order Suzuki symmetrized expansions).",
    [E_INVALID_TROTTER_REPS] = "The number of Trotter repetitions must be >=1.",
    [E_MISMATCHING_QUREG_DIAGONAL_OP_SIZE] = "The qureg must represent an equal number of qubits as that in the applied diagonal operator.",
    [E_DIAGONAL_OP_NOT_INITIALISED] = "The diagonal operator has not been initialised through createDiagonalOperator()."
};

void exitWithError(const char* msg, const char* func) {
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, msg);
    printf("!!!\n");
    printf("exiting..\n");
    exit(1);
}

#pragma weak invalidQuESTInputError
void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
    exitWithError(errMsg, errFunc);
}

void QuESTAssert(int isValid, ErrorCode code, const char* func){
    if (!isValid) invalidQuESTInputError(errorMessages[code], func);
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

#define macro_isMatrixUnitary(m, dim, retVal) { \
    /* elemRe_ and elemIm_ must not exist in caller scope */ \
    qreal elemRe_, elemIm_; \
    retVal = 1; \
    /* check m * ConjugateTranspose(m) == Identity */ \
    for (int r=0; r < (dim); r++) { \
        for (int c=0; c < (dim); c++) { \
            /* m[r][...] * ConjugateTranspose(m)[...][c] */ \
            elemRe_ = 0; \
            elemIm_ = 0; \
            for (int i=0; i < (dim); i++) { \
                /* m[r][i] * conj(m[c][i]) */ \
                elemRe_ += m.real[r][i]*m.real[c][i] + m.imag[r][i]*m.imag[c][i]; \
                elemIm_ += m.imag[r][i]*m.real[c][i] - m.real[r][i]*m.imag[c][i]; \
            } \
            /* check distance from identity */ \
            if ((absReal(elemIm_) > REAL_EPS) || \
                (r == c && absReal(elemRe_ - 1) > REAL_EPS) || \
                (r != c && absReal(elemRe_    ) > REAL_EPS)) { \
                retVal = 0; \
                break; \
            } \
        } \
        if (retVal == 0) \
            break; \
    } \
}
int isMatrix2Unitary(ComplexMatrix2 u) {
    int dim = 2;
    int retVal;
    macro_isMatrixUnitary(u, dim, retVal);
    return retVal;
}
int isMatrix4Unitary(ComplexMatrix4 u) {
    int dim = 4;
    int retVal;
    macro_isMatrixUnitary(u, dim, retVal);
    return retVal;
}
int isMatrixNUnitary(ComplexMatrixN u) {
    int dim = 1 << u.numQubits;
    int retVal;
    macro_isMatrixUnitary(u, dim, retVal);
    return retVal;
}

#define macro_isCompletelyPositiveMap(ops, numOps, opDim) { \
    for (int r=0; r<(opDim); r++) { \
        for (int c=0; c<(opDim); c++) { \
            qreal elemRe_ = 0; \
            qreal elemIm_ = 0; \
            for (int n=0; n<(numOps); n++) { \
                for (int k=0; k<(opDim); k++) { \
                    elemRe_ += ops[n].real[k][r]*ops[n].real[k][c] + ops[n].imag[k][r]*ops[n].imag[k][c]; \
                    elemIm_ += ops[n].real[k][r]*ops[n].imag[k][c] - ops[n].imag[k][r]*ops[n].real[k][c]; \
                } \
            } \
            qreal dist_ = absReal(elemIm_) + absReal(elemRe_ - ((r==c)? 1:0)); \
            if (dist_ > REAL_EPS) \
                return 0; \
        } \
    } \
    return 1; \
}
int isCompletelyPositiveMap2(ComplexMatrix2 *ops, int numOps) {
    macro_isCompletelyPositiveMap(ops, numOps, 2);    
}
int isCompletelyPositiveMap4(ComplexMatrix4 *ops, int numOps) {
    macro_isCompletelyPositiveMap(ops, numOps, 4);
}
int isCompletelyPositiveMapN(ComplexMatrixN *ops, int numOps) {
    int opDim = 1 << ops[0].numQubits;
    macro_isCompletelyPositiveMap(ops, numOps, opDim);
}

int isValidPauliCode(enum pauliOpType code) {
    return (code==PAULI_I || code==PAULI_X || code==PAULI_Y || code==PAULI_Z);
}

int areUniqueQubits(int* qubits, int numQubits) {
    long long int mask = 0;
    long long int bit;
    for (int q=0; q < numQubits; q++) {
        bit = 1LL << qubits[q];
        if (mask & bit)
            return 0;
        mask |= bit;
    }
    return 1;
}

/** returns log2 of numbers which must be gauranteed to be 2^n */
unsigned int calcLog2(long unsigned int num) {
    unsigned int l = 0;
    while (num >>= 1)
        l++;
    return l;
}

void validateNumRanks(int numRanks, const char* caller) {

    /* silly but robust way to determine if numRanks is a power of 2, 
     * in lieu of bit-twiddling (e.g. AND with all-ones) which may be 
     * system / precsision dependent 
     */
    int isValid = 0;
    for (int exp2 = 1; exp2 <= numRanks; exp2 *= 2)
        if (exp2 == numRanks)
            isValid = 1;
    
    QuESTAssert(isValid, E_INVALID_NUM_RANKS, caller);
}

void validateNumQubitsInQureg(int numQubits, int numRanks, const char* caller) {
    QuESTAssert(numQubits>0, E_INVALID_NUM_CREATE_QUBITS, caller);
    
    // mustn't be more amplitudes than can fit in the type
    unsigned int maxQubits = calcLog2(SIZE_MAX);
    QuESTAssert( numQubits <= maxQubits, E_NUM_AMPS_EXCEED_TYPE, caller);
    
    // must be at least one amplitude per node
    long unsigned int numAmps = (1UL<<numQubits);
    QuESTAssert(numAmps >= numRanks, E_DISTRIB_QUREG_TOO_SMALL, caller);
}
 
void validateNumQubitsInMatrix(int numQubits, const char* caller) {
    QuESTAssert(numQubits>0, E_INVALID_NUM_QUBITS, caller);
}

void validateNumQubitsInDiagOp(int numQubits, int numRanks, const char* caller) {
    QuESTAssert(numQubits>0, E_INVALID_NUM_CREATE_QUBITS, caller);
    
    // mustn't be more amplitudes than can fit in the type
    unsigned int maxQubits = calcLog2(SIZE_MAX);
    QuESTAssert( numQubits <= maxQubits, E_NUM_AMPS_EXCEED_TYPE, caller);
    
    // must be at least one amplitude per node
    long unsigned int numAmps = (1UL<<numQubits);
    QuESTAssert(numAmps >= numRanks, E_DISTRIB_DIAG_OP_TOO_SMALL, caller);
}

void validateStateIndex(Qureg qureg, long long int stateInd, const char* caller) {
    long long int stateMax = 1LL << qureg.numQubitsRepresented;
    QuESTAssert(stateInd>=0 && stateInd<stateMax, E_INVALID_STATE_INDEX, caller);
}

void validateAmpIndex(Qureg qureg, long long int ampInd, const char* caller) {
    long long int indMax = 1LL << qureg.numQubitsRepresented;
    QuESTAssert(ampInd>=0 && ampInd<indMax, E_INVALID_AMP_INDEX, caller);
}

void validateNumAmps(Qureg qureg, long long int startInd, long long int numAmps, const char* caller) {
    validateAmpIndex(qureg, startInd, caller);
    QuESTAssert(numAmps >= 0 && numAmps <= qureg.numAmpsTotal, E_INVALID_NUM_AMPS, caller);
    QuESTAssert(numAmps + startInd <= qureg.numAmpsTotal, E_INVALID_OFFSET_NUM_AMPS_QUREG, caller);
}

void validateNumElems(DiagonalOp op, long long int startInd, long long int numElems, const char* caller) {
    long long int indMax = 1LL << op.numQubits;
    QuESTAssert(startInd >= 0 && startInd < indMax, E_INVALID_ELEM_INDEX, caller);
    QuESTAssert(numElems >= 0 && numElems <= indMax, E_INVALID_NUM_ELEMS, caller);
    QuESTAssert(numElems + startInd <= indMax, E_INVALID_OFFSET_NUM_ELEMS_DIAG, caller);
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

void validateNumTargets(Qureg qureg, int numTargetQubits, const char* caller) {
    QuESTAssert(numTargetQubits>0 && numTargetQubits<=qureg.numQubitsRepresented, E_INVALID_NUM_TARGETS, caller);
}

void validateNumControls(Qureg qureg, int numControlQubits, const char* caller) {
    QuESTAssert(numControlQubits>0 && numControlQubits<qureg.numQubitsRepresented, E_INVALID_NUM_CONTROLS, caller);
}

void validateMultiTargets(Qureg qureg, int* targetQubits, int numTargetQubits, const char* caller) {
    validateNumTargets(qureg, numTargetQubits, caller);
    for (int i=0; i < numTargetQubits; i++) 
        validateTarget(qureg, targetQubits[i], caller);
        
    QuESTAssert(areUniqueQubits(targetQubits, numTargetQubits), E_TARGETS_NOT_UNIQUE, caller);
}

void validateMultiControls(Qureg qureg, int* controlQubits, int numControlQubits, const char* caller) {
    validateNumControls(qureg, numControlQubits, caller);
    for (int i=0; i < numControlQubits; i++)
        validateControl(qureg, controlQubits[i], caller);
        
    QuESTAssert(areUniqueQubits(controlQubits, numControlQubits), E_CONTROLS_NOT_UNIQUE, caller);
}

void validateMultiQubits(Qureg qureg, int* qubits, int numQubits, const char* caller) {
    QuESTAssert(numQubits>0 && numQubits<=qureg.numQubitsRepresented, E_INVALID_NUM_QUBITS, caller);
    for (int i=0; i < numQubits; i++)
        QuESTAssert(qubits[i]>=0 && qubits[i]<qureg.numQubitsRepresented, E_INVALID_QUBIT_INDEX, caller);
        
    QuESTAssert(areUniqueQubits(qubits, numQubits), E_QUBITS_NOT_UNIQUE, caller);
}

void validateMultiControlsTarget(Qureg qureg, int* controlQubits, int numControlQubits, int targetQubit, const char* caller) {
    validateTarget(qureg, targetQubit, caller);
    validateMultiControls(qureg, controlQubits, numControlQubits, caller);
    for (int i=0; i < numControlQubits; i++)
        QuESTAssert(controlQubits[i] != targetQubit, E_TARGET_IN_CONTROLS, caller);
}

void validateMultiControlsMultiTargets(Qureg qureg, int* controlQubits, int numControlQubits, int* targetQubits, int numTargetQubits, const char* caller) {
    validateMultiControls(qureg, controlQubits, numControlQubits, caller);
    validateMultiTargets(qureg, targetQubits, numTargetQubits, caller);
    long long int ctrlMask = getQubitBitMask(controlQubits, numControlQubits);
    long long int targMask = getQubitBitMask(targetQubits, numTargetQubits);
    int overlap = ctrlMask & targMask;
    QuESTAssert(!overlap, E_CONTROL_TARGET_COLLISION, caller);
}

void validateControlState(int* controlState, int numControlQubits, const char* caller) {
    for (int i=0; i < numControlQubits; i++)
        QuESTAssert(controlState[i] == 0 || controlState[i] == 1, E_INVALID_CONTROLS_BIT_STATE, caller);
}

void validateMultiQubitMatrixFitsInNode(Qureg qureg, int numTargets, const char* caller) {
    QuESTAssert(qureg.numAmpsPerChunk >= (1LL << numTargets), E_CANNOT_FIT_MULTI_QUBIT_MATRIX, caller);
}

void validateOneQubitUnitaryMatrix(ComplexMatrix2 u, const char* caller) {
    QuESTAssert(isMatrix2Unitary(u), E_NON_UNITARY_MATRIX, caller);
}

void validateTwoQubitUnitaryMatrix(Qureg qureg, ComplexMatrix4 u, const char* caller) {
    validateMultiQubitMatrixFitsInNode(qureg, 2, caller);
    QuESTAssert(isMatrix4Unitary(u), E_NON_UNITARY_MATRIX, caller);
}

void validateMatrixInit(ComplexMatrixN matr, const char* caller) {
    
    /* note that for (most) compilers which don't automatically initialise 
     * pointers to NULL, this can only be used to check the mallocs in createComplexMatrixN
     * succeeded. It can not be used to differentiate whether a user actually attempted 
     * to initialise or create their ComplexMatrixN instance.
     */
    QuESTAssert(matr.real != NULL && matr.imag != NULL, E_COMPLEX_MATRIX_NOT_INIT, caller);
}

void validateMultiQubitMatrix(Qureg qureg, ComplexMatrixN u, int numTargs, const char* caller) { 
    validateMatrixInit(u, caller);
    validateMultiQubitMatrixFitsInNode(qureg, numTargs, caller);
    QuESTAssert(numTargs == u.numQubits, E_INVALID_UNITARY_SIZE, caller);
}

void validateMultiQubitUnitaryMatrix(Qureg qureg, ComplexMatrixN u, int numTargs, const char* caller) { 
    validateMultiQubitMatrix(qureg, u, numTargs, caller);
    QuESTAssert(isMatrixNUnitary(u), E_NON_UNITARY_MATRIX, caller);
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

void validateFileOpened(int opened, char* fn, const char* caller) {
    if (!opened) {
        
        sprintf(errMsgBuffer, errorMessages[E_CANNOT_OPEN_FILE], fn);
        invalidQuESTInputError(errMsgBuffer, caller);
    }
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

void validateOneQubitPauliProbs(qreal probX, qreal probY, qreal probZ, const char* caller) {
    validateProb(probX, caller);
    validateProb(probY, caller);
    validateProb(probZ, caller);
    
    qreal probNoError = 1 - probX - probY - probZ;
    QuESTAssert(probX <= probNoError, E_INVALID_ONE_QUBIT_PAULI_PROBS, caller);
    QuESTAssert(probY <= probNoError, E_INVALID_ONE_QUBIT_PAULI_PROBS, caller);
    QuESTAssert(probZ <= probNoError, E_INVALID_ONE_QUBIT_PAULI_PROBS, caller);
}

void validatePauliCodes(enum pauliOpType* pauliCodes, int numPauliCodes, const char* caller) {
    for (int i=0; i < numPauliCodes; i++) {
        enum pauliOpType code = pauliCodes[i];
        QuESTAssert(isValidPauliCode(code), E_INVALID_PAULI_CODE, caller);
    }
}

void validateNumPauliSumTerms(int numTerms, const char* caller) {
    QuESTAssert(numTerms > 0, E_INVALID_NUM_SUM_TERMS, caller);
}

void validateOneQubitKrausMap(Qureg qureg, ComplexMatrix2* ops, int numOps, const char* caller) {
    int opNumQubits = 1;
    int superOpNumQubits = 2*opNumQubits;
    int maxNumOps = superOpNumQubits*superOpNumQubits;
    QuESTAssert(numOps > 0 && numOps <= maxNumOps, E_INVALID_NUM_ONE_QUBIT_KRAUS_OPS, caller);
    
    validateMultiQubitMatrixFitsInNode(qureg, superOpNumQubits, caller);
    
    int isPos = isCompletelyPositiveMap2(ops, numOps);
    QuESTAssert(isPos, E_INVALID_KRAUS_OPS, caller);
}

void validateTwoQubitKrausMap(Qureg qureg, ComplexMatrix4* ops, int numOps, const char* caller) {
    int opNumQubits = 2;
    int superOpNumQubits = 2*opNumQubits;
    int maxNumOps = superOpNumQubits*superOpNumQubits;
    QuESTAssert(numOps > 0 && numOps <= maxNumOps, E_INVALID_NUM_TWO_QUBIT_KRAUS_OPS, caller);
    
    validateMultiQubitMatrixFitsInNode(qureg, superOpNumQubits, caller);

    int isPos = isCompletelyPositiveMap4(ops, numOps);
    QuESTAssert(isPos, E_INVALID_KRAUS_OPS, caller);
}

void validateMultiQubitKrausMap(Qureg qureg, int numTargs, ComplexMatrixN* ops, int numOps, const char* caller) {
    int opNumQubits = numTargs;
    int superOpNumQubits = 2*opNumQubits;
    int maxNumOps = superOpNumQubits*superOpNumQubits;
    QuESTAssert(numOps>0 && numOps <= maxNumOps, E_INVALID_NUM_N_QUBIT_KRAUS_OPS, caller);
        
    for (int n=0; n<numOps; n++) {
        validateMatrixInit(ops[n], __func__);
        QuESTAssert(ops[n].numQubits == numTargs, E_MISMATCHING_NUM_TARGS_KRAUS_SIZE, caller);    
    }
    
    validateMultiQubitMatrixFitsInNode(qureg, superOpNumQubits, caller);
    
    int isPos = isCompletelyPositiveMapN(ops, numOps);
    QuESTAssert(isPos, E_INVALID_KRAUS_OPS, caller);
}

void validateHamilParams(int numQubits, int numTerms, const char* caller) {
    QuESTAssert(numQubits > 0 && numTerms > 0, E_INVALID_PAULI_HAMIL_PARAMS, caller);
}

void validatePauliHamil(PauliHamil hamil, const char* caller) {
    validateHamilParams(hamil.numQubits, hamil.numSumTerms, caller);
    validatePauliCodes(hamil.pauliCodes, hamil.numSumTerms*hamil.numQubits, caller);
}

void validateMatchingQuregPauliHamilDims(Qureg qureg, PauliHamil hamil, const char* caller) {
    QuESTAssert(hamil.numQubits == qureg.numQubitsRepresented, E_MISMATCHING_PAULI_HAMIL_QUREG_NUM_QUBITS, caller);
}

void validateHamilFileParams(int numQubits, int numTerms, FILE* file, char* fn, const char* caller) {
    if (!(numQubits > 0 && numTerms > 0)) {
        fclose(file);

        sprintf(errMsgBuffer, errorMessages[E_INVALID_PAULI_HAMIL_FILE_PARAMS], fn);
        invalidQuESTInputError(errMsgBuffer, caller);
    }
}

void validateHamilFileCoeffParsed(int parsed, PauliHamil h, FILE* file, char* fn, const char* caller) {
    if (!parsed) {
        destroyPauliHamil(h);
        fclose(file);
        
        sprintf(errMsgBuffer, errorMessages[E_CANNOT_PARSE_PAULI_HAMIL_FILE_COEFF], fn);
        invalidQuESTInputError(errMsgBuffer, caller);
    }
}

void validateHamilFilePauliParsed(int parsed, PauliHamil h, FILE* file, char* fn, const char* caller) {
    if (!parsed) {
        destroyPauliHamil(h);
        fclose(file);
        
        sprintf(errMsgBuffer, errorMessages[E_CANNOT_PARSE_PAULI_HAMIL_FILE_PAULI], fn);
        invalidQuESTInputError(errMsgBuffer, caller);
    }
}

void validateHamilFilePauliCode(enum pauliOpType code, PauliHamil h, FILE* file, char* fn, const char* caller) {
    if (!isValidPauliCode(code)) {
        destroyPauliHamil(h);
        fclose(file);
        
        sprintf(errMsgBuffer, errorMessages[E_INVALID_PAULI_HAMIL_FILE_PAULI_CODE], fn, code);
        invalidQuESTInputError(errMsgBuffer, caller);
    }
}

void validateTrotterParams(int order, int reps, const char* caller) {
    int isEven = (order % 2) == 0;
    QuESTAssert(order > 0 && (isEven || order==1), E_INVALID_TROTTER_ORDER, caller);
    QuESTAssert(reps > 0, E_INVALID_TROTTER_REPS, caller);
}

void validateDiagOpInit(DiagonalOp op, const char* caller) {
    QuESTAssert(op.real != NULL && op.imag != NULL, E_DIAGONAL_OP_NOT_INITIALISED, caller);
}

void validateDiagonalOp(Qureg qureg, DiagonalOp op, const char* caller) {
    validateDiagOpInit(op, caller);
    QuESTAssert(qureg.numQubitsRepresented == op.numQubits, E_MISMATCHING_QUREG_DIAGONAL_OP_SIZE, caller);
}

#ifdef __cplusplus
}
#endif
