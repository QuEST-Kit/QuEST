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
    E_CONTROLS_NOT_UNIQUE,
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
    E_INVALID_KRAUS_OPS
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
    [E_TARGETS_NOT_UNIQUE] = "The target qubits must be unique.",
    [E_CONTROLS_NOT_UNIQUE] = "The control qubits should be unique.",
    [E_INVALID_NUM_TARGETS] = "Invalid number of target qubits. Must be >0 and <=numQubits.",
    [E_INVALID_NUM_CONTROLS] = "Invalid number of control qubits. Must be >0 and <numQubits.",
    [E_NON_UNITARY_MATRIX] = "Matrix is not unitary.",
    [E_NON_UNITARY_COMPLEX_PAIR] = "Compact matrix formed by given complex numbers is not unitary.",
    [E_ZERO_VECTOR] = "Invalid axis vector. Must be non-zero.",
    [E_SYS_TOO_BIG_TO_PRINT] = "Invalid system size. Cannot print output for systems greater than 5 qubits.",
    [E_COLLAPSE_STATE_ZERO_PROB] = "Can't collapse to state with zero probability.",
    [E_INVALID_QUBIT_OUTCOME] = "Invalid measurement outcome -- must be either 0 or 1.",
    [E_CANNOT_OPEN_FILE] = "Could not open file.",
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
    [E_INVALID_PAULI_CODE] = "Invalid Pauli code. Codes must be 0 (or PAULI_I), 1 (PAULI_X), 2 (PAULI_Y) or 3 (PAULI_Z) to indicate the identity, X, Y and Z gates respectively.",
    [E_INVALID_NUM_SUM_TERMS] = "Invalid number of terms in the Pauli sum. The number of terms must be >0.",
    [E_CANNOT_FIT_MULTI_QUBIT_MATRIX] = "The specified matrix targets too many qubits; the batches of amplitudes to modify cannot all fit in a single distributed node's memory allocation.",
    [E_INVALID_UNITARY_SIZE] = "The matrix size does not match the number of target qubits.",
    [E_COMPLEX_MATRIX_NOT_INIT] = "The ComplexMatrixN wasn't initialised with createComplexMatrix().",
    [E_INVALID_NUM_ONE_QUBIT_KRAUS_OPS] = "At least 1 and at most 4 single qubit Kraus operators may be specified.",
    [E_INVALID_NUM_TWO_QUBIT_KRAUS_OPS] = "At least 1 and at most 16 two-qubit Kraus operators may be specified.",
    [E_INVALID_KRAUS_OPS] = "The specified Kraus map is not a completely positive, trace preserving map."
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

int isMatrix2Unitary(ComplexMatrix2 u) {
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

Complex getMatrixProductElement(Complex* row, Complex* col, int dim) {
    Complex elem = {.real = 0, .imag = 0};
    for (int i=0; i < dim; i++) {
        elem.real += row[i].real*col[i].real - row[i].imag*col[i].imag;
        elem.imag += row[i].imag*col[i].real + row[i].real*col[i].imag;
    }
    return elem;
}

ComplexMatrix2 getMatrix2Product(ComplexMatrix2 a, ComplexMatrix2 b) {
    ComplexMatrix2 prod;
    Complex r0[2] = {a.r0c0,a.r0c1};
    Complex r1[2] = {a.r1c0,a.r1c1};
    Complex c0[2] = {b.r0c0,b.r1c0};
    Complex c1[2] = {b.r0c1,b.r1c1};
    
    prod.r0c0 = getMatrixProductElement(r0,c0,2);
    prod.r0c1 = getMatrixProductElement(r0,c1,2);
    prod.r1c0 = getMatrixProductElement(r1,c0,2);
    prod.r1c1 = getMatrixProductElement(r1,c1,2);
    
    return prod;
}

ComplexMatrix4 getMatrix4Product(ComplexMatrix4 a, ComplexMatrix4 b) {
    ComplexMatrix4 prod;
    Complex r0[4] = {a.r0c0,a.r0c1,a.r0c2,a.r0c3};
    Complex r1[4] = {a.r1c0,a.r1c1,a.r1c2,a.r1c3};
    Complex r2[4] = {a.r2c0,a.r2c1,a.r2c2,a.r2c3};
    Complex r3[4] = {a.r3c0,a.r3c1,a.r3c2,a.r3c3};
    Complex c0[4] = {b.r0c0,b.r1c0,b.r2c0,b.r3c0};
    Complex c1[4] = {b.r0c1,b.r1c1,b.r2c1,b.r3c1};
    Complex c2[4] = {b.r0c2,b.r1c2,b.r2c2,b.r3c2};
    Complex c3[4] = {b.r0c3,b.r1c3,b.r2c3,b.r3c3};
    
    prod.r0c0 = getMatrixProductElement(r0,c0,4);
    prod.r0c1 = getMatrixProductElement(r0,c1,4);
    prod.r0c2 = getMatrixProductElement(r0,c2,4);
    prod.r0c3 = getMatrixProductElement(r0,c3,4);
    
    prod.r1c0 = getMatrixProductElement(r1,c0,4);
    prod.r1c1 = getMatrixProductElement(r1,c1,4);
    prod.r1c2 = getMatrixProductElement(r1,c2,4);
    prod.r1c3 = getMatrixProductElement(r1,c3,4);
    
    prod.r2c0 = getMatrixProductElement(r2,c0,4);
    prod.r2c1 = getMatrixProductElement(r2,c1,4);
    prod.r2c2 = getMatrixProductElement(r2,c2,4);
    prod.r2c3 = getMatrixProductElement(r2,c3,4);
    
    prod.r3c0 = getMatrixProductElement(r3,c0,4);
    prod.r3c1 = getMatrixProductElement(r3,c1,4);
    prod.r3c2 = getMatrixProductElement(r3,c2,4);
    prod.r3c3 = getMatrixProductElement(r3,c3,4);
    
    return prod;
}

ComplexMatrix2 getConjugateTransposeMatrix2(ComplexMatrix2 u) {
    ComplexMatrix2 c = getConjugateMatrix2(u);
    ComplexMatrix2 t;
    t.r0c0=c.r0c0;  t.r0c1=c.r1c0;
    t.r1c0=c.r0c1;  t.r1c1=c.r1c1;
    return t;
}

ComplexMatrix4 getConjugateTransposeMatrix4(ComplexMatrix4 u) {
    ComplexMatrix4 c = getConjugateMatrix4(u);
    ComplexMatrix4 t;
    t.r0c0=c.r0c0;  t.r0c1=c.r1c0;  t.r0c2=c.r2c0;  t.r0c3=c.r3c0;
    t.r1c0=c.r0c1;  t.r1c1=c.r1c1;  t.r1c2=c.r2c1;  t.r1c3=c.r3c1;
    t.r2c0=c.r0c2;  t.r2c1=c.r1c2;  t.r2c2=c.r2c2;  t.r2c3=c.r3c2;
    t.r3c0=c.r0c3;  t.r3c1=c.r1c3;  t.r3c2=c.r2c3;  t.r3c3=c.r3c3;
    return t;
}

void addToMatrix2(ComplexMatrix2* dest, ComplexMatrix2 add) {
    dest->r0c0.real += add.r0c0.real; dest->r0c0.imag += add.r0c0.imag;
    dest->r0c1.real += add.r0c1.real; dest->r0c1.imag += add.r0c1.imag;
    dest->r1c0.real += add.r1c0.real; dest->r1c0.imag += add.r1c0.imag;
    dest->r1c1.real += add.r1c1.real; dest->r1c1.imag += add.r1c1.imag;
}

void addToMatrix4(ComplexMatrix4* dest, ComplexMatrix4 add) {
    dest->r0c0.real += add.r0c0.real; dest->r0c0.imag += add.r0c0.imag;
    dest->r0c1.real += add.r0c1.real; dest->r0c1.imag += add.r0c1.imag;
    dest->r0c2.real += add.r0c2.real; dest->r0c2.imag += add.r0c2.imag;
    dest->r0c3.real += add.r0c3.real; dest->r0c3.imag += add.r0c3.imag;
    dest->r1c0.real += add.r1c0.real; dest->r1c0.imag += add.r1c0.imag;
    dest->r1c1.real += add.r1c1.real; dest->r1c1.imag += add.r1c1.imag;
    dest->r1c2.real += add.r1c2.real; dest->r1c2.imag += add.r1c2.imag;
    dest->r1c3.real += add.r1c3.real; dest->r1c3.imag += add.r1c3.imag;
    dest->r2c0.real += add.r2c0.real; dest->r2c0.imag += add.r2c0.imag;
    dest->r2c1.real += add.r2c1.real; dest->r2c1.imag += add.r2c1.imag;
    dest->r2c2.real += add.r2c2.real; dest->r2c2.imag += add.r2c2.imag;
    dest->r2c3.real += add.r2c3.real; dest->r2c3.imag += add.r2c3.imag;
    dest->r3c0.real += add.r3c0.real; dest->r3c0.imag += add.r3c0.imag;
    dest->r3c1.real += add.r3c1.real; dest->r3c1.imag += add.r3c1.imag;
    dest->r3c2.real += add.r3c2.real; dest->r3c2.imag += add.r3c2.imag;
    dest->r3c3.real += add.r3c3.real; dest->r3c3.imag += add.r3c3.imag;
}

/* returns |a - b|^2 */ 
qreal getComplexDist(Complex a, Complex b) {
    qreal reDif = absReal(a.real - b.real);
    qreal imDif = absReal(a.imag - b.imag);
    
    return reDif*reDif + imDif*imDif;
}

qreal getHilbertSchmidtDistFromIdentity2(ComplexMatrix2 a) {
    ComplexMatrix2 iden = {0};
    iden.r0c0.real=1; iden.r1c1.real=1;

    qreal dist = 
        getComplexDist(a.r0c0, iden.r0c0) +
        getComplexDist(a.r0c1, iden.r0c1) +
        getComplexDist(a.r1c0, iden.r1c0) +
        getComplexDist(a.r1c1, iden.r1c1);
    return dist;
}

qreal getHilbertSchmidtDistFromIdentity4(ComplexMatrix4 a) {
    ComplexMatrix4 iden = {0};
    iden.r0c0.real=1; iden.r1c1.real=1; iden.r2c2.real=1; iden.r3c3.real=1;
    
    qreal dist = 
        getComplexDist(a.r0c0, iden.r0c0) + 
        getComplexDist(a.r0c1, iden.r0c1) + 
        getComplexDist(a.r0c2, iden.r0c2) + 
        getComplexDist(a.r0c3, iden.r0c3) + 
        getComplexDist(a.r1c0, iden.r1c0) + 
        getComplexDist(a.r1c1, iden.r1c1) + 
        getComplexDist(a.r1c2, iden.r1c2) + 
        getComplexDist(a.r1c3, iden.r1c3) + 
        getComplexDist(a.r2c0, iden.r2c0) + 
        getComplexDist(a.r2c1, iden.r2c1) + 
        getComplexDist(a.r2c2, iden.r2c2) + 
        getComplexDist(a.r2c3, iden.r2c3) + 
        getComplexDist(a.r3c0, iden.r3c0) + 
        getComplexDist(a.r3c1, iden.r3c1) + 
        getComplexDist(a.r3c2, iden.r3c2) + 
        getComplexDist(a.r3c3, iden.r3c3);
    return dist;
}

int isMatrix4Unitary(ComplexMatrix4 u) {
    
    // make I = U U^dagger
    ComplexMatrix4 identity = getMatrix4Product(u, getConjugateTransposeMatrix4(u));

    // diagonals must be 1 (real)
    if (absReal(identity.r0c0.real - 1) > REAL_EPS) return 0;
    if (absReal(identity.r0c0.imag    ) > REAL_EPS) return 0;
    if (absReal(identity.r1c1.real - 1) > REAL_EPS) return 0;
    if (absReal(identity.r1c1.imag    ) > REAL_EPS) return 0;
    if (absReal(identity.r2c2.real - 1) > REAL_EPS) return 0;
    if (absReal(identity.r2c2.imag    ) > REAL_EPS) return 0;
    if (absReal(identity.r3c3.real - 1) > REAL_EPS) return 0;
    if (absReal(identity.r3c3.imag    ) > REAL_EPS) return 0;
    
    // off diagonals must be 0
    if (absReal(identity.r0c1.real) > REAL_EPS) return 0;
    if (absReal(identity.r0c1.imag) > REAL_EPS) return 0;
    if (absReal(identity.r0c2.real) > REAL_EPS) return 0;
    if (absReal(identity.r0c2.imag) > REAL_EPS) return 0;
    if (absReal(identity.r0c3.real) > REAL_EPS) return 0;
    if (absReal(identity.r0c3.imag) > REAL_EPS) return 0;
    
    if (absReal(identity.r1c0.real) > REAL_EPS) return 0;
    if (absReal(identity.r1c0.imag) > REAL_EPS) return 0;
    if (absReal(identity.r1c2.real) > REAL_EPS) return 0;
    if (absReal(identity.r1c2.imag) > REAL_EPS) return 0;
    if (absReal(identity.r1c3.real) > REAL_EPS) return 0;
    if (absReal(identity.r1c3.imag) > REAL_EPS) return 0;
    
    if (absReal(identity.r2c0.real) > REAL_EPS) return 0;
    if (absReal(identity.r2c0.imag) > REAL_EPS) return 0;
    if (absReal(identity.r2c1.real) > REAL_EPS) return 0;
    if (absReal(identity.r2c1.imag) > REAL_EPS) return 0;
    if (absReal(identity.r2c3.real) > REAL_EPS) return 0;
    if (absReal(identity.r2c3.imag) > REAL_EPS) return 0;
    
    if (absReal(identity.r3c0.real) > REAL_EPS) return 0;
    if (absReal(identity.r3c0.imag) > REAL_EPS) return 0;
    if (absReal(identity.r3c1.real) > REAL_EPS) return 0;
    if (absReal(identity.r3c1.imag) > REAL_EPS) return 0;
    if (absReal(identity.r3c2.real) > REAL_EPS) return 0;
    if (absReal(identity.r3c2.imag) > REAL_EPS) return 0;
    
    return 1;
}

int isMatrixNUnitary(ComplexMatrixN u) {
    long long int r, c, i, dim;
    dim = u.numRows;
    
    // check u * ConjTrans(u) = Identity
    for (r=0; r < dim; r++) {
        for (c=0; c < dim; c++) {
            
            // u[r][...] * ConjTrans(u)[...][c]
            qreal elemRe = 0;
            qreal elemIm = 0;
            for (i=0; i < dim; i++) {
                // u.elems[r][i] * conj(u.elems[c][i])
                elemRe += u.elems[r][i].real*u.elems[c][i].real + u.elems[r][i].imag*u.elems[c][i].imag;
                elemIm += u.elems[r][i].imag*u.elems[c][i].real - u.elems[r][i].real*u.elems[c][i].imag;
            }
            
            if (absReal(elemIm) > REAL_EPS)
                return 0;
            if (r == c && absReal(elemRe - 1) > REAL_EPS)
                return 0;
            if (r != c && absReal(elemRe    ) > REAL_EPS)
                return 0;
        }
    }
    
    return 1;
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

void validateNumTargets(Qureg qureg, const int numTargetQubits, const char* caller) {
    QuESTAssert(numTargetQubits>0 && numTargetQubits<=qureg.numQubitsRepresented, E_INVALID_NUM_TARGETS, caller);
}

void validateNumControls(Qureg qureg, const int numControlQubits, const char* caller) {
    QuESTAssert(numControlQubits>0 && numControlQubits<=qureg.numQubitsRepresented, E_INVALID_NUM_CONTROLS, caller);
}

void validateMultiTargets(Qureg qureg, int* targetQubits, const int numTargetQubits, const char* caller) {
    validateNumTargets(qureg, numTargetQubits, caller);
    for (int i=0; i < numTargetQubits; i++) 
        validateTarget(qureg, targetQubits[i], caller);
        
    QuESTAssert(areUniqueQubits(targetQubits, numTargetQubits), E_TARGETS_NOT_UNIQUE, caller);
}

void validateMultiControls(Qureg qureg, int* controlQubits, const int numControlQubits, const char* caller) {
    validateNumControls(qureg, numControlQubits, caller);
    for (int i=0; i < numControlQubits; i++)
        validateControl(qureg, controlQubits[i], caller);
        
    QuESTAssert(areUniqueQubits(controlQubits, numControlQubits), E_CONTROLS_NOT_UNIQUE, caller);
}

void validateMultiControlsTarget(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, const char* caller) {
    validateTarget(qureg, targetQubit, caller);
    validateMultiControls(qureg, controlQubits, numControlQubits, caller);
    for (int i=0; i < numControlQubits; i++)
        QuESTAssert(controlQubits[i] != targetQubit, E_TARGET_IN_CONTROLS, caller);
}

void validateMultiControlsMultiTargets(Qureg qureg, int* controlQubits, const int numControlQubits, int* targetQubits, const int numTargetQubits, const char* caller) {
    validateMultiControls(qureg, controlQubits, numControlQubits, caller);
    validateMultiTargets(qureg, targetQubits, numTargetQubits, caller);
    long long int ctrlMask = getQubitBitMask(controlQubits, numControlQubits);
    long long int targMask = getQubitBitMask(targetQubits, numTargetQubits);
    int overlap = ctrlMask & targMask;
    QuESTAssert(!overlap, E_TARGET_IN_CONTROLS, caller);
}

void validateControlState(int* controlState, const int numControlQubits, const char* caller) {
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
    QuESTAssert(matr.elems != NULL, E_COMPLEX_MATRIX_NOT_INIT, caller);
}

void validateMultiQubitUnitaryMatrix(Qureg qureg, ComplexMatrixN u, int numTargs, const char* caller) { 
    validateMatrixInit(u, caller);
    validateMultiQubitMatrixFitsInNode(qureg, numTargs, caller);
    QuESTAssert(numTargs == u.numQubits, E_INVALID_UNITARY_SIZE, caller);
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
        int code = pauliCodes[i];
        QuESTAssert(
            code==PAULI_I || code==PAULI_X || code==PAULI_Y || code==PAULI_Z, 
            E_INVALID_PAULI_CODE, caller);
    }
}

void validateNumSumTerms(int numTerms, const char* caller) {
    QuESTAssert(numTerms > 0, E_INVALID_NUM_SUM_TERMS, caller);
}

void validateOneQubitKrausMap(Qureg qureg, ComplexMatrix2* ops, int numOps, const char* caller) {
    validateMultiQubitMatrixFitsInNode(qureg, 2, caller);
    QuESTAssert(numOps > 0 && numOps < 4, E_INVALID_NUM_ONE_QUBIT_KRAUS_OPS, caller);
    
    // sum of conjTrans(op) * op
    ComplexMatrix2 sumConjProd = {0};
    for (int n=0; n < numOps; n++)
        addToMatrix2(&sumConjProd, getMatrix2Product(getConjugateTransposeMatrix2(ops[n]), ops[n]));
    
    qreal idenDist = getHilbertSchmidtDistFromIdentity2(sumConjProd);
    QuESTAssert(idenDist < REAL_EPS, E_INVALID_KRAUS_OPS, caller);
}

void validateTwoQubitKrausMap(Qureg qureg, ComplexMatrix4* ops, int numOps, const char* caller) {
    validateMultiQubitMatrixFitsInNode(qureg, 4, caller);
    QuESTAssert(numOps > 0 && numOps < 16, E_INVALID_NUM_TWO_QUBIT_KRAUS_OPS, caller);
    
    // sum of conjTrans(op) * op
    ComplexMatrix4 sumConjProd = {0};
    for (int n=0; n < numOps; n++)
        addToMatrix4(&sumConjProd, getMatrix4Product(getConjugateTransposeMatrix4(ops[n]), ops[n]));
        
    qreal idenDist = getHilbertSchmidtDistFromIdentity4(sumConjProd);
    QuESTAssert(idenDist < REAL_EPS, E_INVALID_KRAUS_OPS, caller);
}


#ifdef __cplusplus
}
#endif
