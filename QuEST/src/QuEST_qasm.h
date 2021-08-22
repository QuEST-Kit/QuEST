// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Functions for generating QASM output from QuEST circuits
 *
 * @author Tyson Jones
 */

# ifndef QUEST_QASM_H
# define QUEST_QASM_H

# include "QuEST.h"
# include "QuEST_precision.h"

# ifdef __cplusplus
extern "C" {
# endif

/**! Identifiers of single-target gates */
typedef enum {
    GATE_SIGMA_X,
    GATE_SIGMA_Y,
    GATE_SIGMA_Z,
    GATE_T,
    GATE_S,
    GATE_HADAMARD,
    GATE_ROTATE_X,
    GATE_ROTATE_Y,
    GATE_ROTATE_Z,
    GATE_ROTATE_AROUND_AXIS,
    GATE_UNITARY,
    GATE_PHASE_SHIFT,
    GATE_SWAP,
    GATE_SQRT_SWAP
} TargetGate;

void qasm_setup(Qureg* qureg);

void qasm_startRecording(Qureg qureg);

void qasm_stopRecording(Qureg qureg);

void qasm_recordGate(Qureg qureg, TargetGate gate, int targetQubit);

void qasm_recordParamGate(Qureg qureg, TargetGate gate, int targetQubit, qreal param);

void qasm_recordCompactUnitary(Qureg qureg, Complex alpha, Complex beta, int targetQubit);

void qasm_recordUnitary(Qureg qureg, ComplexMatrix2 u, int targetQubit);

void qasm_recordAxisRotation(Qureg qureg, qreal angle, Vector axis, int targetQubit);

void qasm_recordControlledGate(Qureg qureg, TargetGate gate, int controlQubit, int targetQubit);

void qasm_recordControlledParamGate(Qureg qureg, TargetGate gate, int controlQubit, int targetQubit, qreal param);

void qasm_recordControlledCompactUnitary(Qureg qureg, Complex alpha, Complex beta, int controlQubit, int targetQubit);

void qasm_recordControlledUnitary(Qureg qureg, ComplexMatrix2 u, int controlQubit, int targetQubit);

void qasm_recordControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int controlQubit, int targetQubit);

void qasm_recordMultiControlledGate(Qureg qureg,
TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit);

void qasm_recordMultiControlledParamGate(Qureg qureg, TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit, qreal param);

void qasm_recordMultiControlledUnitary(Qureg qureg, ComplexMatrix2 u, int* controlQubits, int numControlQubits, int targetQubit);

void qasm_recordMultiStateControlledUnitary(Qureg qureg, ComplexMatrix2 u, int* controlQubits, int* controlState, int numControlQubits, int targetQubit);

void qasm_recordMultiControlledMultiQubitNot(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs);

/* not actually used. D'oh!
void qasm_recordMultiControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int* controlQubits, int numControlQubits, int targetQubit);\
*/

void qasm_recordMeasurement(Qureg qureg, int measureQubit);

void qasm_recordInitZero(Qureg qureg);

void qasm_recordInitPlus(Qureg qureg);

void qasm_recordInitClassical(Qureg qureg, long long int stateInd);

void qasm_recordPhaseFunc(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, qreal* overridePhases, int numOverrides);

void qasm_recordMultiVarPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg, long long int* overrideInds, qreal* overridePhases, int numOverrides);

void qasm_recordNamedPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams, long long int* overrideInds, qreal* overridePhases, int numOverrides);

void qasm_recordComment(Qureg qureg, char* comment, ...);

void qasm_clearRecorded(Qureg qureg);

void qasm_printRecorded(Qureg qureg);

int qasm_writeRecordedToFile(Qureg qureg, char* filename);

void qasm_free(Qureg qureg);

# ifdef __cplusplus
}
# endif

# endif // QUEST_QASM_H
