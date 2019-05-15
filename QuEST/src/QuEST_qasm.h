// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Functions for generating QASM output from QuEST circuits
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
    GATE_PHASE_SHIFT
} TargetGate;

void qasm_setup(Qureg* qureg);

void qasm_startRecording(Qureg qureg);

void qasm_stopRecording(Qureg qureg);

void qasm_recordGate(Qureg qureg, TargetGate gate, int targetQubit);

void qasm_recordParamGate(Qureg qureg, TargetGate gate, int targetQubit, qreal param);

void qasm_recordCompactUnitary(Qureg qureg, Complex alpha, Complex beta, int targetQubit);

void qasm_recordUnitary(Qureg qureg, ComplexMatrix2 u, int targetQubit);

void qasm_recordAxisRotation(Qureg qureg, qreal angle, Vector axis, const int targetQubit);

void qasm_recordControlledGate(Qureg qureg, TargetGate gate, int controlQubit, int targetQubit);

void qasm_recordControlledParamGate(Qureg qureg, TargetGate gate, int controlQubit, int targetQubit, qreal param);

void qasm_recordControlledCompactUnitary(Qureg qureg, Complex alpha, Complex beta, int controlQubit, int targetQubit);

void qasm_recordControlledUnitary(Qureg qureg, ComplexMatrix2 u, int controlQubit, int targetQubit);

void qasm_recordControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int controlQubit, int targetQubit);

void qasm_recordMultiControlledGate(Qureg qureg,
TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit);

void qasm_recordMultiControlledParamGate(Qureg qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit, qreal param);

void qasm_recordMultiControlledUnitary(Qureg qureg, ComplexMatrix2 u, int* controlQubits, const int numControlQubits, const int targetQubit);

/* not actually used. D'oh!
void qasm_recordMultiControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int* controlQubits, const int numControlQubits, const int targetQubit);\
*/

void qasm_recordMeasurement(Qureg qureg, const int measureQubit);

void qasm_recordInitZero(Qureg qureg);

void qasm_recordInitPlus(Qureg qureg);

void qasm_recordInitClassical(Qureg qureg, long long int stateInd);

void qasm_recordComment(Qureg qureg, char* comment);

void qasm_clearRecorded(Qureg qureg);

void qasm_printRecorded(Qureg qureg);

int qasm_writeRecordedToFile(Qureg qureg, char* filename);

void qasm_free(Qureg qureg);

# ifdef __cplusplus
}
# endif

# endif // QUEST_QASM_H
