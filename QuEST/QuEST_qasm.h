// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Functions for generating QASM output from QuEST circuits
 */

# ifndef QuEST_QASM
# define QuEST_QASM

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

void qasm_setup(QubitRegister* qureg);

void qasm_startRecording(QubitRegister qureg);

void qasm_stopRecording(QubitRegister qureg);

void qasm_recordGate(QubitRegister qureg, TargetGate gate, int targetQubit);

void qasm_recordParamGate(QubitRegister qureg, TargetGate gate, int targetQubit, REAL param);

void qasm_recordCompactUnitary(QubitRegister qureg, Complex alpha, Complex beta, int targetQubit);

void qasm_recordUnitary(QubitRegister qureg, ComplexMatrix2 u, int targetQubit);

void qasm_recordAxisRotation(QubitRegister qureg, REAL angle, Vector axis, const int targetQubit);

void qasm_recordControlledGate(QubitRegister qureg, TargetGate gate, int controlQubit, int targetQubit);

void qasm_recordControlledParamGate(QubitRegister qureg, TargetGate gate, int controlQubit, int targetQubit, REAL param);

void qasm_recordControlledCompactUnitary(QubitRegister qureg, Complex alpha, Complex beta, int controlQubit, int targetQubit);

void qasm_recordControlledUnitary(QubitRegister qureg, ComplexMatrix2 u, int controlQubit, int targetQubit);

void qasm_recordControlledAxisRotation(QubitRegister qureg, REAL angle, Vector axis, int controlQubit, int targetQubit);

void qasm_recordMultiControlledGate(QubitRegister qureg,
TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit);

void qasm_recordMultiControlledParamGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit, REAL param);

void qasm_recordMultiControlledUnitary(QubitRegister qureg, ComplexMatrix2 u, int* controlQubits, const int numControlQubits, const int targetQubit);

/* not actually used. D'oh!
void qasm_recordMultiControlledAxisRotation(QubitRegister qureg, REAL angle, Vector axis, int* controlQubits, const int numControlQubits, const int targetQubit);\
*/

void qasm_recordMeasurement(QubitRegister qureg, const int measureQubit);

void qasm_recordInitZero(QubitRegister qureg);

void qasm_recordInitPlus(QubitRegister qureg);

void qasm_recordInitClassical(QubitRegister qureg, long long int stateInd);

void qasm_recordComment(QubitRegister qureg, char* comment);

void qasm_clearRecorded(QubitRegister qureg);

void qasm_printRecorded(QubitRegister qureg);

int qasm_writeRecordedToFile(QubitRegister qureg, char* filename);

void qasm_free(QubitRegister qureg);

# ifdef __cplusplus
}
# endif

# endif // QuEST_QASM