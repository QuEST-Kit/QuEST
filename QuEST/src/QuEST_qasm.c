// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Functions for generating QASM output from QuEST circuits
 */

/** @TODO
 * - allow user-set decimal precision (useful for when QASM is passed to a plotter)
 * - sort out fixing global phase in controlledPhaseShift to controlledRotateZ plug
 * - add functions to directly add comments to QASM by user
 * - add abilitiy for user to directly add strings to QASM buffer??
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"
# include "QuEST_qasm.h"

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# define QUREG_LABEL "q"        // QASM var-name for the quantum register
# define MESREG_LABEL "c"       // QASM var-name for the classical measurement register
# define CTRL_LABEL_PREF "c"    // QASM syntax which prefixes gates when controlled
# define MEASURE_CMD "measure"  // QASM cmd for measurement operation
# define INIT_ZERO_CMD "reset"  // QASM cmd for setting state 0
# define COMMENT_PREF "//"     // QASM syntax for a comment ;)

# define MAX_LINE_LEN 200       // maximum length (#chars) of a single QASM instruction
# define BUF_INIT_SIZE 1000     // initial size of the QASM buffer (#chars)
# define BUF_GROW_FAC 2         // growth factor when buffer dynamically resizes

static const char* qasmGateLabels[] = {
    [GATE_SIGMA_X] = "x",
    [GATE_SIGMA_Y] = "y",
    [GATE_SIGMA_Z] = "z",
    [GATE_T] = "t",
    [GATE_S] = "s",
    [GATE_HADAMARD] = "h",
    [GATE_ROTATE_X] = "Rx",
    [GATE_ROTATE_Y] = "Ry",
    [GATE_ROTATE_Z] = "Rz",
    [GATE_UNITARY] = "U",     // needs phase fix when controlled
    [GATE_PHASE_SHIFT] = "Rz" // needs phase fix when controlled
};

// @TODO make a proper internal error thing
void bufferOverflow() {
    printf("!!!\nINTERNAL ERROR: QASM line buffer filled!\n!!!");
    exit(1);
}

void qasm_setup(Qureg* qureg) {
    
    // populate and attach QASM logger
    QASMLogger *qasmLog = malloc(sizeof *qasmLog);
    qureg->qasmLog = qasmLog;
    if (qasmLog == NULL)
        bufferOverflow();
    
    qasmLog->isLogging = 0;
    qasmLog->bufferSize = BUF_INIT_SIZE;
    qasmLog->buffer = malloc(qasmLog->bufferSize * sizeof *(qasmLog->buffer));
    if (qasmLog->buffer == NULL)
        bufferOverflow();
    
    // add headers and quantum / classical register creation
    qasmLog->bufferFill = snprintf(
        qasmLog->buffer, qasmLog->bufferSize,
        "OPENQASM 2.0;\nqreg %s[%d];\ncreg %s[%d];\n", 
        QUREG_LABEL, qureg->numQubitsRepresented,
        MESREG_LABEL, qureg->numQubitsRepresented);
    if (qasmLog->bufferFill >= qasmLog->bufferSize)
        bufferOverflow();
}

void qasm_startRecording(Qureg qureg) {
    qureg.qasmLog->isLogging = 1;
}

void qasm_stopRecording(Qureg qureg) {
    qureg.qasmLog->isLogging = 0;
}

void addStringToQASM(Qureg qureg, char line[], int lineLen) {
    
    char* buf = qureg.qasmLog->buffer;
    int bufSize = qureg.qasmLog->bufferSize;
    int bufFill = qureg.qasmLog->bufferFill;
    
    // grow QASM buffer if necessary
    if (lineLen + bufFill > bufSize) {
                
        int newBufSize = BUF_GROW_FAC * bufSize;
        if (lineLen + bufFill > newBufSize)
            bufferOverflow();
        
        char* newBuffer = malloc(newBufSize * sizeof *newBuffer);
        sprintf(newBuffer, "%s", buf);
        free(buf);
        
        qureg.qasmLog->bufferSize = newBufSize;
        qureg.qasmLog->buffer = newBuffer;
        bufSize = newBufSize;
        buf = newBuffer;
    }
        
    // add new str
    int addedChars = snprintf(buf+bufFill, bufSize-bufFill, "%s", line);
    qureg.qasmLog->bufferFill += addedChars;
}

void qasm_recordComment(Qureg qureg, char* comment) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    int len = snprintf(line, MAX_LINE_LEN, "%s %s\n", COMMENT_PREF, comment);
    addStringToQASM(qureg, line, len);
}

void addGateToQASM(Qureg qureg, TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit, qreal* params, int numParams) {
    
    int len = 0;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    
    // add control labels
    for (int i=0; i < numControlQubits; i++)
        len += snprintf(line+len, MAX_LINE_LEN-len, "%s", CTRL_LABEL_PREF);
    
    // add target gate
    len += snprintf(line+len, MAX_LINE_LEN-len, "%s", qasmGateLabels[gate]);
    
    // add parameters
    if (numParams > 0) {
        len += snprintf(line+len, MAX_LINE_LEN-len, "(");
        for (int i=0; i < numParams; i++) {
            len += snprintf(line+len, MAX_LINE_LEN-len, REAL_QASM_FORMAT, params[i]);
            if (i != numParams - 1)
                len += snprintf(line+len, MAX_LINE_LEN-len, ",");
        }
        len += snprintf(line+len, MAX_LINE_LEN-len, ")");
    }
    
    // add space
    len += snprintf(line+len, MAX_LINE_LEN-len, " ");
    
    // add control qubits
    for (int i=0; i < numControlQubits; i++)
        len += snprintf(line+len, MAX_LINE_LEN-len, "%s[%d],", QUREG_LABEL, controlQubits[i]);
    
    // add target qubit, colon and newline
    len += snprintf(line+len, MAX_LINE_LEN-len, "%s[%d];\n", QUREG_LABEL, targetQubit);
    
    // check whether we overflowed buffer
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
        
    addStringToQASM(qureg, line, len);
}

void qasm_recordGate(Qureg qureg, TargetGate gate, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, NULL, 0, targetQubit, NULL, 0);
}

void qasm_recordParamGate(Qureg qureg, TargetGate gate, int targetQubit, qreal param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    qreal params[1] = {param};
    addGateToQASM(qureg, gate, NULL, 0, targetQubit, params, 1);
}

void qasm_recordCompactUnitary(Qureg qureg, Complex alpha, Complex beta, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, NULL, 0, targetQubit, params, 3);
}

void qasm_recordUnitary(Qureg qureg, ComplexMatrix2 u, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    qreal discardedGlobalPhase;
    getComplexPairAndPhaseFromUnitary(u, &alpha, &beta, &discardedGlobalPhase);
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, NULL, 0, targetQubit, params, 3);
}

void qasm_recordAxisRotation(Qureg qureg, qreal angle, Vector axis, const int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, NULL, 0, targetQubit, params, 3);
}

void qasm_recordControlledGate(Qureg qureg, TargetGate gate, int controlQubit, int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int controls[1] = {controlQubit};
    addGateToQASM(qureg, gate, controls, 1, targetQubit, 0, 0);    
}

void qasm_recordControlledParamGate(Qureg qureg, TargetGate gate, int controlQubit, int targetQubit, qreal param) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int controls[1] = {controlQubit};
    qreal params[1] = {param};
    addGateToQASM(qureg, gate, controls, 1, targetQubit, params, 1);
    
    // correct the global phase of controlled phase shifts
    if (gate == GATE_PHASE_SHIFT) {
        qasm_recordComment(qureg, "Restoring the discarded global phase of the previous controlled phase gate");
        qreal phaseFix[1] = {param/2.0};
        addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
    }
}

void qasm_recordControlledCompactUnitary(Qureg qureg, Complex alpha, Complex beta, int controlQubit, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    int controls[1] = {controlQubit};
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controls, 1, targetQubit, params, 3);
}

/** additionally performs Rz on target to restore the global phase lost from u in QASM U(a,b,c) */
void qasm_recordControlledUnitary(Qureg qureg, ComplexMatrix2 u, int controlQubit, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    qreal globalPhase;
    getComplexPairAndPhaseFromUnitary(u, &alpha, &beta, &globalPhase);
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    int controls[1] = {controlQubit};
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controls, 1, targetQubit, params, 3);
    
    // add Rz
    qasm_recordComment(qureg, "Restoring the discarded global phase of the previous controlled unitary");
    qreal phaseFix[1] = {globalPhase};
    addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
}

void qasm_recordControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int controlQubit, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    int controls[1] = {controlQubit};
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controls, 1, targetQubit, params, 3);
}

void qasm_recordMultiControlledGate(Qureg qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, NULL, 0);
}

void qasm_recordMultiControlledParamGate(Qureg qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit, qreal param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    qreal params[1] = {param};
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, params, 1);
    
    // correct the global phase of controlled phase shifts
    if (gate == GATE_PHASE_SHIFT) {
        qasm_recordComment(qureg, "Restoring the discarded global phase of the previous multicontrolled phase gate");
        qreal phaseFix[1] = {param/2.0};
        addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
    }
}

/** additionally performs Rz on target to restore the global phase lost from u in QASM U(a,b,c) */
void qasm_recordMultiControlledUnitary(Qureg qureg, ComplexMatrix2 u, int* controlQubits, const int numControlQubits, const int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    qreal globalPhase;
    getComplexPairAndPhaseFromUnitary(u, &alpha, &beta, &globalPhase);
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controlQubits, numControlQubits, targetQubit, params, 3);
    
    // add Rz
    qreal phaseFix[1] = {globalPhase};
    addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
}

/* not actually used, D'Oh!
void qasm_recordMultiControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int* controlQubits, const int numControlQubits, const int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    
    qreal rz2, ry, rz1;
    getZYZRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    qreal params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controlQubits, numControlQubits, targetQubit, params, 3);
}
*/

void qasm_recordMeasurement(Qureg qureg, const int measureQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    int len = snprintf(
        line, MAX_LINE_LEN, "%s %s[%d] -> %s[%d];\n",
        MEASURE_CMD, QUREG_LABEL, measureQubit, MESREG_LABEL, measureQubit);
        
    // check whether we overflowed buffer
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
    
    addStringToQASM(qureg, line, len);
}

void qasm_recordInitZero(Qureg qureg) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    int len = snprintf(line, MAX_LINE_LEN, "%s %s;\n", INIT_ZERO_CMD, QUREG_LABEL);
    
    // check whether we overflowed buffer
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
    
    addStringToQASM(qureg, line, len);
}

void qasm_recordInitPlus(Qureg qureg) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    // add an explanatory comment
    char buf[MAX_LINE_LEN+1];
    sprintf(buf, "Initialising state |+>");
    qasm_recordComment(qureg, buf);
    
    // it's valid QASM to h the register (I think)
    // |+> = H |0>
    qasm_recordInitZero(qureg);
    int charsWritten = snprintf(
        buf, MAX_LINE_LEN, "%s %s;\n", 
        qasmGateLabels[GATE_HADAMARD], QUREG_LABEL);
    if (charsWritten >= MAX_LINE_LEN)
        bufferOverflow();
    addStringToQASM(qureg, buf, charsWritten);
    
    // old code (before above QASM shortcut)
    /*
    qasm_recordInitZero(qureg);
    for (int q=0; q < qureg.numQubitsRepresented; q++)
        qasm_recordGate(qureg, GATE_HADAMARD, q);
    */
}

void qasm_recordInitClassical(Qureg qureg, long long int stateInd) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    // add an explanatory comment
    char cmt[MAX_LINE_LEN+1];
    sprintf(cmt, "Initialising state |%lld>", stateInd);
    qasm_recordComment(qureg, cmt);
    
    // start in |0>
    qasm_recordInitZero(qureg);
    
    // NOT the 1 bits in stateInd
    for (int q=0; q < qureg.numQubitsRepresented; q++) 
        if ((stateInd >> q) & 1)
            qasm_recordGate(qureg, GATE_SIGMA_X, q);
}

void qasm_clearRecorded(Qureg qureg) {
    
    // maintains current buffer size
    (qureg.qasmLog->buffer)[0] = '\0';
    qureg.qasmLog->bufferFill = 0;
}

void qasm_printRecorded(Qureg qureg) {
    printf("%s", qureg.qasmLog->buffer);
}

/** returns success of file write */
int qasm_writeRecordedToFile(Qureg qureg, char* filename) {
    
    FILE *file = fopen(filename, "w");
    if (file == NULL)
        return 0;
    
    fprintf(file, "%s", qureg.qasmLog->buffer);
    fclose(file);
    return 1;
}

void qasm_free(Qureg qureg) {
    
    free(qureg.qasmLog->buffer);
    free(qureg.qasmLog);
}
