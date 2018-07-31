// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Functions for generating QASM output from QuEST circuits
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_qasm.h"

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# define QUREG_LABEL "q"        // QASM var-name for the quantum register
# define MESREG_LABEL "c"       // QASM var-name for the classical measurement register
# define CTRL_LABEL_PREF "c"    // QASM syntax which prefixes gates when controlled
# define MEASURE_CMD "measure"  // QASM label for measurement operation

# define MAX_LINE_LEN 100       // maximum length (#chars) of a single QASM instruction
# define BUF_INIT_SIZE 1000     // initial size of the QASM buffer (#chars)
# define BUF_GROW_FAC 2         // growth factor when buffer dynamically resizes

/* @TODO
    - shrink param strings where possible
*/

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
    [GATE_UNITARY] = "U"
    //[GATE_ROTATE_AROUND_AXIS] = ,
    //[GATE_PHASE_SHIFT] = 
};

void qasm_setup(QubitRegister* qureg) {
    
    // populate and attach QASM logger
    QASMLogger *qasmLog = malloc(sizeof qasmLog);
    qasmLog->isLogging = 0;
    qasmLog->bufferSize = BUF_INIT_SIZE;
    qasmLog->buffer = malloc(qasmLog->bufferSize * sizeof *(qasmLog->buffer));
    qasmLog->bufferFill = sprintf(qasmLog->buffer, "OPENQASM 2.0;\nqreg %s[%d];\ncreg %s[%d];\n", 
        QUREG_LABEL, qureg->numQubitsRepresented,
        MESREG_LABEL, qureg->numQubitsRepresented);

    qureg->qasmLog = qasmLog;
}

void qasm_startRecording(QubitRegister qureg) {
    qureg.qasmLog->isLogging = 1;
}

void qasm_stopRecording(QubitRegister qureg) {
    qureg.qasmLog->isLogging = 0;
}

// @TODO make a proper internal error thing
void bufferOverflow() {
    printf("!!!\nINTERNAL ERROR: QASM line buffer filled!\n!!!");
    exit(1);
}

/** maps U(alpha, beta) to Rz(rz2) Ry(ry) Rz(rz1) */
void getRotAnglesFromComplexPair(Complex alpha, Complex beta, REAL* rz2, REAL* ry, REAL* rz1) {
    
    REAL alphaMag = sqrt(alpha.real*alpha.real + alpha.imag*alpha.imag);
    *ry = 2.0 * acos(alphaMag);
    
    REAL alphaPhase = atan2(alpha.imag, alpha.real);
    REAL betaPhase  = atan2(beta.imag,  beta.real);
    *rz2 = - alphaPhase + betaPhase;
    *rz1 = - alphaPhase - betaPhase;
}

/** maps U(r0c0, r0c1, r1c0, r1c1) to exp(i globalPhase) U(alpha, beta) */
void getComplexPairAndPhaseFromUnitary(ComplexMatrix2 u, Complex* alpha, Complex* beta, REAL* globalPhase) {
    
    REAL r0c0Phase = atan2(u.r0c0.imag, u.r0c0.real);
    REAL r1c1Phase = atan2(u.r1c1.imag, u.r1c1.real);
    *globalPhase = (r0c0Phase + r1c1Phase)/2.0;
    
    REAL cosPhase = cos(*globalPhase);
    REAL sinPhase = sin(*globalPhase);
    alpha->real = u.r0c0.real*cosPhase + u.r0c0.imag*sinPhase;
    alpha->imag = u.r0c0.imag*cosPhase - u.r0c0.real*sinPhase;
    beta->real = u.r1c0.real*cosPhase + u.r1c0.imag*sinPhase;
    beta->imag = u.r1c0.imag*cosPhase - u.r1c0.real*sinPhase;
}

void addStringToQASM(QubitRegister qureg, char line[], int lineLen) {
    
    int bufSize = qureg.qasmLog->bufferSize;
    int bufFill = qureg.qasmLog->bufferFill;
    
    // grow QASM buffer if necessary
    if (lineLen + bufFill > bufSize) {
                
        int newBufSize = BUF_GROW_FAC * bufSize;
        if (lineLen + bufFill > newBufSize)
            bufferOverflow();
        
        char* newBuffer = malloc(newBufSize * sizeof *newBuffer);
        sprintf(newBuffer, "%s", qureg.qasmLog->buffer);
        free(qureg.qasmLog->buffer);
        
        qureg.qasmLog->bufferSize = newBufSize;
        qureg.qasmLog->buffer = newBuffer;
    }
        
    // add new str
    sprintf(qureg.qasmLog->buffer + qureg.qasmLog->bufferFill, "%s", line);
    qureg.qasmLog->bufferFill += lineLen;    
}

void addGateToQASM(QubitRegister qureg, TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit, REAL* params, int numParams) {
    
    int len = 0;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    
    // add control labels
    for (int i=0; i < numControlQubits; i++)
        len += snprintf(line+len, MAX_LINE_LEN-len, CTRL_LABEL_PREF);
    
    // add target gate
    len += snprintf(line+len, MAX_LINE_LEN-len, qasmGateLabels[gate]);
    
    // add argument if exists
    if (numParams > 0) {
        len += snprintf(line+len, MAX_LINE_LEN-len, "(");
        for (int i=0; i < numParams; i++) {
            len += snprintf(line+len, MAX_LINE_LEN-len, REAL_STRING_FORMAT, params[i]);
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

void qasm_recordGate(QubitRegister qureg, TargetGate gate, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, NULL, 0, targetQubit, NULL, 0);
}

void qasm_recordParamGate(QubitRegister qureg, TargetGate gate, int targetQubit, REAL param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    REAL params[1] = {param};
    addGateToQASM(qureg, gate, NULL, 0, targetQubit, params, 1);
}

void qasm_recordCompactUnitary(QubitRegister qureg, Complex alpha, Complex beta, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    REAL rz2, ry, rz1;
    getRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    REAL params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, NULL, 0, targetQubit, params, 3);
}

void qasm_recordUnitary(QubitRegister qureg, ComplexMatrix2 u, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    REAL discardedGlobalPhase;
    getComplexPairAndPhaseFromUnitary(u, &alpha, &beta, &discardedGlobalPhase);
    
    REAL rz2, ry, rz1;
    getRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    REAL params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, NULL, 0, targetQubit, params, 3);
}

void qasm_recordControlledGate(QubitRegister qureg, TargetGate gate, int controlQubit, int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int controls[1] = {controlQubit};
    addGateToQASM(qureg, gate, controls, 1, targetQubit, 0, 0);    
}

void qasm_recordControlledParamGate(QubitRegister qureg, TargetGate gate, int controlQubit, int targetQubit, REAL param) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int controls[1] = {controlQubit};
    REAL params[1] = {param};
    addGateToQASM(qureg, gate, controls, 1, targetQubit, params, 1);
}

void qasm_recordControlledCompactUnitary(QubitRegister qureg, Complex alpha, Complex beta, int controlQubit, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    REAL rz2, ry, rz1;
    getRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    int controls[1] = {controlQubit};
    REAL params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controls, 1, targetQubit, params, 3);
}

/** additionally performs Rz on target to restore the global phase lost from u in QASM U(a,b,c) */
void qasm_recordControlledUnitary(QubitRegister qureg, ComplexMatrix2 u, int controlQubit, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    REAL globalPhase;
    getComplexPairAndPhaseFromUnitary(u, &alpha, &beta, &globalPhase);
    
    REAL rz2, ry, rz1;
    getRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    int controls[1] = {controlQubit};
    REAL params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controls, 1, targetQubit, params, 3);
    
    // add Rz
    REAL phaseFix[1] = {globalPhase};
    addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
}

void qasm_recordMultiControlledGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, NULL, 0);
}

void qasm_recordMultiControlledParamGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit, REAL param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    REAL params[1] = {param};
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, params, 1);
}

/** additionally performs Rz on target to restore the global phase lost from u in QASM U(a,b,c) */
void qasm_recordMultiControlledUnitary(QubitRegister qureg, ComplexMatrix2 u, int* controlQubits, const int numControlQubits, const int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    Complex alpha, beta;
    REAL globalPhase;
    getComplexPairAndPhaseFromUnitary(u, &alpha, &beta, &globalPhase);
    
    REAL rz2, ry, rz1;
    getRotAnglesFromComplexPair(alpha, beta, &rz2, &ry, &rz1);
    
    REAL params[3] = {rz2, ry, rz1};
    addGateToQASM(qureg, GATE_UNITARY, controlQubits, numControlQubits, targetQubit, params, 3);
    
    // add Rz
    REAL phaseFix[1] = {globalPhase};
    addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
}

void qasm_recordMeasurement(QubitRegister qureg, const int measureQubit) {

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

void qasm_clearRecorded(QubitRegister qureg) {
    
    // maintains current buffer size
    (qureg.qasmLog->buffer)[0] = '\0';
    qureg.qasmLog->bufferFill = 0;
}

void qasm_free(QubitRegister qureg) {
    
    free(qureg.qasmLog->buffer);
    free(qureg.qasmLog);
}


// will need more for specifying general gates (rotations, unitaries, etc)









