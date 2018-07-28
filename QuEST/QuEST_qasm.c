// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Functions for generating QASM output from QuEST circuits
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_qasm.h"

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# define QUREG_LABEL "q"        // QASM var-name for the quantum register
# define CTRL_LABEL_PREF "c"    // QASM syntax which prefixes gates when controlled
# define MAX_LINE_LEN 100       // maximum length (#chars) of a single QASM instruction
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
    [GATE_ROTATE_Z] = "Rz"
    //[GATE_ROTATE_AROUND_AXIS] = ,
    //[GATE_UNITARY] = ,
    //[GATE_PHASE_SHIFT] = 
};

void qasm_setup(QubitRegister* qureg) {
    
    // populate and attach QASM logger
    QASMLogger *qasmLog = malloc(sizeof qasmLog);
    qasmLog->isLogging = 0;
    qasmLog->bufferSize = BUF_INIT_SIZE;
    qasmLog->buffer = malloc(qasmLog->bufferSize * sizeof *(qasmLog->buffer));
    qasmLog->bufferFill = sprintf(qasmLog->buffer, "qreg q[%d];\n", qureg->numQubitsRepresented);

    qureg->qasmLog = qasmLog;
}

void qasm_startRecording(QubitRegister qureg) {
    qureg.qasmLog->isLogging = 1;
}

void qasm_stopRecording(QubitRegister qureg) {
    qureg.qasmLog->isLogging = 0;
}

// make a proper internal error thing
void bufferOverflow() {
    printf("!!!\nINTERNAL ERROR: QASM line buffer filled!\n!!!");
    exit(1);
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

void addGateToQASM(QubitRegister qureg, TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit, int hasParam, REAL param) {
    
    int len = 0;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    
    // add control labels
    for (int i=0; i < numControlQubits; i++)
        len += snprintf(line+len, MAX_LINE_LEN-len, CTRL_LABEL_PREF);
    
    // add target gate
    len += snprintf(line+len, MAX_LINE_LEN-len, qasmGateLabels[gate]);
    
    // add argument if exists
    if (hasParam) {
        len += snprintf(line+len, MAX_LINE_LEN-len, "(");
        len += snprintf(line+len, MAX_LINE_LEN-len, REAL_STRING_FORMAT, param);
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
    
    addGateToQASM(qureg, gate, NULL, 0, targetQubit, 0, 0);
}

void qasm_recordParamGate(QubitRegister qureg, TargetGate gate, int targetQubit, REAL param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, NULL, 0, targetQubit, 1, param);
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
    addGateToQASM(qureg, gate, controls, 1, targetQubit, 1, param);
}

void qasm_recordMultiControlledGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, 0, 0);
}

void qasm_recordMultiControlledParamGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit, REAL param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, 1, param);

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









