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
    qasmLog->bufferFill = sprintf(qasmLog->buffer, "qreg q[%d]\n", qureg->numQubitsRepresented);

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


int addQubitIndsString(char buf[], int bufLen, int* controlQubits, const int numControlQubits, const int targetQubit) {

    for (int i=0; i < numControlQubits; i++)
        bufLen += snprintf(buf+bufLen, MAX_LINE_LEN, "%s[%d],", QUREG_LABEL, controlQubits[i]);
    
    bufLen += snprintf(buf+bufLen, MAX_LINE_LEN, "%s[%d]", QUREG_LABEL, targetQubit);

    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
    if (bufLen >= MAX_LINE_LEN)
        bufferOverflow();
    
    return bufLen;
}

int addParamString(char buf[], int bufLen, REAL param) {
    
    // might want to make sci notation or something
    bufLen += sprintf(buf+bufLen, "(");
    bufLen += sprintf(buf+bufLen, REAL_STRING_FORMAT, param);
    bufLen += sprintf(buf+bufLen, ") ");
    return bufLen;
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


void qasm_recordGate(QubitRegister qureg, TargetGate gate, int targetQubit) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    int len;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    len  = sprintf(line, "%s ", qasmGateLabels[gate]);
    len  = addQubitIndsString(line, len, NULL, 0, targetQubit);
    len += sprintf(line+len, ";\n");
    
    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
        
    addStringToQASM(qureg, line, len);
}

void qasm_recordParamGate(QubitRegister qureg, TargetGate gate, int targetQubit, REAL param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    int len;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    len  = sprintf(line, "%s", qasmGateLabels[gate]);
    len  = addParamString(line, len, param);
    len  = addQubitIndsString(line, len, NULL, 0, targetQubit);
    len += sprintf(line+len, ";\n");
    
    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
        
    addStringToQASM(qureg, line, len);
}

void qasm_recordControlledGate(QubitRegister qureg, TargetGate gate, int controlQubit, int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int len;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    len  = sprintf(line, "%s%s ", CTRL_LABEL_PREF, qasmGateLabels[gate]);
    len  = addQubitIndsString(line, len, (int []) {controlQubit}, 1, targetQubit);
    len += sprintf(line+len, ";\n");
    
    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
        
    addStringToQASM(qureg, line, len);
    
}

void qasm_recordControlledParamGate(QubitRegister qureg, TargetGate gate, int controlQubit, int targetQubit, REAL param) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int len;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    len  = sprintf(line, "%s%s", CTRL_LABEL_PREF, qasmGateLabels[gate]);
    len  = addParamString(line, len, param);
    len  = addQubitIndsString(line, len, (int []) {controlQubit}, 1, targetQubit);
    len += sprintf(line+len, ";\n");
    
    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
        
    addStringToQASM(qureg, line, len);
}

void qasm_recordMultiControlledGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    int len;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    len  = sprintf(line, "%s%s ", CTRL_LABEL_PREF, qasmGateLabels[gate]);
    len  = addQubitIndsString(line, len, controlQubits, numControlQubits, targetQubit);
    len += sprintf(line+len, ";\n");
    
    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
    if (len >= MAX_LINE_LEN)
        bufferOverflow();
        
    addStringToQASM(qureg, line, len);
}

void qasm_recordMultiControlledParamGate(QubitRegister qureg, TargetGate gate, int* controlQubits, const int numControlQubits, const int targetQubit, REAL param) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    int len;
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    len  = sprintf(line, "%s%s", CTRL_LABEL_PREF, qasmGateLabels[gate]);
    len  = addParamString(line, len, param);
    len  = addQubitIndsString(line, len, controlQubits, numControlQubits, targetQubit);
    len += sprintf(line+len, ";\n");
    
    // @TODO will snprintf return -1 if failed/overflow? That would ruin our condition
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









