// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Functions for generating QASM output from QuEST circuits
 *
 * @author Tyson Jones
 */

/** TODO
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
# include <stdarg.h>
# include <string.h>

# define QUREG_LABEL "q"        // QASM var-name for the quantum register
# define MESREG_LABEL "c"       // QASM var-name for the classical measurement register
# define CTRL_LABEL_PREF "c"    // QASM syntax which prefixes gates when controlled
# define MEASURE_CMD "measure"  // QASM cmd for measurement operation
# define INIT_ZERO_CMD "reset"  // QASM cmd for setting state 0
# define COMMENT_PREF "//"     // QASM syntax for a comment ;)

# define MAX_LINE_LEN 1024       // maximum length (#chars) of a single QASM instruction
# define BUF_INIT_SIZE 1024     // initial size of the QASM buffer (#chars)
# define BUF_GROW_FAC 2         // growth factor when buffer dynamically resizes
# define MAX_REG_SYMBS 24        // maximum number of single-char symbols in phase function QASM

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
    [GATE_PHASE_SHIFT] = "Rz",// needs phase fix when controlled
    [GATE_SWAP] = "swap",     // needs decomp into cNOTs?
    [GATE_SQRT_SWAP] = "sqrtswap" // needs decomp into cNOTs and Rx(pi/2)?
};

// @TODO make a proper internal error thing
void bufferOverflow(void) {
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

void qasm_recordComment(Qureg qureg, char* comment, ...) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    // write formatted comment to buff
    va_list argp;
    va_start(argp, comment);
    char buff[MAX_LINE_LEN - 4];
    vsnprintf(buff, MAX_LINE_LEN-5, comment, argp);
    va_end(argp);
    
    // add chars to buff, write to QASM logger
    char line[MAX_LINE_LEN + 1]; // for trailing \0
    int len = snprintf(line, MAX_LINE_LEN, "%s %s\n", COMMENT_PREF, buff);
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

void qasm_recordAxisRotation(Qureg qureg, qreal angle, Vector axis, int targetQubit) {
    
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

void qasm_recordMultiControlledGate(Qureg qureg, TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit) {

    if (!qureg.qasmLog->isLogging)
        return;
    
    addGateToQASM(qureg, gate, controlQubits, numControlQubits, targetQubit, NULL, 0);
}

void qasm_recordMultiControlledParamGate(Qureg qureg, TargetGate gate, int* controlQubits, int numControlQubits, int targetQubit, qreal param) {
    
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
void qasm_recordMultiControlledUnitary(Qureg qureg, ComplexMatrix2 u, int* controlQubits, int numControlQubits, int targetQubit) {
    
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
    qasm_recordComment(qureg, "Restoring the discarded global phase of the previous multicontrolled unitary");
    qreal phaseFix[1] = {globalPhase};
    addGateToQASM(qureg, GATE_ROTATE_Z, NULL, 0, targetQubit, phaseFix, 1);
}

void qasm_recordMultiStateControlledUnitary(
    Qureg qureg, ComplexMatrix2 u, int* controlQubits, int* controlState, int numControlQubits, int targetQubit
) {
    if (!qureg.qasmLog->isLogging)
        return;
    
    qasm_recordComment(qureg, "NOTing some gates so that the subsequent unitary is controlled-on-0");
    for (int i=0; i < numControlQubits; i++)
        if (controlState[i] == 0)
            addGateToQASM(qureg, GATE_SIGMA_X, NULL, 0, controlQubits[i], NULL, 0);
    
    qasm_recordMultiControlledUnitary(qureg, u, controlQubits, numControlQubits, targetQubit);

    qasm_recordComment(qureg, "Undoing the NOTing of the controlled-on-0 qubits of the previous unitary");
    for (int i=0; i < numControlQubits; i++)
        if (controlState[i] == 0)
            addGateToQASM(qureg, GATE_SIGMA_X, NULL, 0, controlQubits[i], NULL, 0);
}

void qasm_recordMultiControlledMultiQubitNot(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs) {
    
    if (!qureg.qasmLog->isLogging)
        return;
    
    qasm_recordComment(qureg, "The following %d gates resulted from a single %s() call", numTargs,
        (numCtrls > 0)? "multiControlledMultiQubitNot" : "multiQubitNot");
    
    for (int t=0; t<numTargs; t++)
        addGateToQASM(qureg, GATE_SIGMA_X, ctrls, numCtrls, targs[t], NULL, 0);
}

/* not actually used, D'Oh!
void qasm_recordMultiControlledAxisRotation(Qureg qureg, qreal angle, Vector axis, int* controlQubits, int numControlQubits, int targetQubit) {
    
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

void qasm_recordMeasurement(Qureg qureg, int measureQubit) {

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

void qasm_recordPhaseFunc(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, qreal* overridePhases, int numOverrides) {

    if (!qureg.qasmLog->isLogging)
        return;

    qasm_recordComment(qureg, "Here, applyPhaseFunc() multiplied a complex scalar of the form");

    // record like: 
    //     exp(i (-.5 x^2 + .5 x^(-1.5) - 1.3 x^4 ))
    char line[MAX_LINE_LEN+1];
    int len = snprintf(line, MAX_LINE_LEN, "//     exp(i (");
    for (int t=0; t<numTerms; t++) {
        len += snprintf(line+len, MAX_LINE_LEN-len, 
                (exponents[t] > 0)? 
                    (REAL_QASM_FORMAT " x^" REAL_QASM_FORMAT) : 
                    (REAL_QASM_FORMAT " x^(" REAL_QASM_FORMAT ")"), 
                (t>0)? 
                    absReal(coeffs[t]) : 
                    coeffs[t], 
                exponents[t]);
        if (t < numTerms-1)
            len += snprintf(line+len, MAX_LINE_LEN-len, (coeffs[t+1] > 0)? " + ":" - ");
    }
    len += snprintf(line+len, MAX_LINE_LEN-len, "))\n");

    if (len >= MAX_LINE_LEN)
        bufferOverflow();
    addStringToQASM(qureg, line, len);

    char encBuf[MAX_LINE_LEN];
    if (encoding == UNSIGNED)           sprintf(encBuf, "an unsigned");
    if (encoding == TWOS_COMPLEMENT)    sprintf(encBuf, "a two's complement");
    qasm_recordComment(qureg, "  upon every substate |x>, informed by qubits (under %s binary encoding)", encBuf);

    // record like:
    //     {0, 3, 2}
    len=0;
    len = snprintf(line, MAX_LINE_LEN, "//     {");
    for (int q=0; q<numQubits; q++)
        len += snprintf(line+len, MAX_LINE_LEN-len, (q < numQubits-1)? "%d, ":"%d}\n", qubits[q]);

    if (len >= MAX_LINE_LEN)
        bufferOverflow();
    addStringToQASM(qureg, line, len);

    if (numOverrides > 0) {
        // optionally record like:
        //      |0> -> exp(i .45)
        //      |1> -> exp(i (-.5))
        qasm_recordComment(qureg, "  though with overrides");
        for (int v=0; v<numOverrides; v++)
            qasm_recordComment(qureg, (overridePhases[v] >= 0)? 
                "    |%lld> -> exp(i " REAL_QASM_FORMAT ")" :
                "    |%lld> -> exp(i (" REAL_QASM_FORMAT "))", 
                overrideInds[v], overridePhases[v]);
    }

    // Here, applyPhaseFunction() multiplied a complex scalar of the form
    //      exp(i (.5 x^2 + .5 x^(-1.5) - 1.3 x^4 ))
    // upon every sub-state |x>, informed by qubits
    //      {0, 1, 2}
    // though with overrides
    //      |0> -> exp(i .45)
    //      |1> -> exp(i (-.5))
}

char getPhaseFuncSymbol(int numSymbs, int ind) {

    static char xyz[7] = {'x', 'y', 'z', 't', 'r', 'v', 'u'};
    if (numSymbs <= 7)
        return xyz[ind];

    static char abc[MAX_REG_SYMBS] = {'a','b','c','d','e','f','g','h','j','k','l','m','n','p','q','r','s','t','u','v','w','x','y','z'}; // no i or o
    if (numSymbs <= MAX_REG_SYMBS)
        return abc[ind];

    // we should never reach here, since caller should handle when numSymbs > 24
    bufferOverflow();
    return 'x';
}

void addMultiVarRegsToQASM(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding) {

    char encBuf[MAX_LINE_LEN];
    if (encoding == UNSIGNED)           sprintf(encBuf, "an unsigned");
    if (encoding == TWOS_COMPLEMENT)    sprintf(encBuf, "a two's complement");
    qasm_recordComment(qureg, "  upon substates informed by qubits (under %s binary encoding)", encBuf);

    char line[MAX_LINE_LEN+1];
    int len = 0;

    // record like:
    //     |x> = {0, 3, 2}
    //     |y> = {1, 2}
    int qInd = 0;
    for (int r=0; r<numRegs; r++) {
        len = 0;
        if (numRegs <= MAX_REG_SYMBS)
            len += snprintf(line+len, MAX_LINE_LEN-len, "//     |%c> = {", getPhaseFuncSymbol(numRegs,r));
        else
            len += snprintf(line+len, MAX_LINE_LEN-len, "//     |x%d> = {", r);
        for (int q=0; q<numQubitsPerReg[r]; q++)
            len += snprintf(line+len, MAX_LINE_LEN-len, (q < numQubitsPerReg[r]-1)? "%d, ":"%d}\n", qubits[qInd++]);

        if (len >= MAX_LINE_LEN)
            bufferOverflow();
        addStringToQASM(qureg, line, len);
    }
}

void addMultiVarOverridesToQASM(Qureg qureg, int numRegs, long long int* overrideInds, qreal* overridePhases, int numOverrides) {

    qasm_recordComment(qureg, "  though with overrides");

    char line[MAX_LINE_LEN+1];
    int len = 0;

    // record like:
    //       |x=0, y=1, z=2> -> exp(i .45)
    //       |x=0, y=1, z=5> -> exp(i (-.5))
    int vInd=0;
    for (int v=0; v<numOverrides; v++) {
        len = 0;
        len += snprintf(line+len, MAX_LINE_LEN-len, "//     |");
        for (int r=0; r<numRegs; r++) {
            if (numRegs <= MAX_REG_SYMBS)
                len += snprintf(line+len, MAX_LINE_LEN-len, 
                    (r<numRegs-1)? 
                        "%c=%lld, " : 
                        "%c=%lld>", 
                    getPhaseFuncSymbol(numRegs,r), 
                    overrideInds[vInd++]);
            else
                len += snprintf(line+len, MAX_LINE_LEN-len, 
                    (r<numRegs-1)? 
                        "x%d=%lld, " : 
                        "x%d=%lld>",
                    r,
                    overrideInds[vInd++]);
        }
        len += snprintf(line+len, MAX_LINE_LEN-len, 
            (overridePhases[v] >= 0)? 
                " -> exp(i " REAL_QASM_FORMAT ")\n" : 
                " -> exp(i (" REAL_QASM_FORMAT "))\n", 
                overridePhases[v]);

        if (len >= MAX_LINE_LEN)
            bufferOverflow();
        addStringToQASM(qureg, line, len);
    }
}

void addShiftValuesToQASM(Qureg qureg, enum phaseFunc funcName, int numRegs, qreal* params) {

    char line[MAX_LINE_LEN+1];
    int len = 0;
    int numDeltas;
    if (funcName == SCALED_INVERSE_SHIFTED_NORM)
        numDeltas = numRegs;
    else if (funcName == SCALED_INVERSE_SHIFTED_DISTANCE)
        numDeltas = numRegs / 2;
    else
        return;

    qasm_recordComment(qureg, "  with the additional parameters");

    for (int k=0; k<numDeltas; k++) {
        len = 0;
        len += snprintf(line+len, MAX_LINE_LEN-len, "//     delta%d = " REAL_QASM_FORMAT "\n", k, params[2+k]);
        if (len >= MAX_LINE_LEN)
            bufferOverflow();
        addStringToQASM(qureg, line, len);
    }

}

void qasm_recordMultiVarPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg, long long int* overrideInds, qreal* overridePhases, int numOverrides) {

    if (!qureg.qasmLog->isLogging)
        return;

    qasm_recordComment(qureg, "Here, applyMultiVarPhaseFunc() multiplied a complex scalar of the form");

    // Here, applyMultiVarPhaseFunction() multiplied a complex scalar of the form 
    //     exp(i (
    //         .5 x^2 + .6 x + x 
    //         - y^2 - 5 y - y
    //         + z^2 + z^3 ))

    qasm_recordComment(qureg, "    exp(i (");
    char line[MAX_LINE_LEN+1];
    int len=0;

    int tFlatInd = 0;
    for (int r=0; r<numRegs; r++) {
        len = snprintf(line, MAX_LINE_LEN, "//         ");

        // manually force sign of first term
        len += snprintf(line+len, MAX_LINE_LEN-len, (coeffs[tFlatInd] > 0)? " + ":" - ");

        for (int t=0; t<numTermsPerReg[r]; t++) {
            if (numRegs <= MAX_REG_SYMBS)
                len += snprintf(line+len, MAX_LINE_LEN-len, 
                    (exponents[tFlatInd] > 0)? 
                        REAL_QASM_FORMAT " %c^" REAL_QASM_FORMAT : 
                        REAL_QASM_FORMAT " %c^(" REAL_QASM_FORMAT ")", 
                    absReal(coeffs[tFlatInd]), 
                    getPhaseFuncSymbol(numRegs,r), 
                    exponents[tFlatInd]);
            else
                len += snprintf(line+len, MAX_LINE_LEN-len, 
                    (exponents[tFlatInd] > 0)? 
                        REAL_QASM_FORMAT " x%d^" REAL_QASM_FORMAT : 
                        REAL_QASM_FORMAT " x%d^(" REAL_QASM_FORMAT ")", 
                    absReal(coeffs[tFlatInd]), r, exponents[tFlatInd]);
            if (t < numTermsPerReg[r]-1)
                len += snprintf(line+len, MAX_LINE_LEN-len, (coeffs[tFlatInd+1] > 0)? " + ":" - ");                
            tFlatInd++;
        }

        if (r < numRegs-1)
            len += snprintf(line+len, MAX_LINE_LEN-len, "\n");
        else
            len += snprintf(line+len, MAX_LINE_LEN-len, " ))\n");

        if (len >= MAX_LINE_LEN)
            bufferOverflow();
        addStringToQASM(qureg, line, len);
    }

    addMultiVarRegsToQASM(qureg, qubits, numQubitsPerReg, numRegs, encoding);

    if (numOverrides > 0)
        addMultiVarOverridesToQASM(qureg, numRegs, overrideInds, overridePhases, numOverrides);
}

void qasm_recordNamedPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc funcName, qreal* params, int numParams, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
    // I apologise to my future self and my esteemed readers for such terrible code

    if (!qureg.qasmLog->isLogging)
        return;

    qasm_recordComment(qureg, "Here, applyNamedPhaseFunc() multiplied a complex scalar of form");
    char line[MAX_LINE_LEN+1];

    int len = snprintf(line, MAX_LINE_LEN, "//     exp(i ");

    // record norm-based function, like:  (-1.0) / sqrt(x^2 + y^2 + z^2 + ...)
    if (funcName == NORM || funcName == SCALED_NORM || 
        funcName == INVERSE_NORM || funcName == SCALED_INVERSE_NORM ||
        funcName == SCALED_INVERSE_SHIFTED_NORM)
    {
        // coefficient
        if (funcName == SCALED_NORM || funcName == SCALED_INVERSE_NORM || funcName == SCALED_INVERSE_SHIFTED_NORM)
            len += snprintf(line+len, MAX_LINE_LEN-len, 
                (params[0]>0)? 
                    REAL_QASM_FORMAT " " : 
                    "(" REAL_QASM_FORMAT ") ", 
                params[0]);

        // sqrt(
        if (funcName == NORM || funcName == SCALED_NORM)
            len += snprintf(line+len, MAX_LINE_LEN-len, "sqrt(");
        else if (funcName == INVERSE_NORM)
            len += snprintf(line+len, MAX_LINE_LEN-len, "1 / sqrt(");
        else if (funcName == SCALED_INVERSE_NORM || funcName == SCALED_INVERSE_SHIFTED_NORM)
            len += snprintf(line+len, MAX_LINE_LEN-len, "/ sqrt(");

        // x^2 + y^2 + ...
        if (numRegs <= MAX_REG_SYMBS)
            for (int r=0; r<numRegs; r++) {
                if (funcName == SCALED_INVERSE_SHIFTED_NORM)
                    len += snprintf(line+len, MAX_LINE_LEN-len,
                        (params[2+r] < 0)?
                            "(%c^2+" REAL_QASM_FORMAT ")" :
                            "(%c^2-" REAL_QASM_FORMAT ")",
                        getPhaseFuncSymbol(numRegs,r), fabs(params[2+r]));
                else
                    len += snprintf(line+len, MAX_LINE_LEN-len, "%c^2", getPhaseFuncSymbol(numRegs,r));
                len += snprintf(line+len, MAX_LINE_LEN-len, (r < numRegs - 1)? " + ":"))\n");
            }
        else {
            if (funcName == SCALED_INVERSE_SHIFTED_NORM)
                len += snprintf(line+len, MAX_LINE_LEN-len, "(x0-delta0)^2 + (x1-delta1)^2 + (x2-delta2)^2... ))\n");
            else
                len += snprintf(line+len, MAX_LINE_LEN-len, "x0^2 + x1^2 + x2^2... ))\n");
        }
    }
    // record product-based, like (-1.0) 1/(x y z) ...
    else if (funcName == PRODUCT || funcName == SCALED_PRODUCT || 
        funcName == INVERSE_PRODUCT || funcName == SCALED_INVERSE_PRODUCT)
    {
        // coefficient
        if (funcName == SCALED_PRODUCT || funcName == SCALED_INVERSE_PRODUCT)
            len += snprintf(line+len, MAX_LINE_LEN-len, 
                (params[0]>0)? 
                    REAL_QASM_FORMAT " ":
                    "(" REAL_QASM_FORMAT ") ", 
                params[0]);

        // reciprocal
        if (funcName == INVERSE_PRODUCT)
            len += snprintf(line+len, MAX_LINE_LEN-len, "1 / (");
        else if (funcName == SCALED_INVERSE_PRODUCT)
            len += snprintf(line+len, MAX_LINE_LEN-len, "/ (");

        // x y z ...
        if (numRegs <= MAX_REG_SYMBS)
            for (int r=0; r<numRegs; r++)
                len += snprintf(line+len, MAX_LINE_LEN-len, (r < numRegs - 1)? "%c ":"%c)", getPhaseFuncSymbol(numRegs,r));
        else
            len += snprintf(line+len, MAX_LINE_LEN-len, "x0 x1 x2 ...)");

        // close reciprocal brackets
        if (funcName == INVERSE_PRODUCT || funcName == SCALED_INVERSE_PRODUCT)
            len += snprintf(line+len, MAX_LINE_LEN-len, ")");
        len += snprintf(line+len, MAX_LINE_LEN-len, "\n");
    }
    // record distance-based, like (-1.0) 1/sqrt((x1-x2)^2 + (y1-y2)^2 + ...)
    else if (funcName == DISTANCE || funcName == SCALED_DISTANCE || 
        funcName == INVERSE_DISTANCE || funcName == SCALED_INVERSE_DISTANCE ||
        funcName == SCALED_INVERSE_SHIFTED_DISTANCE)
    {
        // coefficient
        if (funcName == SCALED_DISTANCE || funcName == SCALED_INVERSE_DISTANCE || funcName == SCALED_INVERSE_SHIFTED_DISTANCE)
            len += snprintf(line+len, MAX_LINE_LEN-len, 
                (params[0]>0)? 
                    REAL_QASM_FORMAT " " :
                    "(" REAL_QASM_FORMAT ") ", 
                params[0]);

        // sqrt(
        if (funcName == DISTANCE || funcName == SCALED_DISTANCE)
            len += snprintf(line+len, MAX_LINE_LEN-len, "sqrt(");
        else if (funcName == INVERSE_DISTANCE)
            len += snprintf(line+len, MAX_LINE_LEN-len, "1 / sqrt(");
        else if (funcName == SCALED_INVERSE_DISTANCE || funcName == SCALED_INVERSE_SHIFTED_DISTANCE)
            len += snprintf(line+len, MAX_LINE_LEN-len, "/ sqrt(");

        // (x-y)^2 + (z-t)^2 + ...
        if (numRegs <= MAX_REG_SYMBS)
            for (int r=0; r<numRegs; r+=2) {
                if (funcName == SCALED_INVERSE_SHIFTED_DISTANCE)
                    len += snprintf(line+len, MAX_LINE_LEN-len,
                        (params[2+r/2] < 0)?
                            "(%c-%c+" REAL_QASM_FORMAT ")^2":
                            "(%c-%c-" REAL_QASM_FORMAT ")^2", 
                        getPhaseFuncSymbol(numRegs,r), getPhaseFuncSymbol(numRegs,r+1), fabs(params[2+r/2]));
                else
                    len += snprintf(line+len, MAX_LINE_LEN-len, "(%c-%c)^2", 
                        getPhaseFuncSymbol(numRegs,r), getPhaseFuncSymbol(numRegs,r+1));
                len += snprintf(line+len, MAX_LINE_LEN-len, (r+1 < numRegs-1)? " + ":"))\n");
            }
        else {
            if (funcName == SCALED_INVERSE_SHIFTED_DISTANCE)
                len += snprintf(line+len, MAX_LINE_LEN-len, "(x0-x1-delta0)^2 + (x2-x3-delta1)^2 + ...))\n");
            else
                len += snprintf(line+len, MAX_LINE_LEN-len, "(x0-x1)^2 + (x2-x3)^2 + ...))\n");
        }
    }

    if (len >= MAX_LINE_LEN)
        bufferOverflow();
    addStringToQASM(qureg, line, len);

    addMultiVarRegsToQASM(qureg, qubits, numQubitsPerReg, numRegs, encoding);
    if (numRegs > MAX_REG_SYMBS && (funcName == SCALED_INVERSE_SHIFTED_NORM || funcName == SCALED_INVERSE_SHIFTED_DISTANCE))
        addShiftValuesToQASM(qureg, funcName, numRegs, params);

    if (numOverrides > 0)
        addMultiVarOverridesToQASM(qureg, numRegs, overrideInds,overridePhases, numOverrides);
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
