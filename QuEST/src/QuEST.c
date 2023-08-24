// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Implements the QuEST.h API (and some debugging functions) in a hardware-agnostic way, 
 * for both pure and mixed states. These functions mostly wrap hardware-specific functions,
 * and should never call eachother.
 *
 * Density matrices rho of N qubits are flattened to appear as state-vectors |s> of 2N qubits.
 * Operations U rho U^dag are implemented as U^* U |s> and make use of the pure state backend,
 * and often don't need to explicitly compute U^*.
 *
 * @author Tyson Jones (architecture, validation, qasm, density matrices)
 * @author Ania Brown (setDensityAmps())
 * @author Balint Koczor (Kraus maps, calcDensityInnerProduct())
 * @author Nicolas Vogt of HQS (one-qubit damping)
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"
# include "QuEST_validation.h"
# include "QuEST_qasm.h"

# include <stdlib.h>
# include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    
/*
 * state-vector management
 */

Qureg createQureg(int numQubits, QuESTEnv env) {
    validateNumQubitsInQureg(numQubits, env.numRanks, __func__);
    
    Qureg qureg;
    statevec_createQureg(&qureg, numQubits, env);
    qureg.isDensityMatrix = 0;
    qureg.numQubitsRepresented = numQubits;
    qureg.numQubitsInStateVec = numQubits;
    
    qasm_setup(&qureg);
    initZeroState(qureg); // safe call to public function
    return qureg;
}

Qureg createDensityQureg(int numQubits, QuESTEnv env) {
    validateNumQubitsInQureg(2*numQubits, env.numRanks, __func__);
    
    Qureg qureg;
    statevec_createQureg(&qureg, 2*numQubits, env);
    qureg.isDensityMatrix = 1;
    qureg.numQubitsRepresented = numQubits;
    qureg.numQubitsInStateVec = 2*numQubits;
    
    qasm_setup(&qureg);
    initZeroState(qureg); // safe call to public function
    return qureg;
}

Qureg createCloneQureg(Qureg qureg, QuESTEnv env) {

    Qureg newQureg;
    statevec_createQureg(&newQureg, qureg.numQubitsInStateVec, env);
    newQureg.isDensityMatrix = qureg.isDensityMatrix;
    newQureg.numQubitsRepresented = qureg.numQubitsRepresented;
    newQureg.numQubitsInStateVec = qureg.numQubitsInStateVec;
    
    qasm_setup(&newQureg);
    statevec_cloneQureg(newQureg, qureg);
    return newQureg;
}

void destroyQureg(Qureg qureg, QuESTEnv env) {
    statevec_destroyQureg(qureg, env);
    qasm_free(qureg);
}


/*
 * QASM
 */

void startRecordingQASM(Qureg qureg) {
    qasm_startRecording(qureg);
}

void stopRecordingQASM(Qureg qureg) {
    qasm_stopRecording(qureg);
}

void clearRecordedQASM(Qureg qureg) {
    qasm_clearRecorded(qureg);
}

void printRecordedQASM(Qureg qureg) {
    qasm_printRecorded(qureg);
}

void writeRecordedQASMToFile(Qureg qureg, char* filename) {
    int success = qasm_writeRecordedToFile(qureg, filename);
    validateFileOpened(success, filename, __func__);
}


/*
 * state initialisation
 */

void initZeroState(Qureg qureg) {
    statevec_initZeroState(qureg); // valid for both statevec and density matrices
    
    qasm_recordInitZero(qureg);
}

void initBlankState(Qureg qureg) {
    statevec_initBlankState(qureg);
    
    qasm_recordComment(qureg, "Here, the register was initialised to an unphysical all-zero-amplitudes 'state'.");
}

void initPlusState(Qureg qureg) {
    if (qureg.isDensityMatrix)
        densmatr_initPlusState(qureg);
    else
        statevec_initPlusState(qureg);
    
    qasm_recordInitPlus(qureg);
}

void initClassicalState(Qureg qureg, long long int stateInd) {
    validateStateIndex(qureg, stateInd, __func__);
    
    if (qureg.isDensityMatrix)
        densmatr_initClassicalState(qureg, stateInd);
    else
        statevec_initClassicalState(qureg, stateInd);
    
    qasm_recordInitClassical(qureg, stateInd);
}

void initPureState(Qureg qureg, Qureg pure) {
    validateSecondQuregStateVec(pure, __func__);
    validateMatchingQuregDims(qureg, pure, __func__);

    if (qureg.isDensityMatrix)
        densmatr_initPureState(qureg, pure);
    else
        statevec_cloneQureg(qureg, pure);
    
    qasm_recordComment(qureg, "Here, the register was initialised to an undisclosed given pure state.");
}

void initStateFromAmps(Qureg qureg, qreal* reals, qreal* imags) {
    
    statevec_setAmps(qureg, 0, reals, imags, qureg.numAmpsTotal);
    
    qasm_recordComment(qureg, "Here, the register was initialised to an undisclosed given pure state.");
}

void cloneQureg(Qureg targetQureg, Qureg copyQureg) {
    validateMatchingQuregTypes(targetQureg, copyQureg, __func__);
    validateMatchingQuregDims(targetQureg, copyQureg, __func__);
    
    statevec_cloneQureg(targetQureg, copyQureg);
}

void setQuregToPauliHamil(Qureg qureg, PauliHamil hamil) {
    validateDensityMatrQureg(qureg, __func__);
    validatePauliHamil(hamil, __func__);
    validateMatchingQuregPauliHamilDims(qureg, hamil, __func__);

    densmatr_setQuregToPauliHamil(qureg, hamil);
}


/*
 * unitary gates
 */

void hadamard(Qureg qureg, int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_hadamard(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_hadamard(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_HADAMARD, targetQubit);
}

void rotateX(Qureg qureg, int targetQubit, qreal angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateX(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateX(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_X, targetQubit, angle);
}

void rotateY(Qureg qureg, int targetQubit, qreal angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateY(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateY(qureg, targetQubit+qureg.numQubitsRepresented, angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_Y, targetQubit, angle);
}

void rotateZ(Qureg qureg, int targetQubit, qreal angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateZ(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateZ(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_Z, targetQubit, angle);
}

void controlledRotateX(Qureg qureg, int controlQubit, int targetQubit, qreal angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateX(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateX(qureg, controlQubit+shift, targetQubit+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_ROTATE_X, controlQubit, targetQubit, angle);
}

void controlledRotateY(Qureg qureg, int controlQubit, int targetQubit, qreal angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateY(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateY(qureg, controlQubit+shift, targetQubit+shift, angle); // rotateY is real
    }

    qasm_recordControlledParamGate(qureg, GATE_ROTATE_Y, controlQubit, targetQubit, angle);
}

void controlledRotateZ(Qureg qureg, int controlQubit, int targetQubit, qreal angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateZ(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateZ(qureg, controlQubit+shift, targetQubit+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_ROTATE_Z, controlQubit, targetQubit, angle);
}

void twoQubitUnitary(Qureg qureg, int targetQubit1, int targetQubit2, ComplexMatrix4 u) {
    validateMultiTargets(qureg, (int []) {targetQubit1, targetQubit2}, 2, __func__);
    validateTwoQubitUnitaryMatrix(qureg, u, __func__);
    
    statevec_twoQubitUnitary(qureg, targetQubit1, targetQubit2, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_twoQubitUnitary(qureg, targetQubit1+shift, targetQubit2+shift, getConjugateMatrix4(u));
    }
    
    qasm_recordComment(qureg, "Here, an undisclosed 2-qubit unitary was applied.");
}

void controlledTwoQubitUnitary(Qureg qureg, int controlQubit, int targetQubit1, int targetQubit2, ComplexMatrix4 u) {
    validateMultiControlsMultiTargets(qureg, (int[]) {controlQubit}, 1, (int[]) {targetQubit1, targetQubit2}, 2, __func__);
    validateTwoQubitUnitaryMatrix(qureg, u, __func__);
    
    statevec_controlledTwoQubitUnitary(qureg, controlQubit, targetQubit1, targetQubit2, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledTwoQubitUnitary(qureg, controlQubit+shift, targetQubit1+shift, targetQubit2+shift, getConjugateMatrix4(u));
    }

    qasm_recordComment(qureg, "Here, an undisclosed controlled 2-qubit unitary was applied.");
}

void multiControlledTwoQubitUnitary(Qureg qureg, int* controlQubits, int numControlQubits, int targetQubit1, int targetQubit2, ComplexMatrix4 u) {
    validateMultiControlsMultiTargets(qureg, controlQubits, numControlQubits, (int[]) {targetQubit1, targetQubit2}, 2, __func__);
    validateTwoQubitUnitaryMatrix(qureg, u, __func__);
    
    long long int ctrlQubitsMask = getQubitBitMask(controlQubits, numControlQubits);
    statevec_multiControlledTwoQubitUnitary(qureg, ctrlQubitsMask, targetQubit1, targetQubit2, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiControlledTwoQubitUnitary(qureg, ctrlQubitsMask<<shift, targetQubit1+shift, targetQubit2+shift, getConjugateMatrix4(u));
    }
    
    qasm_recordComment(qureg, "Here, an undisclosed multi-controlled 2-qubit unitary was applied.");
}

void multiQubitUnitary(Qureg qureg, int* targs, int numTargs, ComplexMatrixN u) {
    validateMultiTargets(qureg, targs, numTargs, __func__);
    validateMultiQubitUnitaryMatrix(qureg, u, numTargs, __func__);
    
    statevec_multiQubitUnitary(qureg, targs, numTargs, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(u);
        statevec_multiQubitUnitary(qureg, targs, numTargs, u);
        shiftIndices(targs, numTargs, -shift);
        setConjugateMatrixN(u);
    }
    
    qasm_recordComment(qureg, "Here, an undisclosed multi-qubit unitary was applied.");
}

void controlledMultiQubitUnitary(Qureg qureg, int ctrl, int* targs, int numTargs, ComplexMatrixN u) {
    validateMultiControlsMultiTargets(qureg, (int[]) {ctrl}, 1, targs, numTargs, __func__);
    validateMultiQubitUnitaryMatrix(qureg, u, numTargs, __func__);
    
    statevec_controlledMultiQubitUnitary(qureg, ctrl, targs, numTargs, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(u);
        statevec_controlledMultiQubitUnitary(qureg, ctrl+shift, targs, numTargs, u);
        shiftIndices(targs, numTargs, -shift);
        setConjugateMatrixN(u);
    }
    
    qasm_recordComment(qureg, "Here, an undisclosed controlled multi-qubit unitary was applied.");
}

void multiControlledMultiQubitUnitary(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, ComplexMatrixN u) {
    validateMultiControlsMultiTargets(qureg, ctrls, numCtrls, targs, numTargs, __func__);
    validateMultiQubitUnitaryMatrix(qureg, u, numTargs, __func__);
    
    long long int ctrlMask = getQubitBitMask(ctrls, numCtrls);
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targs, numTargs, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(u);
        statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask<<shift, targs, numTargs, u);
        shiftIndices(targs, numTargs, -shift);
        setConjugateMatrixN(u);
    }
    
    qasm_recordComment(qureg, "Here, an undisclosed multi-controlled multi-qubit unitary was applied.");
}

void unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u) {
    validateTarget(qureg, targetQubit, __func__);
    validateOneQubitUnitaryMatrix(u, __func__);
    
    statevec_unitary(qureg, targetQubit, u);
    if (qureg.isDensityMatrix) {
        statevec_unitary(qureg, targetQubit+qureg.numQubitsRepresented, getConjugateMatrix2(u));
    }
    
    qasm_recordUnitary(qureg, u, targetQubit);
}

void controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateOneQubitUnitaryMatrix(u, __func__);
    
    statevec_controlledUnitary(qureg, controlQubit, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledUnitary(qureg, controlQubit+shift, targetQubit+shift, getConjugateMatrix2(u));
    }
    
    qasm_recordControlledUnitary(qureg, u, controlQubit, targetQubit);
}

void multiControlledUnitary(Qureg qureg, int* controlQubits, int numControlQubits, int targetQubit, ComplexMatrix2 u) {
    validateMultiControlsTarget(qureg, controlQubits, numControlQubits, targetQubit, __func__);
    validateOneQubitUnitaryMatrix(u, __func__);
    
    long long int ctrlQubitsMask = getQubitBitMask(controlQubits, numControlQubits);
    long long int ctrlFlipMask = 0;
    statevec_multiControlledUnitary(qureg, ctrlQubitsMask, ctrlFlipMask, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiControlledUnitary(qureg, ctrlQubitsMask<<shift, ctrlFlipMask<<shift, targetQubit+shift, getConjugateMatrix2(u));
    }
    
    qasm_recordMultiControlledUnitary(qureg, u, controlQubits, numControlQubits, targetQubit);
}

void multiStateControlledUnitary(Qureg qureg, int* controlQubits, int* controlState, int numControlQubits, int targetQubit, ComplexMatrix2 u) {
    validateMultiControlsTarget(qureg, controlQubits, numControlQubits, targetQubit, __func__);
    validateOneQubitUnitaryMatrix(u, __func__);
    validateControlState(controlState, numControlQubits, __func__);

    long long int ctrlQubitsMask = getQubitBitMask(controlQubits, numControlQubits);
    long long int ctrlFlipMask = getControlFlipMask(controlQubits, controlState, numControlQubits);
    statevec_multiControlledUnitary(qureg, ctrlQubitsMask, ctrlFlipMask, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiControlledUnitary(qureg, ctrlQubitsMask<<shift, ctrlFlipMask<<shift, targetQubit+shift, getConjugateMatrix2(u));
    }
    
    qasm_recordMultiStateControlledUnitary(qureg, u, controlQubits, controlState, numControlQubits, targetQubit);
}

void compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta) {
    validateTarget(qureg, targetQubit, __func__);
    validateUnitaryComplexPair(alpha, beta, __func__);
    
    statevec_compactUnitary(qureg, targetQubit, alpha, beta);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_compactUnitary(qureg, targetQubit+shift, getConjugateScalar(alpha), getConjugateScalar(beta));
    }

    qasm_recordCompactUnitary(qureg, alpha, beta, targetQubit);
}

void controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateUnitaryComplexPair(alpha, beta, __func__);
    
    statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledCompactUnitary(qureg, 
            controlQubit+shift, targetQubit+shift, 
            getConjugateScalar(alpha), getConjugateScalar(beta));
    }
    
    qasm_recordControlledCompactUnitary(qureg, alpha, beta, controlQubit, targetQubit);
}

void pauliX(Qureg qureg, int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_pauliX(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_pauliX(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_X, targetQubit);
}

void pauliY(Qureg qureg, int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_pauliY(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_pauliYConj(qureg, targetQubit + qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_Y, targetQubit);
}

void pauliZ(Qureg qureg, int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_pauliZ(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_pauliZ(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_Z, targetQubit);
}

void sGate(Qureg qureg, int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sGate(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_S, targetQubit);
}

void tGate(Qureg qureg, int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_tGate(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_tGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_T, targetQubit);
}

void phaseShift(Qureg qureg, int targetQubit, qreal angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_phaseShift(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_phaseShift(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_PHASE_SHIFT, targetQubit, angle);
}

void controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle) {
    validateControlTarget(qureg, idQubit1, idQubit2, __func__);
    
    statevec_controlledPhaseShift(qureg, idQubit1, idQubit2, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPhaseShift(qureg, idQubit1+shift, idQubit2+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_PHASE_SHIFT, idQubit1, idQubit2, angle);
}

void multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle) {
    validateMultiQubits(qureg, controlQubits, numControlQubits, __func__);
    
    statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, -angle);
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledParamGate(qureg, GATE_PHASE_SHIFT, controlQubits, numControlQubits-1, controlQubits[numControlQubits-1], angle);
}

void controlledNot(Qureg qureg, int controlQubit, int targetQubit) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledNot(qureg, controlQubit, targetQubit);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledNot(qureg, controlQubit+shift, targetQubit+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_X, controlQubit, targetQubit);
}

void multiQubitNot(Qureg qureg, int* targs, int numTargs) {
    validateMultiTargets(qureg, targs, numTargs, __func__);
    
    long long int targMask = getQubitBitMask(targs, numTargs);
    statevec_multiControlledMultiQubitNot(qureg, 0, targMask);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiControlledMultiQubitNot(qureg, 0, targMask<<shift);
    }
    
    qasm_recordMultiControlledMultiQubitNot(qureg, NULL, 0, targs, numTargs);
}

void multiControlledMultiQubitNot(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs) {
    validateMultiControlsMultiTargets(qureg, ctrls, numCtrls, targs, numTargs, __func__);
    
    long long int ctrlMask = getQubitBitMask(ctrls, numCtrls);
    long long int targMask = getQubitBitMask(targs, numTargs);
    statevec_multiControlledMultiQubitNot(qureg, ctrlMask, targMask);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiControlledMultiQubitNot(qureg, ctrlMask<<shift, targMask<<shift);
    }
    
    qasm_recordMultiControlledMultiQubitNot(qureg, ctrls, numCtrls, targs, numTargs);
}

void controlledPauliY(Qureg qureg, int controlQubit, int targetQubit) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledPauliY(qureg, controlQubit, targetQubit);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPauliYConj(qureg, controlQubit+shift, targetQubit+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_Y, controlQubit, targetQubit);
}

void controlledPhaseFlip(Qureg qureg, int idQubit1, int idQubit2) {
    validateControlTarget(qureg, idQubit1, idQubit2, __func__);
    
    statevec_controlledPhaseFlip(qureg, idQubit1, idQubit2);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPhaseFlip(qureg, idQubit1+shift, idQubit2+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_Z, idQubit1, idQubit2);
}

void multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits) {
    validateMultiQubits(qureg, controlQubits, numControlQubits, __func__);
    
    statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledGate(qureg, GATE_SIGMA_Z, controlQubits, numControlQubits-1, controlQubits[numControlQubits-1]);
}

void rotateAroundAxis(Qureg qureg, int rotQubit, qreal angle, Vector axis) {
    validateTarget(qureg, rotQubit, __func__);
    validateVector(axis, __func__);
    
    statevec_rotateAroundAxis(qureg, rotQubit, angle, axis);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_rotateAroundAxisConj(qureg, rotQubit+shift, angle, axis);
    }
    
    qasm_recordAxisRotation(qureg, angle, axis, rotQubit);
}

void controlledRotateAroundAxis(Qureg qureg, int controlQubit, int targetQubit, qreal angle, Vector axis) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateVector(axis, __func__);
    
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateAroundAxisConj(qureg, controlQubit+shift, targetQubit+shift, angle, axis);
    }
    
    qasm_recordControlledAxisRotation(qureg, angle, axis, controlQubit, targetQubit);
}

void swapGate(Qureg qureg, int qb1, int qb2) {
    validateUniqueTargets(qureg, qb1, qb2, __func__);

    statevec_swapQubitAmps(qureg, qb1, qb2);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_swapQubitAmps(qureg, qb1+shift, qb2+shift);
    }

    qasm_recordControlledGate(qureg, GATE_SWAP, qb1, qb2);
}

void sqrtSwapGate(Qureg qureg, int qb1, int qb2) {
    validateUniqueTargets(qureg, qb1, qb2, __func__);
    validateMultiQubitMatrixFitsInNode(qureg, 2, __func__); // uses 2qb unitary in QuEST_common

    statevec_sqrtSwapGate(qureg, qb1, qb2);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_sqrtSwapGateConj(qureg, qb1+shift, qb2+shift);
    }

    qasm_recordControlledGate(qureg, GATE_SQRT_SWAP, qb1, qb2);
}

void multiRotateZ(Qureg qureg, int* qubits, int numQubits, qreal angle) {
    validateMultiTargets(qureg, qubits, numQubits, __func__);
    
    long long int mask = getQubitBitMask(qubits, numQubits);
    statevec_multiRotateZ(qureg, mask, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiRotateZ(qureg, mask << shift, -angle);
    }
    
    // @TODO: create actual QASM
    qasm_recordComment(qureg, 
        "Here a %d-qubit multiRotateZ of angle " REAL_QASM_FORMAT " was performed (QASM not yet implemented)",
        numQubits, angle);
}

void multiControlledMultiRotateZ(Qureg qureg, int* controlQubits, int numControls, int* targetQubits, int numTargets, qreal angle) {
    validateMultiControlsMultiTargets(qureg, controlQubits, numControls, targetQubits, numTargets, __func__);
    
    long long int ctrlMask = getQubitBitMask(controlQubits, numControls);
    long long int targMask = getQubitBitMask(targetQubits, numTargets);
    statevec_multiControlledMultiRotateZ(qureg, ctrlMask, targMask, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_multiControlledMultiRotateZ(qureg, ctrlMask<<shift, targMask<<shift, - angle);
    }
    
    // @TODO: create actual QASM
    qasm_recordComment(qureg, 
        "Here a %d-control %d-target multiControlledMultiRotateZ of angle " REAL_QASM_FORMAT " was performed (QASM not yet implemented)",
        numControls, numTargets, angle);
}

void multiRotatePauli(Qureg qureg, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle) {
    validateMultiTargets(qureg, targetQubits, numTargets, __func__);
    validatePauliCodes(targetPaulis, numTargets, __func__);
    
    int conj=0;
    statevec_multiRotatePauli(qureg, targetQubits, targetPaulis, numTargets, angle, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targetQubits, numTargets, shift);
        statevec_multiRotatePauli(qureg, targetQubits, targetPaulis, numTargets, angle, conj);
        shiftIndices(targetQubits, numTargets, -shift);
    }
    
    // @TODO: create actual QASM
    qasm_recordComment(qureg, 
        "Here a %d-qubit multiRotatePauli of angle " REAL_QASM_FORMAT " was performed (QASM not yet implemented)",
        numTargets, angle);
}

void multiControlledMultiRotatePauli(Qureg qureg, int* controlQubits, int numControls, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle) {
    validateMultiControlsMultiTargets(qureg, controlQubits, numControls, targetQubits, numTargets, __func__);
    validatePauliCodes(targetPaulis, numTargets, __func__);
    
    int conj=0;
    long long int ctrlMask = getQubitBitMask(controlQubits, numControls);
    statevec_multiControlledMultiRotatePauli(qureg, ctrlMask, targetQubits, targetPaulis, numTargets, angle, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targetQubits, numTargets, shift);
        statevec_multiControlledMultiRotatePauli(qureg, ctrlMask<<shift, targetQubits, targetPaulis, numTargets, angle, conj);
        shiftIndices(targetQubits, numTargets, -shift);
    }
    
    // @TODO: create actual QASM
    qasm_recordComment(qureg, 
        "Here a %d-control %d-target multiControlledMultiRotatePauli of angle " REAL_QASM_FORMAT " was performed (QASM not yet implemented)",
        numControls, numTargets, angle);
}

void applyPhaseFunc(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms) {
    validateMultiQubits(qureg, qubits, numQubits, __func__);
    validateBitEncoding(numQubits, encoding, __func__);
    validatePhaseFuncTerms(numQubits, encoding, coeffs, exponents, numTerms, NULL, 0, __func__);

    int conj = 0;
    statevec_applyPhaseFuncOverrides(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, NULL, NULL, 0, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftIndices(qubits, numQubits, qureg.numQubitsRepresented);
        statevec_applyPhaseFuncOverrides(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, NULL, NULL, 0, conj);
        shiftIndices(qubits, numQubits, - qureg.numQubitsRepresented);
    }

    qasm_recordPhaseFunc(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, NULL, NULL, 0);
}

void applyPhaseFuncOverrides(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
    validateMultiQubits(qureg, qubits, numQubits, __func__);
    validateBitEncoding(numQubits, encoding, __func__);
    validatePhaseFuncOverrides(numQubits, encoding, overrideInds, numOverrides, __func__);
    validatePhaseFuncTerms(numQubits, encoding, coeffs, exponents, numTerms, overrideInds, numOverrides, __func__);

    int conj = 0;
    statevec_applyPhaseFuncOverrides(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, overrideInds, overridePhases, numOverrides, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftIndices(qubits, numQubits, qureg.numQubitsRepresented);
        statevec_applyPhaseFuncOverrides(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, overrideInds, overridePhases, numOverrides, conj);
        shiftIndices(qubits, numQubits, - qureg.numQubitsRepresented);
    }

    qasm_recordPhaseFunc(qureg, qubits, numQubits, encoding, coeffs, exponents, numTerms, overrideInds, overridePhases, numOverrides);
}

void applyMultiVarPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg) {
    validateQubitSubregs(qureg, qubits, numQubitsPerReg, numRegs, __func__);
    validateMultiRegBitEncoding(numQubitsPerReg, numRegs, encoding, __func__);
    validateMultiVarPhaseFuncTerms(numQubitsPerReg, numRegs, encoding, exponents, numTermsPerReg, __func__);

    int conj = 0;
    statevec_applyMultiVarPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, NULL, NULL, 0, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, qureg.numQubitsRepresented);
        statevec_applyMultiVarPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, NULL, NULL, 0, conj);
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, - qureg.numQubitsRepresented);
    }

    qasm_recordMultiVarPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, NULL, NULL, 0);
}

void applyMultiVarPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
    validateQubitSubregs(qureg, qubits, numQubitsPerReg, numRegs, __func__);
    validateMultiRegBitEncoding(numQubitsPerReg, numRegs, encoding, __func__);
    validateMultiVarPhaseFuncTerms(numQubitsPerReg, numRegs, encoding, exponents, numTermsPerReg, __func__);
    validateMultiVarPhaseFuncOverrides(numQubitsPerReg, numRegs, encoding, overrideInds, numOverrides, __func__);

    int conj = 0;
    statevec_applyMultiVarPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases, numOverrides, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, qureg.numQubitsRepresented);
        statevec_applyMultiVarPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases, numOverrides, conj);
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, - qureg.numQubitsRepresented);
    }

    qasm_recordMultiVarPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, coeffs, exponents, numTermsPerReg, overrideInds, overridePhases, numOverrides);
}

void applyNamedPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode) {
    validateQubitSubregs(qureg, qubits, numQubitsPerReg, numRegs, __func__);
    validateMultiRegBitEncoding(numQubitsPerReg, numRegs, encoding, __func__);
    validatePhaseFuncName(functionNameCode, numRegs, 0, __func__);

    int conj = 0;
    statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, NULL, 0, NULL, NULL, 0, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, qureg.numQubitsRepresented);
        statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, NULL, 0, NULL, NULL, 0, conj);
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, - qureg.numQubitsRepresented);
    }

    qasm_recordNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, NULL, 0, NULL, NULL, 0);
}

void applyNamedPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
    validateQubitSubregs(qureg, qubits, numQubitsPerReg, numRegs, __func__);
    validateMultiRegBitEncoding(numQubitsPerReg, numRegs, encoding, __func__);
    validatePhaseFuncName(functionNameCode, numRegs, 0, __func__);
    validateMultiVarPhaseFuncOverrides(numQubitsPerReg, numRegs, encoding, overrideInds, numOverrides, __func__);

    int conj = 0;
    statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, NULL, 0, overrideInds, overridePhases, numOverrides, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, qureg.numQubitsRepresented);
        statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, NULL, 0, overrideInds, overridePhases, numOverrides, conj);
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, - qureg.numQubitsRepresented);
    }

    qasm_recordNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, NULL, 0, overrideInds, overridePhases, numOverrides);
}

void applyParamNamedPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams) {
    validateQubitSubregs(qureg, qubits, numQubitsPerReg, numRegs, __func__);
    validateMultiRegBitEncoding(numQubitsPerReg, numRegs, encoding, __func__);
    validatePhaseFuncName(functionNameCode, numRegs, numParams, __func__);

    int conj = 0;
    statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, NULL, NULL, 0, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, qureg.numQubitsRepresented);
        statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, NULL, NULL, 0, conj);
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, - qureg.numQubitsRepresented);
    }

    qasm_recordNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, NULL, NULL, 0);
}

void applyParamNamedPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams, long long int* overrideInds, qreal* overridePhases, int numOverrides) {
    validateQubitSubregs(qureg, qubits, numQubitsPerReg, numRegs, __func__);
    validateMultiRegBitEncoding(numQubitsPerReg, numRegs, encoding, __func__);
    validatePhaseFuncName(functionNameCode, numRegs, numParams, __func__);
    validateMultiVarPhaseFuncOverrides(numQubitsPerReg, numRegs, encoding, overrideInds, numOverrides, __func__);

    int conj = 0;
    statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, overrideInds, overridePhases, numOverrides, conj);
    if (qureg.isDensityMatrix) {
        conj = 1;
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, qureg.numQubitsRepresented);
        statevec_applyParamNamedPhaseFuncOverrides(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, overrideInds, overridePhases, numOverrides, conj);
        shiftSubregIndices(qubits, numQubitsPerReg, numRegs, - qureg.numQubitsRepresented);
    }

    qasm_recordNamedPhaseFunc(qureg, qubits, numQubitsPerReg, numRegs, encoding, functionNameCode, params, numParams, overrideInds, overridePhases, numOverrides);
}

void applyQFT(Qureg qureg, int* qubits, int numQubits) {
    validateMultiTargets(qureg, qubits, numQubits, __func__);
    
    qasm_recordComment(qureg, "Beginning of QFT circuit");
        
    agnostic_applyQFT(qureg, qubits, numQubits);

    qasm_recordComment(qureg, "End of QFT circuit");
}

void applyFullQFT(Qureg qureg) {

    qasm_recordComment(qureg, "Beginning of QFT circuit");
        
    int qubits[100];
    for (int i=0; i<qureg.numQubitsRepresented; i++)
        qubits[i] = i;
    agnostic_applyQFT(qureg, qubits, qureg.numQubitsRepresented);

    qasm_recordComment(qureg, "End of QFT circuit");
}

void applyProjector(Qureg qureg, int qubit, int outcome) {
    validateTarget(qureg, qubit, __func__);
    validateOutcome(outcome, __func__);
     
    qreal renorm = 1;
    
    if (qureg.isDensityMatrix)
        densmatr_collapseToKnownProbOutcome(qureg, qubit, outcome, renorm);
    else
        statevec_collapseToKnownProbOutcome(qureg, qubit, outcome, renorm);
    
    qasm_recordComment(qureg, "Here, qubit %d was un-physically projected into outcome %d", qubit, outcome);
}

void diagonalUnitary(Qureg qureg, int* qubits, int numQubits, SubDiagonalOp op) {
    validateSubDiagOpTargets(qureg, qubits, numQubits, op, __func__);
    validateUnitarySubDiagOp(op, __func__);
    
    int shift = 0;
    int conj = 0;
    statevec_applySubDiagonalOp(qureg, qubits, op, conj);
    
    if (qureg.isDensityMatrix) {
        conj = 1;
        shift = qureg.numQubitsRepresented;
        
        shiftIndices(qubits, numQubits, shift);
        statevec_applySubDiagonalOp(qureg, qubits, op, conj);
        shiftIndices(qubits, numQubits, - shift);
    }
    
    qasm_recordComment(qureg, "Here, the register was modified by an undisclosed diagonal unitary (via diagonalUnitary).");
}



/*
 * register attributes
 */

int getNumQubits(Qureg qureg) {
    return qureg.numQubitsRepresented;
}

long long int getNumAmps(Qureg qureg) {
    validateStateVecQureg(qureg, __func__);
    
    return qureg.numAmpsTotal;
}

qreal getRealAmp(Qureg qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateAmpIndex(qureg, index, __func__);
    
    return statevec_getRealAmp(qureg, index);
}

qreal getImagAmp(Qureg qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateAmpIndex(qureg, index, __func__);
    
    return statevec_getImagAmp(qureg, index);
}

qreal getProbAmp(Qureg qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateAmpIndex(qureg, index, __func__);
    
    return statevec_getProbAmp(qureg, index);
}

Complex getAmp(Qureg qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateAmpIndex(qureg, index, __func__);
    
    Complex amp;
    amp.real = statevec_getRealAmp(qureg, index);
    amp.imag = statevec_getImagAmp(qureg, index);
    return amp;
}

Complex getDensityAmp(Qureg qureg, long long int row, long long int col) {
    validateDensityMatrQureg(qureg, __func__);
    validateAmpIndex(qureg, row, __func__);
    validateAmpIndex(qureg, col, __func__);
    
    long long ind = row + col*(1LL << qureg.numQubitsRepresented);
    Complex amp;
    amp.real = statevec_getRealAmp(qureg, ind);
    amp.imag = statevec_getImagAmp(qureg, ind);
    return amp;
}


/*
 * non-unitary actions
 */

qreal collapseToOutcome(Qureg qureg, int measureQubit, int outcome) {
    validateTarget(qureg, measureQubit, __func__);
    validateOutcome(outcome, __func__);
    
    qreal outcomeProb;
    if (qureg.isDensityMatrix) {
        outcomeProb = densmatr_calcProbOfOutcome(qureg, measureQubit, outcome);
        validateMeasurementProb(outcomeProb, __func__);
        densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
    } else {
        outcomeProb = statevec_calcProbOfOutcome(qureg, measureQubit, outcome);
        validateMeasurementProb(outcomeProb, __func__);
        statevec_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
    }
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcomeProb;
}

int measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb) {
    validateTarget(qureg, measureQubit, __func__);

    int outcome;
    if (qureg.isDensityMatrix)
        outcome = densmatr_measureWithStats(qureg, measureQubit, outcomeProb);
    else
        outcome = statevec_measureWithStats(qureg, measureQubit, outcomeProb);
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcome;
}

int measure(Qureg qureg, int measureQubit) {
    validateTarget(qureg, measureQubit, __func__);
    
    int outcome;
    qreal discardedProb;
    if (qureg.isDensityMatrix)
        outcome = densmatr_measureWithStats(qureg, measureQubit, &discardedProb);
    else
        outcome = statevec_measureWithStats(qureg, measureQubit, &discardedProb);
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcome;
}

void mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg) {
    validateDensityMatrQureg(combineQureg, __func__);
    validateDensityMatrQureg(otherQureg, __func__);
    validateMatchingQuregDims(combineQureg, otherQureg, __func__);
    validateProb(otherProb, __func__);
    
    densmatr_mixDensityMatrix(combineQureg, otherProb, otherQureg);
}

void setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps) {
    validateStateVecQureg(qureg, __func__);
    validateNumAmps(qureg, startInd, numAmps, __func__);
    
    statevec_setAmps(qureg, startInd, reals, imags, numAmps);
    
    qasm_recordComment(qureg, "Here, some amplitudes in the statevector were manually edited.");
}

void setDensityAmps(Qureg qureg, long long int startRow, long long int startCol, qreal* reals, qreal* imags, long long int numAmps) {
    validateDensityMatrQureg(qureg, __func__);
    validateNumDensityAmps(qureg, startRow, startCol, numAmps, __func__);

    long long int startInd = startRow + startCol*(1 << qureg.numQubitsRepresented);
    statevec_setAmps(qureg, startInd, reals, imags, numAmps);
    
    qasm_recordComment(qureg, "Here, some amplitudes in the density matrix were manually edited.");
}

void setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out) {
    validateMatchingQuregTypes(qureg1, qureg2, __func__);
    validateMatchingQuregTypes(qureg1, out, __func__);
    validateMatchingQuregDims(qureg1, qureg2,  __func__);
    validateMatchingQuregDims(qureg1, out, __func__);

    statevec_setWeightedQureg(fac1, qureg1, fac2, qureg2, facOut, out);

    qasm_recordComment(out, "Here, the register was modified to an undisclosed and possibly unphysical state (setWeightedQureg).");
} 

void applyPauliSum(Qureg inQureg, enum pauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg outQureg) {
    validateMatchingQuregTypes(inQureg, outQureg, __func__);
    validateMatchingQuregDims(inQureg, outQureg, __func__);
    validateNumPauliSumTerms(numSumTerms, __func__);
    validatePauliCodes(allPauliCodes, numSumTerms*inQureg.numQubitsRepresented, __func__);
    
    statevec_applyPauliSum(inQureg, allPauliCodes, termCoeffs, numSumTerms, outQureg);
    
    qasm_recordComment(outQureg, "Here, the register was modified to an undisclosed and possibly unphysical state (applyPauliSum).");
}

void applyPauliHamil(Qureg inQureg, PauliHamil hamil, Qureg outQureg) {
    validateMatchingQuregTypes(inQureg, outQureg, __func__);
    validateMatchingQuregDims(inQureg, outQureg, __func__);
    validatePauliHamil(hamil, __func__);
    validateMatchingQuregPauliHamilDims(inQureg, hamil, __func__);
    
    statevec_applyPauliSum(inQureg, hamil.pauliCodes, hamil.termCoeffs, hamil.numSumTerms, outQureg);
    
    qasm_recordComment(outQureg, "Here, the register was modified to an undisclosed and possibly unphysical state (applyPauliHamil).");
}

void applyTrotterCircuit(Qureg qureg, PauliHamil hamil, qreal time, int order, int reps) {
    validateTrotterParams(order, reps, __func__);
    validatePauliHamil(hamil, __func__);
    validateMatchingQuregPauliHamilDims(qureg, hamil, __func__);
    
    qasm_recordComment(qureg, 
        "Beginning of Trotter circuit (time " REAL_QASM_FORMAT ", order %d, %d repetitions).",
        time, order, reps);
        
    agnostic_applyTrotterCircuit(qureg, hamil, time, order, reps);

    qasm_recordComment(qureg, "End of Trotter circuit");
}

void applyMatrix2(Qureg qureg, int targetQubit, ComplexMatrix2 u) {
    validateTarget(qureg, targetQubit, __func__);
    
    // actually just left-multiplies any complex matrix
    statevec_unitary(qureg, targetQubit, u);

    qasm_recordComment(qureg, "Here, an undisclosed 2-by-2 matrix (possibly non-unitary) was multiplied onto qubit %d", targetQubit);
}

void applyMatrix4(Qureg qureg, int targetQubit1, int targetQubit2, ComplexMatrix4 u) {
    validateMultiTargets(qureg, (int []) {targetQubit1, targetQubit2}, 2, __func__);
    validateMultiQubitMatrixFitsInNode(qureg, 2, __func__);
    
    // actually just left-multiplies any complex matrix
    statevec_twoQubitUnitary(qureg, targetQubit1, targetQubit2, u);

    qasm_recordComment(qureg, "Here, an undisclosed 4-by-4 matrix (possibly non-unitary) was multiplied onto qubits %d and %d", targetQubit1, targetQubit2);
}

void applyMatrixN(Qureg qureg, int* targs, int numTargs, ComplexMatrixN u) {
    validateMultiTargets(qureg, targs, numTargs, __func__);
    validateMultiQubitMatrix(qureg, u, numTargs, __func__);
    
    // actually just left-multiplies any complex matrix
    statevec_multiQubitUnitary(qureg, targs, numTargs, u);
    
    int dim = (1 << numTargs);
    qasm_recordComment(qureg, "Here, an undisclosed %d-by-%d matrix (possibly non-unitary) was multiplied onto %d undisclosed qubits", dim, dim, numTargs);
}

void applyGateMatrixN(Qureg qureg, int* targs, int numTargs, ComplexMatrixN u) {
    validateMultiTargets(qureg, targs, numTargs, __func__);
    validateMultiQubitMatrix(qureg, u, numTargs, __func__);
    
    statevec_multiQubitUnitary(qureg, targs, numTargs, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(u);
        statevec_multiQubitUnitary(qureg, targs, numTargs, u);
        shiftIndices(targs, numTargs, -shift);
        setConjugateMatrixN(u);
    }
    
    int dim = (1 << numTargs);
    qasm_recordComment(qureg, "Here, an undisclosed %d-by-%d gate matrix (possibly non-unitary) was applied to %d undisclosed qubits", dim, dim, numTargs);
}

void applyMultiControlledGateMatrixN(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, ComplexMatrixN m) {
    validateMultiControlsMultiTargets(qureg, ctrls, numCtrls, targs, numTargs, __func__);
    validateMultiQubitMatrix(qureg, m, numTargs, __func__);
    
    long long int ctrlMask = getQubitBitMask(ctrls, numCtrls);
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targs, numTargs, m);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(targs, numTargs, shift);
        setConjugateMatrixN(m);
        statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask<<shift, targs, numTargs, m);
        shiftIndices(targs, numTargs, -shift);
        setConjugateMatrixN(m);
    }
    
    int dim = (1 << numTargs);
    qasm_recordComment(qureg, "Here, an undisclosed %d-controlled %d-by-%d gate matrix (possibly non-unitary) was applied to %d undisclosed qubits", numCtrls, dim, dim, numTargs);
}

void applyMultiControlledMatrixN(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, ComplexMatrixN u) {
    validateMultiControlsMultiTargets(qureg, ctrls, numCtrls, targs, numTargs, __func__);
    validateMultiQubitMatrix(qureg, u, numTargs, __func__);
    
    // actually just left-multiplies any complex matrix
    long long int ctrlMask = getQubitBitMask(ctrls, numCtrls);
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targs, numTargs, u);
    
    int numTot = numTargs + numCtrls;
    int dim = (1 << numTot );
    qasm_recordComment(qureg, "Here, an undisclosed %d-by-%d matrix (possibly non-unitary, and including %d controlled qubits) was multiplied onto %d undisclosed qubits", dim, dim, numCtrls, numTot);
}

void applyDiagonalOp(Qureg qureg, DiagonalOp op) {
    validateDiagonalOp(qureg, op, __func__);

    if (qureg.isDensityMatrix)
        densmatr_applyDiagonalOp(qureg, op);
    else
        statevec_applyDiagonalOp(qureg, op);

    qasm_recordComment(qureg, "Here, the register was modified to an undisclosed and possibly unphysical state (via applyDiagonalOp).");
}

void applySubDiagonalOp(Qureg qureg, int* qubits, int numQubits, SubDiagonalOp op) {
    validateSubDiagOpTargets(qureg, qubits, numQubits, op, __func__);
    
    int conj = 0;
    statevec_applySubDiagonalOp(qureg, qubits, op, conj);
    
    qasm_recordComment(qureg, "Here, the register was modified to an undisclosed and possibly unphysical state (via applySubDiagonalOp).");
}

void applyGateSubDiagonalOp(Qureg qureg, int* qubits, int numQubits, SubDiagonalOp op) {
    validateSubDiagOpTargets(qureg, qubits, numQubits, op, __func__);
    
    int shift = 0;
    int conj = 0;
    statevec_applySubDiagonalOp(qureg, qubits, op, conj);
    
    if (qureg.isDensityMatrix) {
        conj = 1;
        shift = qureg.numQubitsRepresented;
        
        shiftIndices(qubits, numQubits, shift);
        statevec_applySubDiagonalOp(qureg, qubits, op, conj);
        shiftIndices(qubits, numQubits, - shift);
    }
    
    qasm_recordComment(qureg, "Here, the register was modified by an undisclosed sub-diagonal unitary, though which did not enforce numerical unitarity.");
}


/*
 * calculations
 */

qreal calcTotalProb(Qureg qureg) {
    if (qureg.isDensityMatrix)  
            return densmatr_calcTotalProb(qureg);
        else
            return statevec_calcTotalProb(qureg);
}

Complex calcInnerProduct(Qureg bra, Qureg ket) {
    validateStateVecQureg(bra, __func__);
    validateStateVecQureg(ket, __func__);
    validateMatchingQuregDims(bra, ket,  __func__);
    
    return statevec_calcInnerProduct(bra, ket);
}

qreal calcDensityInnerProduct(Qureg rho1, Qureg rho2) {
    validateDensityMatrQureg(rho1, __func__);
    validateDensityMatrQureg(rho2, __func__);
    validateMatchingQuregDims(rho1, rho2, __func__);
    
    return densmatr_calcInnerProduct(rho1, rho2);
}

qreal calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome) {
    validateTarget(qureg, measureQubit, __func__);
    validateOutcome(outcome, __func__);
    
    if (qureg.isDensityMatrix)
        return densmatr_calcProbOfOutcome(qureg, measureQubit, outcome);
    else
        return statevec_calcProbOfOutcome(qureg, measureQubit, outcome);
}

void calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits) {
    validateMultiTargets(qureg, qubits, numQubits, __func__);

    if (qureg.isDensityMatrix)
        densmatr_calcProbOfAllOutcomes(retProbs, qureg, qubits, numQubits);
    else
        statevec_calcProbOfAllOutcomes(retProbs, qureg, qubits, numQubits);
}

qreal calcPurity(Qureg qureg) {
    validateDensityMatrQureg(qureg, __func__);
    
    return densmatr_calcPurity(qureg);
}

qreal calcFidelity(Qureg qureg, Qureg pureState) {
    validateSecondQuregStateVec(pureState, __func__);
    validateMatchingQuregDims(qureg, pureState, __func__);
    
    if (qureg.isDensityMatrix)
        return densmatr_calcFidelity(qureg, pureState);
    else
        return statevec_calcFidelity(qureg, pureState);
}

qreal calcExpecPauliProd(Qureg qureg, int* targetQubits, enum pauliOpType* pauliCodes, int numTargets, Qureg workspace) {
    validateMultiTargets(qureg, targetQubits, numTargets, __func__);
    validatePauliCodes(pauliCodes, numTargets, __func__);
    validateMatchingQuregTypes(qureg, workspace, __func__);
    validateMatchingQuregDims(qureg, workspace, __func__);
    
    return statevec_calcExpecPauliProd(qureg, targetQubits, pauliCodes, numTargets, workspace);
}

qreal calcExpecPauliSum(Qureg qureg, enum pauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg workspace) {
    validateNumPauliSumTerms(numSumTerms, __func__);
    validatePauliCodes(allPauliCodes, numSumTerms*qureg.numQubitsRepresented, __func__);
    validateMatchingQuregTypes(qureg, workspace, __func__);
    validateMatchingQuregDims(qureg, workspace, __func__);
    
    return statevec_calcExpecPauliSum(qureg, allPauliCodes, termCoeffs, numSumTerms, workspace);
}

qreal calcExpecPauliHamil(Qureg qureg, PauliHamil hamil, Qureg workspace) {
    validateMatchingQuregTypes(qureg, workspace, __func__);
    validateMatchingQuregDims(qureg, workspace, __func__);
    validatePauliHamil(hamil, __func__);
    validateMatchingQuregPauliHamilDims(qureg, hamil, __func__);
    
    return statevec_calcExpecPauliSum(qureg, hamil.pauliCodes, hamil.termCoeffs, hamil.numSumTerms, workspace);
}

Complex calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {
    validateDiagonalOp(qureg, op, __func__);
    
    if (qureg.isDensityMatrix)
        return densmatr_calcExpecDiagonalOp(qureg, op);
    else
        return statevec_calcExpecDiagonalOp(qureg, op);
}

qreal calcHilbertSchmidtDistance(Qureg a, Qureg b) {
    validateDensityMatrQureg(a, __func__);
    validateDensityMatrQureg(b, __func__);
    validateMatchingQuregDims(a, b, __func__);
    
    return densmatr_calcHilbertSchmidtDistance(a, b);
}


/*
 * decoherence
 */

void mixDephasing(Qureg qureg, int targetQubit, qreal prob) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateOneQubitDephaseProb(prob, __func__);
    
    densmatr_mixDephasing(qureg, targetQubit, 2*prob);
    qasm_recordComment(qureg, 
        "Here, a phase (Z) error occured on qubit %d with probability " REAL_QASM_FORMAT, targetQubit, prob);
}

void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    validateDensityMatrQureg(qureg, __func__);
    validateUniqueTargets(qureg, qubit1, qubit2, __func__);
    validateTwoQubitDephaseProb(prob, __func__);

    ensureIndsIncrease(&qubit1, &qubit2);
    densmatr_mixTwoQubitDephasing(qureg, qubit1, qubit2, (4*prob)/3.0);
    qasm_recordComment(qureg,
        "Here, a phase (Z) error occured on either or both of qubits "
        "%d and %d with total probability " REAL_QASM_FORMAT, qubit1, qubit2, prob);
}

void mixDepolarising(Qureg qureg, int targetQubit, qreal prob) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateOneQubitDepolProb(prob, __func__);
    
    densmatr_mixDepolarising(qureg, targetQubit, (4*prob)/3.0);
    qasm_recordComment(qureg,
        "Here, a homogeneous depolarising error (X, Y, or Z) occured on "
        "qubit %d with total probability " REAL_QASM_FORMAT, targetQubit, prob);
}

void mixDamping(Qureg qureg, int targetQubit, qreal prob) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateOneQubitDampingProb(prob, __func__);
    
    densmatr_mixDamping(qureg, targetQubit, prob);
}

void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    validateDensityMatrQureg(qureg, __func__);
    validateUniqueTargets(qureg, qubit1, qubit2, __func__);
    validateTwoQubitDepolProb(prob, __func__);
    
    ensureIndsIncrease(&qubit1, &qubit2);
    densmatr_mixTwoQubitDepolarising(qureg, qubit1, qubit2, (16*prob)/15.0);
    qasm_recordComment(qureg,
        "Here, a homogeneous depolarising error occured on qubits %d and %d "
        "with total probability " REAL_QASM_FORMAT, qubit1, qubit2, prob);
}

void mixPauli(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, qubit, __func__);
    validateOneQubitPauliProbs(probX, probY, probZ, __func__);
    
    densmatr_mixPauli(qureg, qubit, probX, probY, probZ);
    qasm_recordComment(qureg,
        "Here, X, Y and Z errors occured on qubit %d with probabilities "
        REAL_QASM_FORMAT ", " REAL_QASM_FORMAT " and " REAL_QASM_FORMAT " respectively", qubit, probX, probY, probZ);
}

void mixKrausMap(Qureg qureg, int target, ComplexMatrix2 *ops, int numOps) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, target, __func__);
    validateOneQubitKrausMap(qureg, ops, numOps, __func__);
    
    densmatr_mixKrausMap(qureg, target, ops, numOps);
    qasm_recordComment(qureg, 
        "Here, an undisclosed Kraus map was effected on qubit %d", target);
}

void mixTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps) {
    validateDensityMatrQureg(qureg, __func__);
    validateMultiTargets(qureg, (int[]) {target1,target2}, 2, __func__);
    validateTwoQubitKrausMap(qureg, ops, numOps, __func__);
    
    densmatr_mixTwoQubitKrausMap(qureg, target1, target2, ops, numOps);
    qasm_recordComment(qureg, 
        "Here, an undisclosed two-qubit Kraus map was effected on qubits %d and %d", target1, target2);
}

void mixMultiQubitKrausMap(Qureg qureg, int* targets, int numTargets, ComplexMatrixN* ops, int numOps) {
    validateDensityMatrQureg(qureg, __func__);
    validateMultiTargets(qureg, targets, numTargets, __func__);
    validateMultiQubitKrausMap(qureg, numTargets, ops, numOps, __func__);
    
    densmatr_mixMultiQubitKrausMap(qureg, targets, numTargets, ops, numOps);
    qasm_recordComment(qureg,
        "Here, an undisclosed %d-qubit Kraus map was applied to undisclosed qubits", numTargets);
}

void mixNonTPKrausMap(Qureg qureg, int target, ComplexMatrix2 *ops, int numOps) {
    validateDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, target, __func__);
    validateOneQubitKrausMapDimensions(qureg, ops, numOps, __func__);
    
    densmatr_mixKrausMap(qureg, target, ops, numOps);
    qasm_recordComment(qureg, 
        "Here, an undisclosed non-trace-preserving Kraus map was effected on qubit %d", target);
}

void mixNonTPTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps) {
    validateDensityMatrQureg(qureg, __func__);
    validateMultiTargets(qureg, (int[]) {target1,target2}, 2, __func__);
    validateTwoQubitKrausMapDimensions(qureg, ops, numOps, __func__);
    
    densmatr_mixTwoQubitKrausMap(qureg, target1, target2, ops, numOps);
    qasm_recordComment(qureg, 
        "Here, an undisclosed non-trace-preserving two-qubit Kraus map was effected on qubits %d and %d", target1, target2);
}

void mixNonTPMultiQubitKrausMap(Qureg qureg, int* targets, int numTargets, ComplexMatrixN* ops, int numOps) {
    validateDensityMatrQureg(qureg, __func__);
    validateMultiTargets(qureg, targets, numTargets, __func__);
    validateMultiQubitKrausMapDimensions(qureg, numTargets, ops, numOps, __func__);
    
    densmatr_mixMultiQubitKrausMap(qureg, targets, numTargets, ops, numOps);
    qasm_recordComment(qureg,
        "Here, an undisclosed non-trace-preserving %d-qubit Kraus map was applied to undisclosed qubits", numTargets);
}

/*
 * other data structures 
 */
 
 ComplexMatrixN createComplexMatrixN(int numQubits) {
     validateNumQubitsInMatrix(numQubits, __func__);

     int numRows = 1 << numQubits;

     ComplexMatrixN m = {
         .numQubits = numQubits,
         .real = malloc(numRows * sizeof *m.real),
         .imag = malloc(numRows * sizeof *m.imag)};

     for (int n=0; n < 1<<numQubits; n++) {
         m.real[n] = calloc(numRows, sizeof **m.real);
         m.imag[n] = calloc(numRows, sizeof **m.imag);
     }
     
     // error if the ComplexMatrixN was not successfully malloc'ds
     validateMatrixInit(m, __func__);
     
     return m;
 }
 
void destroyComplexMatrixN(ComplexMatrixN m) {
    /* this checks m.real/imag != NULL, which is only ever set when the mallocs 
     * in createComplexMatrixN fail, which already prompts an error. Hence 
     * this check if useless 
     */
    validateMatrixInit(m, __func__);
    
    int numRows = 1 << m.numQubits;
    for (int r=0; r < numRows; r++) {
        free(m.real[r]);
        free(m.imag[r]);
    }
    free(m.real);
    free(m.imag);
}

#ifndef _WIN32
void initComplexMatrixN(ComplexMatrixN m, qreal re[][1<<m.numQubits], qreal im[][1<<m.numQubits]) {
    validateMatrixInit(m, __func__);
    
    int dim = 1 << m.numQubits;
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++) {
            m.real[i][j] = re[i][j];
            m.imag[i][j] = im[i][j];
        }
}
#endif

PauliHamil createPauliHamil(int numQubits, int numSumTerms) {
    validateHamilParams(numQubits, numSumTerms, __func__);
    
    PauliHamil h;
    h.numQubits = numQubits;
    h.numSumTerms = numSumTerms;
    h.termCoeffs = malloc(numSumTerms * sizeof *h.termCoeffs);
    h.pauliCodes = malloc(numQubits*numSumTerms * sizeof *h.pauliCodes);
    
    // initialise pauli codes to identity 
    for (int i=0; i<numQubits*numSumTerms; i++)
        h.pauliCodes[i] = PAULI_I;
    
    return h;
}

void destroyPauliHamil(PauliHamil h) {
    
    free(h.termCoeffs);
    free(h.pauliCodes);
}

PauliHamil createPauliHamilFromFile(char* fn) {
    
    /* The validation in this function must close the file handle and free 
     * allocated memory before raising an error (whether that's a C exit, or
     * an overriden C++ exception).
     */
    
	FILE* file = fopen(fn, "r");
    int success = (file != NULL);
    validateFileOpened(success, fn, __func__);
    
    /* file format: coeff {term} \n where {term} is #numQubits values of 
	 * 0 1 2 3 signifying I X Y Z acting on that qubit index
	 */
	
	// count the number of qubits (ignore trailing whitespace)
	int numQubits = -1; // to exclude coeff at start
	char ch = getc(file);
    char prev = '0'; // anything not space
	while (ch != '\n' && ch != EOF) {
		if (ch == ' ' && prev != ' ') // skip multiple spaces
			numQubits++;
        prev = ch;
        ch = getc(file);
    }
    // edge-case: if we hit EOF/newline without a space
    if (prev != ' ')
        numQubits++;

    /* TODO:
     * The below code may break on Windows where newlines are multiple characters
     */
	
	// count the number of terms (being cautious of trailing newlines)
    int numTerms = 0;
    prev = '\n';
    rewind(file);
	while ((ch=getc(file)) != EOF) {
		if (ch == '\n' && prev != '\n')
			numTerms++;
        prev = ch;
    }
    // edge-case: if we hit EOF without a newline, count that line
    if (prev != '\n')
        numTerms++;

    // validate the inferred number of terms and qubits (closes file if error)
    validateHamilFileParams(numQubits, numTerms, file, fn, __func__);
    
    // allocate space for PauliHamil data
    PauliHamil h = createPauliHamil(numQubits, numTerms);
     
    // specifier for a qreal number then a space 
    char strSpec[50];
    strcpy(strSpec, REAL_SPECIFIER);
    strcat(strSpec, " ");
    
    // collect coefficients and terms
	rewind(file);
	for (int t=0; t<numTerms; t++) {
		
		// record coefficient, and validate (closes file and frees h if error)
        success = fscanf(file, strSpec, &(h.termCoeffs[t])) == 1;
        validateHamilFileCoeffParsed(success, h, file, fn, __func__);

		// record Pauli operations, and validate (closes file and frees h if error)
        for (int q=0; q<numQubits; q++) {
            int i = t*numQubits + q;
            
            // verbose, to avoid type warnings
            int code;
            success = fscanf(file, "%d ", &code) == 1;
            h.pauliCodes[i] = (enum pauliOpType) code;
            validateHamilFilePauliParsed(success, h, file, fn, __func__);
            validateHamilFilePauliCode(h.pauliCodes[i], h, file, fn, __func__);
		}

		// the trailing newline is magically eaten
	}

    fclose(file);
	return h;
}

void initPauliHamil(PauliHamil hamil, qreal* coeffs, enum pauliOpType* codes) {
    validateHamilParams(hamil.numQubits, hamil.numSumTerms, __func__);
    validatePauliCodes(codes, hamil.numSumTerms*hamil.numQubits, __func__);
    
    int i=0;
    for (int t=0; t<hamil.numSumTerms; t++) {
        hamil.termCoeffs[t] = coeffs[t];
        for (int q=0; q<hamil.numQubits; q++) {
            hamil.pauliCodes[i] = codes[i];
            i++;
        }
    }
}

DiagonalOp createDiagonalOp(int numQubits, QuESTEnv env) {
    validateNumQubitsInDiagOp(numQubits, env.numRanks, __func__);
    
    return agnostic_createDiagonalOp(numQubits, env);
}

void destroyDiagonalOp(DiagonalOp op, QuESTEnv env) {
    // env accepted for API consistency
    validateDiagOpInit(op, __func__);
    
    agnostic_destroyDiagonalOp(op);
}

void syncDiagonalOp(DiagonalOp op) {
    validateDiagOpInit(op, __func__);
    
    agnostic_syncDiagonalOp(op);
}

void initDiagonalOp(DiagonalOp op, qreal* real, qreal* imag) {
    validateDiagOpInit(op, __func__);
    
    agnostic_setDiagonalOpElems(op, 0, real, imag, 1LL << op.numQubits);
}

void setDiagonalOpElems(DiagonalOp op, long long int startInd, qreal* real, qreal* imag, long long int numElems) {
    validateDiagOpInit(op, __func__);
    validateNumElems(op, startInd, numElems, __func__);
    
    agnostic_setDiagonalOpElems(op, startInd, real, imag, numElems);
}

void initDiagonalOpFromPauliHamil(DiagonalOp op, PauliHamil hamil) {
    validateHamilParams(hamil.numQubits, hamil.numSumTerms, __func__);
    validateDiagOpInit(op, __func__);
    validateDiagPauliHamil(op, hamil, __func__);
    
    agnostic_initDiagonalOpFromPauliHamil(op, hamil);
}

DiagonalOp createDiagonalOpFromPauliHamilFile(char* fn, QuESTEnv env) {
    PauliHamil h = createPauliHamilFromFile(fn); // validates fn
    validateDiagPauliHamilFromFile(h, env.numRanks, __func__);  // destroys h if invalid

    DiagonalOp op = agnostic_createDiagonalOp(h.numQubits, env);
    agnostic_initDiagonalOpFromPauliHamil(op, h);
    
    destroyPauliHamil(h);
    return op;
}

SubDiagonalOp createSubDiagonalOp(int numQubits) {
    validateNumQubitsInSubDiagOp(numQubits, __func__);

    SubDiagonalOp op;
    op.numQubits = numQubits;
    op.numElems = (1LL << numQubits);
    
    // allocate CPU memory (initialised to zero)
    op.real = (qreal*) calloc(op.numElems, sizeof(qreal));
    op.imag = (qreal*) calloc(op.numElems, sizeof(qreal));
    
    return op;
}

void destroySubDiagonalOp(SubDiagonalOp op) {
    free(op.real);
    free(op.imag);
}

/*
 * debug
 */

void initDebugState(Qureg qureg) {
    statevec_initDebugState(qureg);
}

void reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank)  {
    statevec_reportStateToScreen(qureg, env, reportRank);
}

void reportPauliHamil(PauliHamil hamil) {
    validatePauliHamil(hamil, __func__);
    
    for (int t=0; t<hamil.numSumTerms; t++) {
        printf(REAL_QASM_FORMAT, hamil.termCoeffs[t]);
        printf("\t");
        for (int q=0; q<hamil.numQubits; q++)
            printf("%d ", (int) hamil.pauliCodes[q+t*hamil.numQubits]);
        printf("\n");
    }
}

int  getQuEST_PREC(void) {
  return sizeof(qreal)/4;
}

void seedQuESTDefault(QuESTEnv* env) {
    
    // seed Mersenne Twister random number generator with two keys -- time and pid
    unsigned long int keys[2];
    getQuESTDefaultSeedKey(keys);
    seedQuEST(env, keys, 2);
}

void getQuESTSeeds(QuESTEnv env, unsigned long int** seeds, int* numSeeds) {
    *seeds = env.seeds;
    *numSeeds = env.numSeeds;
}

void copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps) {
    validateNumAmps(qureg, startInd, numAmps, __func__);
    
    statevec_copySubstateToGPU(qureg, startInd, numAmps);
}

void copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps) {
    validateNumAmps(qureg, startInd, numAmps, __func__);
    
    statevec_copySubstateFromGPU(qureg, startInd, numAmps);
}

  

#ifdef __cplusplus
}
#endif
