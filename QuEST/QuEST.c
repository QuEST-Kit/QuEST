// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Implements the QuEST.h API (and some debugging functions) in a hardware-agnostic way, 
 * for both pure and mixed states. These functions mostly wrap hardware-specific functions,
 * and should never call eachother.
 *
 * Density matrices rho of N qubits are flattened to appear as state-vectors |s> of 2N qubits.
 * Operations U rho U^dag are implemented as U^* U |s> and make use of the pure state backend,
 * and often don't need to explicitly compute U^*.
 */

// @TODO unit test the density functionality of all below methods
// @TODO for initPureState:
//      - densmatr_initPureStateDistributed on CPU

// @TODO add controlled Hadamard

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"
# include "QuEST_validation.h"
# include "QuEST_ops.h"
# include "QuEST_qasm.h"

#ifdef __cplusplus
extern "C" {
#endif

void createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
    validateCreateNumQubits(numQubits, __func__);
    
    statevec_createQubitRegister(qureg, numQubits, env);
    qureg->isDensityMatrix = 0;
    qureg->numQubitsRepresented = numQubits;
    qureg->numQubitsInStateVec = numQubits;
    
    qasm_setup(qureg);
}

void createDensityQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
    validateCreateNumQubits(numQubits, __func__);
    
    statevec_createQubitRegister(qureg, 2*numQubits, env);
    qureg->isDensityMatrix = 1;
    qureg->numQubitsRepresented = numQubits;
    qureg->numQubitsInStateVec = 2*numQubits;
    
    qasm_setup(qureg);
}

void destroyQubitRegister(QubitRegister qureg, QuESTEnv env) {
    statevec_destroyQubitRegister(qureg, env);
    qasm_free(qureg);
}

void startRecordingQASM(QubitRegister qureg) {
    qasm_startRecording(qureg);
}

void stopRecordingQASM(QubitRegister qureg) {
    qasm_stopRecording(qureg);
}

void clearRecordedQASM(QubitRegister qureg) {
    qasm_clearRecorded(qureg);
}

void printRecordedQASM(QubitRegister qureg) {
    qasm_printRecorded(qureg);
}

void writeRecordedQASMToFile(QubitRegister qureg, char* filename) {
    int success = qasm_writeRecordedToFile(qureg, filename);
    validateFileOpened(success, __func__);
}

void initStateZero(QubitRegister qureg) {
    statevec_initStateZero(qureg); // valid for both statevec and density matrices
    
    // @TODO QASM?
}

void initStatePlus(QubitRegister qureg) {
    if (qureg.isDensityMatrix)
        densmatr_initStatePlus(qureg);
    else
        statevec_initStatePlus(qureg);
    
    // @TODO QASM?
}

void initClassicalState(QubitRegister qureg, long long int stateInd) {
    validateStateIndex(qureg, stateInd, __func__);
    
    if (qureg.isDensityMatrix)
        densmatr_initClassicalState(qureg, stateInd);
    else
        statevec_initClassicalState(qureg, stateInd);
    
    // @TODO QASM?
}

void hadamard(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_hadamard(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_hadamard(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_HADAMARD, targetQubit);
}

void rotateX(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateX(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateX(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_X, targetQubit, angle);
}

void rotateY(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateY(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateY(qureg, targetQubit+qureg.numQubitsRepresented, angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_Y, targetQubit, angle);
}

void rotateZ(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_rotateZ(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_rotateZ(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_ROTATE_Z, targetQubit, angle);
}

void controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateX(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateX(qureg, controlQubit+shift, targetQubit+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_ROTATE_X, controlQubit, targetQubit, angle);
}

void controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateY(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateY(qureg, controlQubit+shift, targetQubit+shift, angle);
    }

    qasm_recordControlledParamGate(qureg, GATE_ROTATE_Y, controlQubit, targetQubit, angle);
}

void controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledRotateZ(qureg, controlQubit, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateZ(qureg, controlQubit+shift, targetQubit+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_ROTATE_Z, controlQubit, targetQubit, angle);
}

void unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u) {
    validateTarget(qureg, targetQubit, __func__);
    validateUnitaryMatrix(u, __func__);
    
    statevec_unitary(qureg, targetQubit, u);
    if (qureg.isDensityMatrix) {
        statevec_unitary(qureg, targetQubit+qureg.numQubitsRepresented, getConjugateMatrix(u));
    }
    
    qasm_recordUnitary(qureg, u, targetQubit);
}

void controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateUnitaryMatrix(u, __func__);
    
    statevec_controlledUnitary(qureg, controlQubit, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledUnitary(qureg, controlQubit+shift, targetQubit+shift, getConjugateMatrix(u));
    }
    
    qasm_recordControlledUnitary(qureg, u, controlQubit, targetQubit);
}

void multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) {
    validateMultiControlsTarget(qureg, controlQubits, numControlQubits, targetQubit, __func__);
    validateUnitaryMatrix(u, __func__);
    
    statevec_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit, u);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit+shift, getConjugateMatrix(u));
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledUnitary(qureg, u, controlQubits, numControlQubits, targetQubit);
}

void compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) {
    validateTarget(qureg, targetQubit, __func__);
    validateUnitaryComplexPair(alpha, beta, __func__);
    
    statevec_compactUnitary(qureg, targetQubit, alpha, beta);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_compactUnitary(qureg, targetQubit+shift, getConjugateScalar(alpha), getConjugateScalar(beta));
    }

    qasm_recordCompactUnitary(qureg, alpha, beta, targetQubit);
}

void controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) {
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

void sigmaX(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sigmaX(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sigmaX(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_X, targetQubit);
}

void sigmaY(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sigmaY(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sigmaYConj(qureg, targetQubit + qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_Y, targetQubit);
}

void sigmaZ(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sigmaZ(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sigmaZ(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_SIGMA_Z, targetQubit);
}

void sGate(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_sGate(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_sGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_S, targetQubit);
}

void tGate(QubitRegister qureg, const int targetQubit) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_tGate(qureg, targetQubit);
    if (qureg.isDensityMatrix) {
        statevec_tGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
    }
    
    qasm_recordGate(qureg, GATE_T, targetQubit);
}

void phaseShift(QubitRegister qureg, const int targetQubit, REAL angle) {
    validateTarget(qureg, targetQubit, __func__);
    
    statevec_phaseShift(qureg, targetQubit, angle);
    if (qureg.isDensityMatrix) {
        statevec_phaseShift(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
    }
    
    qasm_recordParamGate(qureg, GATE_PHASE_SHIFT, targetQubit, angle);
}

void controlledPhaseShift(QubitRegister qureg, const int idQubit1, const int idQubit2, REAL angle) {
    validateControlTarget(qureg, idQubit1, idQubit2, __func__); // a little bit semantically dodgy
    
    statevec_controlledPhaseShift(qureg, idQubit1, idQubit2, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPhaseShift(qureg, idQubit1+shift, idQubit2+shift, -angle);
    }
    
    qasm_recordControlledParamGate(qureg, GATE_PHASE_SHIFT, idQubit1, idQubit2, angle);
}

void multiControlledPhaseShift(QubitRegister qureg, int *controlQubits, int numControlQubits, REAL angle) {
    validateMultiControls(qureg, controlQubits, numControlQubits, __func__);
    
    statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledParamGate(qureg, GATE_PHASE_SHIFT, controlQubits, numControlQubits-1, controlQubits[numControlQubits-1], angle);
}

void controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledNot(qureg, controlQubit, targetQubit);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledNot(qureg, controlQubit+shift, targetQubit+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_X, controlQubit, targetQubit);
}

void controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    
    statevec_controlledSigmaY(qureg, controlQubit, targetQubit);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledSigmaYConj(qureg, controlQubit+shift, targetQubit+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_Y, controlQubit, targetQubit);
}

void controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2) {
    validateControlTarget(qureg, idQubit1, idQubit2, __func__); // a little bit semantically dodgy
    
    statevec_controlledPhaseFlip(qureg, idQubit1, idQubit2);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledPhaseFlip(qureg, idQubit1+shift, idQubit2+shift);
    }
    
    qasm_recordControlledGate(qureg, GATE_SIGMA_Z, idQubit1, idQubit2);
}

void multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits) {
    validateMultiControls(qureg, controlQubits, numControlQubits, __func__);
    
    statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        shiftIndices(controlQubits, numControlQubits, shift);
        statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
        shiftIndices(controlQubits, numControlQubits, -shift);
    }
    
    qasm_recordMultiControlledGate(qureg, GATE_SIGMA_Z, controlQubits, numControlQubits-1, controlQubits[numControlQubits-1]);
}

void rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis) {
    validateTarget(qureg, rotQubit, __func__);
    validateVector(axis, __func__);
    
    statevec_rotateAroundAxis(qureg, rotQubit, angle, axis);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_rotateAroundAxisConj(qureg, rotQubit+shift, angle, axis);
    }
    
    qasm_recordAxisRotation(qureg, angle, axis, rotQubit);
}

void controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis) {
    validateControlTarget(qureg, controlQubit, targetQubit, __func__);
    validateVector(axis, __func__);
    
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis);
    if (qureg.isDensityMatrix) {
        int shift = qureg.numQubitsRepresented;
        statevec_controlledRotateAroundAxisConj(qureg, controlQubit+shift, targetQubit+shift, angle, axis);
    }
    
    qasm_recordControlledAxisRotation(qureg, angle, axis, controlQubit, targetQubit);
}

int getNumQubits(QubitRegister qureg) {
    return qureg.numQubitsRepresented;
}

int getNumAmps(QubitRegister qureg) {
    validateStateVecQureg(qureg, __func__);
    
    return qureg.numAmpsTotal;
}

REAL getRealAmpEl(QubitRegister qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateStateIndex(qureg, index, __func__);
    
    return statevec_getRealAmpEl(qureg, index);
}

REAL getImagAmpEl(QubitRegister qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateStateIndex(qureg, index, __func__);
    
    return statevec_getImagAmpEl(qureg, index);
}

REAL getProbEl(QubitRegister qureg, long long int index) {
    validateStateVecQureg(qureg, __func__);
    validateStateIndex(qureg, index, __func__);
    
    return statevec_getProbEl(qureg, index);
}

int compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision) {
    validateStateVecQureg(mq1, __func__);
    validateStateVecQureg(mq2, __func__);
    
    return statevec_compareStates(mq1, mq2, precision);
}

REAL calcTotalProbability(QubitRegister qureg) {
    if (qureg.isDensityMatrix)  
            return densmatr_calcTotalProbability(qureg);
        else
            return statevec_calcTotalProbability(qureg);
}

REAL findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
    validateTarget(qureg, measureQubit, __func__); // should rename? meh
    validateOutcome(outcome, __func__);
    
    if (qureg.isDensityMatrix)
        return densmatr_findProbabilityOfOutcome(qureg, measureQubit, outcome);
    else
        return statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
}

void cloneQubitRegister(QubitRegister targetQureg, QubitRegister copyQureg) {
    validateMatchingQuregTypes(targetQureg, copyQureg, __func__);
    validateMatchingQuregDims(targetQureg, copyQureg, __func__);
    
    statevec_cloneQubitRegister(targetQureg, copyQureg);
}





// @TODO add density copying to distributed CPU
void initPureState(QubitRegister qureg, QubitRegister pure) {
    validateSecondQuregStateVec(pure, __func__);
    validateMatchingQuregDims(qureg, pure, __func__);

    if (qureg.isDensityMatrix)
        densmatr_initPureState(qureg, pure);
    else
        statevec_cloneQubitRegister(qureg, pure);
    
    // @TODO: QASM?
}









// @TODO implement CPU (local & MPI) densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
REAL collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
    validateTarget(qureg, measureQubit, __func__); // should rename? eh
    validateOutcome(outcome, __func__);
    
    REAL outcomeProb;
    if (qureg.isDensityMatrix) {
        outcomeProb = densmatr_findProbabilityOfOutcome(qureg, measureQubit, outcome);
        validateMeasurementProb(outcomeProb, __func__);
        densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
    } else {
        outcomeProb = statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
        validateMeasurementProb(outcomeProb, __func__);
        statevec_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
    }
    
    // @TODO: add QASM for post-selecting the outcome?
    qasm_recordMeasurement(qureg, measureQubit);

    return outcomeProb;
}

int measureWithStats(QubitRegister qureg, int measureQubit, REAL *outcomeProb) {
    validateTarget(qureg, measureQubit, __func__); // should rename? eh

    int outcome;
    if (qureg.isDensityMatrix)
        outcome = densmatr_measureWithStats(qureg, measureQubit, outcomeProb);
    else
        outcome = statevec_measureWithStats(qureg, measureQubit, outcomeProb);
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcome;
}

int measure(QubitRegister qureg, int measureQubit) {
    validateTarget(qureg, measureQubit, __func__); // should rename? eh
    
    int outcome;
    REAL discardedProb;
    if (qureg.isDensityMatrix)
        outcome = densmatr_measureWithStats(qureg, measureQubit, &discardedProb);
    else
        outcome = statevec_measureWithStats(qureg, measureQubit, &discardedProb);
    
    qasm_recordMeasurement(qureg, measureQubit);
    return outcome;
}





// new experimental dephasing functions

// @TODO add to CPU local and distributed
void oneQubitDephase(QubitRegister qureg, const int targetQubit, REAL dephase) {
    validatDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateNoise(dephase, __func__);
    
    densmatr_oneQubitDephase(qureg, targetQubit, dephase);
}

// @TODO add to CPU local and distributed
void twoQubitDephase(QubitRegister qureg, const int qubit1, const int qubit2, REAL dephase) {
    validatDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, qubit1, __func__);
    validateTarget(qureg, qubit2, __func__);
    validateNoise(dephase, __func__);

    densmatr_twoQubitDephase(qureg, qubit1, qubit2, dephase);
}

// @TODO add to CPU local and distributed
void oneQubitDepolarise(QubitRegister qureg, const int targetQubit, REAL depolLevel) {
    validatDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, targetQubit, __func__);
    validateNoise(depolLevel, __func__);
    
    densmatr_oneQubitDepolarise(qureg, targetQubit, depolLevel);
}

// @TODO add to CPU local and distributed
void twoQubitDepolarise(QubitRegister qureg, const int qubit1, const int qubit2, REAL depolLevel) {
    validatDensityMatrQureg(qureg, __func__);
    validateTarget(qureg, qubit1, __func__);
    validateTarget(qureg, qubit2, __func__);
    validateNoise(depolLevel, __func__);
    
    densmatr_twoQubitDepolarise(qureg, qubit1, qubit2, depolLevel);
}

// @TODO add to CPU local and distributed
void combineDensityMatrices(REAL combineProb, QubitRegister combineQureg, REAL otherProb, QubitRegister otherQureg) {
    validatDensityMatrQureg(combineQureg, __func__);
    validatDensityMatrQureg(otherQureg, __func__);
    validateMatchingQuregDims(combineQureg, otherQureg, __func__);
    validateNormProbs(combineProb, otherProb, __func__);
    
    densmatr_combineDensityMatrices(combineProb, combineQureg, otherProb, otherQureg);
}






// @TODO
void initStateDebug(QubitRegister qureg) {
    statevec_initStateDebug(qureg);
}

// @TODO
void initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env) {
    
    int success = 0;
    
    // @TODO allow density matrix loading from file
    if (qureg->isDensityMatrix)
        validateStateVecQureg(*qureg, __func__);
    
    else
        success = statevec_initStateFromSingleFile(qureg, filename, env);
    
    validateFileOpened(success, __func__);
}

// @TODO
void initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome) {
    return statevec_initStateOfSingleQubit(qureg, qubitId, outcome);
}

// @TODO
void reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank)  {
    statevec_reportStateToScreen(qureg, env, reportRank);
}






#ifdef __cplusplus
}
#endif