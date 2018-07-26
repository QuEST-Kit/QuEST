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
// 		- densmatr_initPureStateDistributed on CPU

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "QuEST_ops.h"


// just for debug printing and exiting: can remove later
# include <stdio.h>
# include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif

void createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
	validateCreateNumQubits(numQubits);
	
	statevec_createQubitRegister(qureg, numQubits, env);
	qureg->isDensityMatrix = 0;
	qureg->numQubitsRepresented = numQubits;
	qureg->numQubitsInStateVec = numQubits;
}

void createDensityQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
	validateCreateNumQubits(numQubits);
	
	statevec_createQubitRegister(qureg, 2*numQubits, env);
	qureg->isDensityMatrix = 1;
	qureg->numQubitsRepresented = numQubits;
	qureg->numQubitsInStateVec = 2*numQubits;
}

void destroyQubitRegister(QubitRegister qureg, QuESTEnv env) {
	statevec_destroyQubitRegister(qureg, env);
}

void initStateZero(QubitRegister qureg) {
	statevec_initStateZero(qureg); // valid for both statevec and density matrices
}

void initStatePlus(QubitRegister qureg) {
	if (qureg.isDensityMatrix)
		densmatr_initStatePlus(qureg);
	else
		statevec_initStatePlus(qureg);
}

void initClassicalState(QubitRegister qureg, long long int stateInd) {
	validateStateIndex(qureg, stateInd);
	
	if (qureg.isDensityMatrix)
		densmatr_initClassicalState(qureg, stateInd);
	else
		statevec_initClassicalState(qureg, stateInd);
}

void hadamard(QubitRegister qureg, const int targetQubit) {
	validateTarget(qureg, targetQubit);

	statevec_hadamard(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		statevec_hadamard(qureg, targetQubit+qureg.numQubitsRepresented);
	}
}

void rotateX(QubitRegister qureg, const int targetQubit, REAL angle) {
	validateTarget(qureg, targetQubit);
	
	statevec_rotateX(qureg, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		statevec_rotateX(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
	}
}

void rotateY(QubitRegister qureg, const int targetQubit, REAL angle) {
	validateTarget(qureg, targetQubit);
	
	statevec_rotateY(qureg, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		statevec_rotateY(qureg, targetQubit+qureg.numQubitsRepresented, angle);
	}
}

void rotateZ(QubitRegister qureg, const int targetQubit, REAL angle) {
	validateTarget(qureg, targetQubit);
	
	statevec_rotateZ(qureg, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		statevec_rotateZ(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
	}
}

void controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	
	statevec_controlledRotateX(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledRotateX(qureg, controlQubit+shift, targetQubit+shift, -angle);
	}
}

void controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	
	statevec_controlledRotateY(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledRotateY(qureg, controlQubit+shift, targetQubit+shift, angle);
	}
}

void controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	
	statevec_controlledRotateZ(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledRotateZ(qureg, controlQubit+shift, targetQubit+shift, -angle);
	}
}

void unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u) {
	validateTarget(qureg, targetQubit);
	validateUnitaryMatrix(u);
	
	statevec_unitary(qureg, targetQubit, u);
	if (qureg.isDensityMatrix) {
		statevec_unitary(qureg, targetQubit+qureg.numQubitsRepresented, getConjugateMatrix(u));
	}
}

void controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	validateUnitaryMatrix(u);
	
	statevec_controlledUnitary(qureg, controlQubit, targetQubit, u);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledUnitary(qureg, controlQubit+shift, targetQubit+shift, getConjugateMatrix(u));
	}
}

void multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) {
	validateMultiControlsTarget(qureg, controlQubits, numControlQubits, targetQubit);
	validateUnitaryMatrix(u);
	
	statevec_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit, u);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		shiftIndices(controlQubits, numControlQubits, shift);
		statevec_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit+shift, getConjugateMatrix(u));
		shiftIndices(controlQubits, numControlQubits, -shift);
	}
}

void compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) {
	validateTarget(qureg, targetQubit);
	validateUnitaryComplexPair(alpha, beta);
	
	statevec_compactUnitary(qureg, targetQubit, alpha, beta);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_compactUnitary(qureg, targetQubit+shift, getConjugateScalar(alpha), getConjugateScalar(beta));
	}
}

void controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	validateUnitaryComplexPair(alpha, beta);
	
	statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledCompactUnitary(qureg, 
			controlQubit+shift, targetQubit+shift, 
			getConjugateScalar(alpha), getConjugateScalar(beta));
	}
}

void sigmaX(QubitRegister qureg, const int targetQubit) {
	validateTarget(qureg, targetQubit);
	
	statevec_sigmaX(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		statevec_sigmaX(qureg, targetQubit+qureg.numQubitsRepresented);
	}
}

void sigmaY(QubitRegister qureg, const int targetQubit) {
	validateTarget(qureg, targetQubit);
	
	statevec_sigmaY(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		statevec_sigmaYConj(qureg, targetQubit + qureg.numQubitsRepresented);
	}
}

void sigmaZ(QubitRegister qureg, const int targetQubit) {
	validateTarget(qureg, targetQubit);
	
	statevec_sigmaZ(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		statevec_sigmaZ(qureg, targetQubit+qureg.numQubitsRepresented);
	}
}

void sGate(QubitRegister qureg, const int targetQubit) {
	validateTarget(qureg, targetQubit);
	
	statevec_sGate(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		statevec_sGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
	}
}

void tGate(QubitRegister qureg, const int targetQubit) {
	validateTarget(qureg, targetQubit);
	
	statevec_tGate(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		statevec_tGateConj(qureg, targetQubit+qureg.numQubitsRepresented);
	}
}

void phaseShift(QubitRegister qureg, const int targetQubit, REAL angle) {
	validateTarget(qureg, targetQubit);
	
	statevec_phaseShift(qureg, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		statevec_phaseShift(qureg, targetQubit+qureg.numQubitsRepresented, -angle);
	}
}

void controlledPhaseShift(QubitRegister qureg, const int idQubit1, const int idQubit2, REAL angle) {
	validateControlTarget(qureg, idQubit1, idQubit2); // a little bit semantically dodgy
	
	statevec_controlledPhaseShift(qureg, idQubit1, idQubit2, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledPhaseShift(qureg, idQubit1+shift, idQubit2+shift, -angle);
	}
}

void multiControlledPhaseShift(QubitRegister qureg, int *controlQubits, int numControlQubits, REAL angle) {
	validateMultiControls(qureg, controlQubits, numControlQubits);
	
	statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		shiftIndices(controlQubits, numControlQubits, shift);
		statevec_multiControlledPhaseShift(qureg, controlQubits, numControlQubits, angle);
		shiftIndices(controlQubits, numControlQubits, -shift);
	}
}

void controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	
	statevec_controlledNot(qureg, controlQubit, targetQubit);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledNot(qureg, controlQubit+shift, targetQubit+shift);
	}
}

void controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	
	statevec_controlledSigmaY(qureg, controlQubit, targetQubit);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledSigmaYConj(qureg, controlQubit+shift, targetQubit+shift);
	}
}

void controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2) {
	validateControlTarget(qureg, idQubit1, idQubit2); // a little bit semantically dodgy
	
	statevec_controlledPhaseFlip(qureg, idQubit1, idQubit2);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledPhaseFlip(qureg, idQubit1+shift, idQubit2+shift);
	}
}

void multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits) {
	validateMultiControls(qureg, controlQubits, numControlQubits);
	
	statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		shiftIndices(controlQubits, numControlQubits, shift);
		statevec_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
		shiftIndices(controlQubits, numControlQubits, -shift);
	}
}

void rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis) {
	validateTarget(qureg, rotQubit);
	validateVector(axis);
	
	statevec_rotateAroundAxis(qureg, rotQubit, angle, axis);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_rotateAroundAxisConj(qureg, rotQubit+shift, angle, axis);
	}
}

void controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis) {
	validateControlTarget(qureg, controlQubit, targetQubit);
	validateVector(axis);
	
	statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numQubitsRepresented;
		statevec_controlledRotateAroundAxisConj(qureg, controlQubit+shift, targetQubit+shift, angle, axis);
	}
}

int getNumQubits(QubitRegister qureg) {
	return qureg.numQubitsRepresented;
}

int getNumAmps(QubitRegister qureg) {
	validateStateVecQureg(qureg);
	
	return qureg.numAmpsTotal;
}

REAL getRealAmpEl(QubitRegister qureg, long long int index) {
	validateStateVecQureg(qureg);
	validateStateIndex(qureg, index);
	
	return statevec_getRealAmpEl(qureg, index);
}

REAL getImagAmpEl(QubitRegister qureg, long long int index) {
	validateStateVecQureg(qureg);
	validateStateIndex(qureg, index);
	
	return statevec_getImagAmpEl(qureg, index);
}

REAL getProbEl(QubitRegister qureg, long long int index) {
	validateStateVecQureg(qureg);
	validateStateIndex(qureg, index);
	
	return statevec_getProbEl(qureg, index);
}

int compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision) {
	validateStateVecQureg(mq1);
	validateStateVecQureg(mq2);
	
	return statevec_compareStates(mq1, mq2, precision);
}

REAL calcTotalProbability(QubitRegister qureg) {
	if (qureg.isDensityMatrix)	
			return densmatr_calcTotalProbability(qureg);
		else
			return statevec_calcTotalProbability(qureg);
}

REAL findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	validateTarget(qureg, measureQubit); // should rename? meh
	validateOutcome(outcome);
	
	if (qureg.isDensityMatrix)
		return densmatr_findProbabilityOfOutcome(qureg, measureQubit, outcome);
	else
		return statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
}









// @TODO add density copying to distributed CPU
void initPureState(QubitRegister qureg, QubitRegister pure) {
	
	// validation here that second arg is a state-vector
	
	if (qureg.isDensityMatrix) {
		QuESTAssert(qureg.numQubitsRepresented==pure.numQubitsInStateVec, 13, __func__);
		densmatr_initPureState(qureg, pure);
		
	} else {
		QuESTAssert(qureg.numQubitsInStateVec==pure.numQubitsInStateVec, 13, __func__);
		statevec_initPureState(qureg, pure);
	}
}





/* simons requested func for cloning density matrices */





// @TODO implement densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
REAL collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	validateTarget(qureg, measureQubit); // should rename? eh
	validateOutcome(outcome);
	
	REAL outcomeProb;
	if (qureg.isDensityMatrix) {
		outcomeProb = densmatr_findProbabilityOfOutcome(qureg, measureQubit, outcome);
		validateMeasurementProb(outcomeProb);
		densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
	} else {
		outcomeProb = statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
		validateMeasurementProb(outcomeProb);
		statevec_collapseToKnownProbOutcome(qureg, measureQubit, outcome, outcomeProb);
	}

	return outcomeProb;
}

int measureWithStats(QubitRegister qureg, int measureQubit, REAL *outcomeProb) {
	validateTarget(qureg, measureQubit); // should rename? eh

	if (qureg.isDensityMatrix)
		return densmatr_measureWithStats(qureg, measureQubit, outcomeProb);
	else
		return statevec_measureWithStats(qureg, measureQubit, outcomeProb);
}


int measure(QubitRegister qureg, int measureQubit) {
	validateTarget(qureg, measureQubit); // should rename? eh
	
	REAL discardedProb;
	if (qureg.isDensityMatrix)
		return densmatr_measureWithStats(qureg, measureQubit, &discardedProb);
	else
		return statevec_measureWithStats(qureg, measureQubit, &discardedProb);
}









// @TODO
void initStateDebug(QubitRegister qureg) {
	statevec_initStateDebug(qureg);
}

// @TODO
void initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env) {
	statevec_initStateFromSingleFile(qureg, filename, env);
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