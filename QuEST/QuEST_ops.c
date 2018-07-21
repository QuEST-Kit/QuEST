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
// 		- pure_initPureState on GPU
// 		- mixed_initPureState on GPU
// 		- mixed_initPureStateDistributed on CPU

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "QuEST_ops_pure.h"
# include "QuEST_ops_mixed.h"

#ifdef __cplusplus
extern "C" {
#endif

void createDensityQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
	pure_createQubitRegister(qureg, 2*numQubits, env);
	qureg->isDensityMatrix = 1;
	qureg->numDensityQubits = numQubits;
}

void createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
	pure_createQubitRegister(qureg, numQubits, env);
	qureg->isDensityMatrix = 0;
	qureg->numDensityQubits = -1;
}

void destroyQubitRegister(QubitRegister qureg, QuESTEnv env) {
	pure_destroyQubitRegister(qureg, env);
}

void initStateZero(QubitRegister qureg) {
	pure_initStateZero(qureg); // valid for |0> and |0><0|
}

void initStatePlus(QubitRegister qureg) {
	if (qureg.isDensityMatrix)
		mixed_initStatePlus(qureg);
	else
		pure_initStatePlus(qureg);
}

void initClassicalState(QubitRegister qureg, long long int stateInd) {
	if (qureg.isDensityMatrix)
		mixed_initClassicalState(qureg, stateInd);
	else
		pure_initClassicalState(qureg, stateInd);
}

// @TODO add pure copying to GPU
// @TODO add density copying to distributed CPU
// @TODO add density copying to GPU
void initPureState(QubitRegister qureg, QubitRegister pure) {
	QuESTAssert(!pure.isDensityMatrix, 12, __func__);
	
	if (qureg.isDensityMatrix) {
		QuESTAssert(qureg.numDensityQubits==pure.numQubits, 13, __func__);
		mixed_initPureState(qureg, pure);
		
	} else {
		QuESTAssert(qureg.numQubits==pure.numQubits, 13, __func__);
		pure_initPureState(qureg, pure);
	}
}

void hadamard(QubitRegister qureg, const int targetQubit) {
	pure_hadamard(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		pure_hadamard(qureg, targetQubit+qureg.numDensityQubits);
	}
}

void rotateX(QubitRegister qureg, const int rotQubit, REAL angle) {
	pure_rotateX(qureg, rotQubit, angle);
	if (qureg.isDensityMatrix) {
		pure_rotateX(qureg, rotQubit+qureg.numDensityQubits, -angle);
	}
}

void rotateY(QubitRegister qureg, const int rotQubit, REAL angle) {
	pure_rotateY(qureg, rotQubit, angle);
	if (qureg.isDensityMatrix) {
		pure_rotateY(qureg, rotQubit+qureg.numDensityQubits, angle);
	}
}

void rotateZ(QubitRegister qureg, const int rotQubit, REAL angle) {
	pure_rotateZ(qureg, rotQubit, angle);
	if (qureg.isDensityMatrix) {
		pure_rotateZ(qureg, rotQubit+qureg.numDensityQubits, -angle);
	}
}

void controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateX(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledRotateX(qureg, controlQubit+shift, targetQubit+shift, -angle);
	}
}

void controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateY(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledRotateY(qureg, controlQubit+shift, targetQubit+shift, angle);
	}
}

void controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateZ(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledRotateZ(qureg, controlQubit+shift, targetQubit+shift, -angle);
	}
}

void unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u) {
	pure_unitary(qureg, targetQubit, u);
	if (qureg.isDensityMatrix) {
		pure_unitary(qureg, targetQubit+qureg.numDensityQubits, getConjugateMatrix(u));
	}
}

void controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) {
	pure_controlledUnitary(qureg, controlQubit, targetQubit, u);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledUnitary(qureg, controlQubit+shift, targetQubit+shift, getConjugateMatrix(u));
	}
}

void multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) {
	pure_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit, u);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		shiftIndices(controlQubits, numControlQubits, shift);
		pure_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit+shift, getConjugateMatrix(u));
		shiftIndices(controlQubits, numControlQubits, -shift);
	}
}

void compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) {
	pure_compactUnitary(qureg, targetQubit, alpha, beta);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_compactUnitary(qureg, targetQubit+shift, getConjugateScalar(alpha), getConjugateScalar(beta));
	}
}

void controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) {
	pure_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledCompactUnitary(qureg, 
			controlQubit+shift, targetQubit+shift, 
			getConjugateScalar(alpha), getConjugateScalar(beta));
	}
}

void sigmaX(QubitRegister qureg, const int targetQubit) {
	pure_sigmaX(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		pure_sigmaX(qureg, targetQubit+qureg.numDensityQubits);
	}
}

void sigmaZ(QubitRegister qureg, const int targetQubit) {
	pure_sigmaZ(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		pure_sigmaZ(qureg, targetQubit+qureg.numDensityQubits);
	}
}

void phaseShift(QubitRegister qureg, const int targetQubit, REAL angle) {
	pure_phaseShift(qureg, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		pure_phaseShift(qureg, targetQubit+qureg.numDensityQubits, -angle);
	}
}

void sGate(QubitRegister qureg, const int targetQubit) {
	pure_sGate(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		pure_sGateConj(qureg, targetQubit+qureg.numDensityQubits);
	}
}

void tGate(QubitRegister qureg, const int targetQubit) {
	pure_tGate(qureg, targetQubit);
	if (qureg.isDensityMatrix) {
		pure_tGateConj(qureg, targetQubit+qureg.numDensityQubits);
	}
}

void controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit) {
	pure_controlledNot(qureg, controlQubit, targetQubit);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledNot(qureg, controlQubit+shift, targetQubit+shift);
	}
}

void controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2) {
	pure_controlledPhaseFlip(qureg, idQubit1, idQubit2);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledPhaseFlip(qureg, idQubit1+shift, idQubit2+shift);
	}
}

void multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits) {
	pure_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		shiftIndices(controlQubits, numControlQubits, shift);
		pure_multiControlledPhaseFlip(qureg, controlQubits, numControlQubits);
		shiftIndices(controlQubits, numControlQubits, -shift);
	}
}

int getNumQubits(QubitRegister qureg) {
	if (qureg.isDensityMatrix)
		return qureg.numDensityQubits;
	else
		return pure_getNumQubits(qureg);
}

int compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision) {
	QuESTAssert(!mq1.isDensityMatrix && !mq2.isDensityMatrix, 15, __func__);
	return pure_compareStates(mq1, mq2, precision);
}

int getNumAmps(QubitRegister qureg) {
	QuESTAssert(!qureg.isDensityMatrix, 14, __func__);
	return pure_getNumAmps(qureg);
}

REAL getRealAmpEl(QubitRegister qureg, long long int index) {
	QuESTAssert(!qureg.isDensityMatrix, 14, __func__);
	return pure_getRealAmpEl(qureg, index);
}

REAL getImagAmpEl(QubitRegister qureg, long long int index) {
	QuESTAssert(!qureg.isDensityMatrix, 14, __func__);
	return pure_getImagAmpEl(qureg, index);
}

REAL getProbEl(QubitRegister qureg, long long int index) {
	QuESTAssert(!qureg.isDensityMatrix, 14, __func__);
	return pure_getProbEl(qureg, index);
}


// @TODO very dubious
void controlledPhaseShift(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledPhaseShift(qureg, controlQubit, targetQubit, angle);
	if (qureg.isDensityMatrix) {
		int shift = qureg.numDensityQubits;
		pure_controlledPhaseShift(qureg, controlQubit+shift, targetQubit+shift, -angle);
	}
}






// @TODO
void sigmaY(QubitRegister qureg, const int targetQubit) {
	pure_sigmaY(qureg, targetQubit);
}





// @TODO
void initStateDebug(QubitRegister qureg) {
	pure_initStateDebug(qureg);
}



// @TODO
void initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env) {
	pure_initStateFromSingleFile(qureg, filename, env);
}

// @TODO
void initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome) {
	return pure_initStateOfSingleQubit(qureg, qubitId, outcome);
}

// @TODO
void reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank)  {
	pure_reportStateToScreen(qureg, env, reportRank);
}





// @TODO
void rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis) {
	pure_rotateAroundAxis(qureg, rotQubit, angle, axis);
}

// @TODO
void controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis) {
	pure_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis);
}

// @TODO
REAL calcTotalProbability(QubitRegister qureg) {
	return pure_calcTotalProbability(qureg);
}

// @TODO
REAL findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	return pure_findProbabilityOfOutcome(qureg, measureQubit, outcome);
}

// @TODO
REAL collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	return pure_collapseToOutcome(qureg, measureQubit, outcome);
}

// @TODO
int measure(QubitRegister qureg, int measureQubit) {
	return pure_measure(qureg, measureQubit);
}

// @TODO
int measureWithStats(QubitRegister qureg, int measureQubit, REAL *stateProb) {
	return pure_measureWithStats(qureg, measureQubit, stateProb);
}






#ifdef __cplusplus
}
#endif