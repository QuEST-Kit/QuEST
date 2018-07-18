// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Implements the QuEST.h API, in a hardware-agnostic way, for both pure and mixed states
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_ops_pure.h"

#ifdef __cplusplus
extern "C" {
#endif



void createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env) {
	pure_createQubitRegister(qureg, numQubits, env);
}

void destroyQubitRegister(QubitRegister qureg, QuESTEnv env) {
	pure_destroyQubitRegister(qureg, env);
}

void initStateZero(QubitRegister qureg) {
	pure_initStateZero(qureg);
}

void initStatePlus(QubitRegister qureg) {
	pure_initStatePlus(qureg);
}

void initStateDebug(QubitRegister qureg) {
	pure_initStateDebug(qureg);
}

void initClassicalState(QubitRegister qureg, long long int stateInd) {
	pure_initClassicalState(qureg, stateInd);
}

void initializeStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env) {
	pure_initializeStateFromSingleFile(qureg, filename, env);
}

void reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank)  {
	pure_reportStateToScreen(qureg, env, reportRank);
}

void phaseGate(QubitRegister qureg, const int targetQubit, enum phaseGateType type) {
	pure_phaseGate(qureg, targetQubit, type);
}

void multiControlledPhaseGate(QubitRegister qureg, int *controlQubits, int numControlQubits) {
	pure_multiControlledPhaseGate(qureg, controlQubits, numControlQubits);
}

void controlledPhaseGate (QubitRegister qureg, const int idQubit1, const int idQubit2) {
	pure_controlledPhaseGate (qureg, idQubit1, idQubit2);
}

void sGate(QubitRegister qureg, const int targetQubit) {
	pure_sGate(qureg, targetQubit);
}

void tGate(QubitRegister qureg, const int targetQubit) {
	pure_tGate(qureg, targetQubit);
}

void compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) {
	pure_compactUnitary(qureg, targetQubit, alpha, beta);
}

void unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u) {
	pure_unitary(qureg, targetQubit, u);
}

void rotateX(QubitRegister qureg, const int rotQubit, REAL angle) {
	pure_rotateX(qureg, rotQubit, angle);
}

void rotateY(QubitRegister qureg, const int rotQubit, REAL angle) {
	pure_rotateY(qureg, rotQubit, angle);
}

void rotateZ(QubitRegister qureg, const int rotQubit, REAL angle) {
	pure_rotateZ(qureg, rotQubit, angle);
}

void rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis) {
	pure_rotateAroundAxis(qureg, rotQubit, angle, axis);
}

void controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateX(qureg, controlQubit, targetQubit, angle);
}

void controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateY(qureg, controlQubit, targetQubit, angle);
}

void controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateZ(qureg, controlQubit, targetQubit, angle);
}

void controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis) {
	pure_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, axis);
}

void controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) {
	pure_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
}

void controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) {
	pure_controlledUnitary(qureg, controlQubit, targetQubit, u);
}

void multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) {
	pure_multiControlledUnitary(qureg, controlQubits, numControlQubits, targetQubit, u);
}

void sigmaX(QubitRegister qureg, const int targetQubit) {
	pure_sigmaX(qureg, targetQubit);
}

void sigmaY(QubitRegister qureg, const int targetQubit) {
	pure_sigmaY(qureg, targetQubit);
}

void sigmaZ(QubitRegister qureg, const int targetQubit) {
	pure_sigmaZ(qureg, targetQubit);
}

void hadamard(QubitRegister qureg, const int targetQubit) {
	pure_hadamard(qureg, targetQubit);
}

void controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit) {
	pure_controlledNot(qureg, controlQubit, targetQubit);
}


int compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision){
	return pure_compareStates(mq1, mq2, precision);
}
void initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome) {
	return pure_initStateOfSingleQubit(qureg, qubitId, outcome);
}

int getNumQubits(QubitRegister qureg){
	return pure_getNumQubits(qureg);
}

int getNumAmps(QubitRegister qureg) {
	return pure_getNumAmps(qureg);
}

REAL getRealAmpEl(QubitRegister qureg, long long int index) {
	return pure_getRealAmpEl(qureg, index);
}

REAL getImagAmpEl(QubitRegister qureg, long long int index) {
	return pure_getImagAmpEl(qureg, index);
}

REAL getProbEl(QubitRegister qureg, long long int index) {
	return pure_getProbEl(qureg, index);
}

REAL calcTotalProbability(QubitRegister qureg) {
	return pure_calcTotalProbability(qureg);
}

REAL findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	return pure_findProbabilityOfOutcome(qureg, measureQubit, outcome);
}

REAL collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	return pure_collapseToOutcome(qureg, measureQubit, outcome);
}

int measure(QubitRegister qureg, int measureQubit) {
	return pure_measure(qureg, measureQubit);
}

int measureWithStats(QubitRegister qureg, int measureQubit, REAL *stateProb) {
	return pure_measureWithStats(qureg, measureQubit, stateProb);
}



#ifdef __cplusplus
}
#endif