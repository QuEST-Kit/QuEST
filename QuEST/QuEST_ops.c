
# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_ops_pure.h"

#ifdef __cplusplus
extern "C" {
#endif


	/*
for e in y:
	sig = e[5:-2]
	name, allargs = sig.replace(')','').strip().split('(')
	args = allargs.split(', ')
	params = [arg.split(' ')[-1].replace('*','') for arg in args]
	call = 'pure_' + name + '(' + ', '.join(params) + ');'
	print(e[:-2] + ' {\n\t' + call + '\n}\n')
	*/

	
int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision){
	return pure_compareStates(mq1, mq2, precision);
}
void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome) {
	return pure_initStateOfSingleQubit(multiQubit, qubitId, outcome);
}

int getNumQubits(MultiQubit multiQubit){
	return pure_getNumQubits(multiQubit);
}

int getNumAmps(MultiQubit multiQubit) {
	return pure_getNumAmps(multiQubit);
}

REAL getRealAmpEl(MultiQubit multiQubit, long long int index) {
	return pure_getRealAmpEl(multiQubit, index);
}

REAL getImagAmpEl(MultiQubit multiQubit, long long int index) {
	return pure_getImagAmpEl(multiQubit, index);
}

REAL getProbEl(MultiQubit multiQubit, long long int index) {
	return pure_getProbEl(multiQubit, index);
}

REAL calcTotalProbability(MultiQubit multiQubit) {
	return pure_calcTotalProbability(multiQubit);
}

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome) {
	return pure_findProbabilityOfOutcome(multiQubit, measureQubit, outcome);
}

REAL collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome) {
	return pure_collapseToOutcome(multiQubit, measureQubit, outcome);
}

int measure(MultiQubit multiQubit, int measureQubit) {
	return pure_measure(multiQubit, measureQubit);
}

int measureWithStats(MultiQubit multiQubit, int measureQubit, REAL *stateProb) {
	return pure_measureWithStats(multiQubit, measureQubit, stateProb);
}










void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env) {
	pure_initializeStateFromSingleFile(multiQubit, filename, env);
}

void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank)  {
	pure_reportStateToScreen(multiQubit, env, reportRank);
}

void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type) {
	pure_phaseGate(multiQubit, targetQubit, type);
}

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env) {
	pure_createMultiQubit(multiQubit, numQubits, env);
}

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env) {
	pure_destroyMultiQubit(multiQubit, env);
}

void initStateZero(MultiQubit multiQubit) {
	pure_initStateZero(multiQubit);
}

void initStatePlus(MultiQubit multiQubit) {
	pure_initStatePlus(multiQubit);
}

void initStateDebug(MultiQubit multiQubit) {
	pure_initStateDebug(multiQubit);
}

void initClassicalState(MultiQubit multiQubit, long long int stateInd) {
	pure_initClassicalState(multiQubit, stateInd);
}

void multiControlledPhaseGate(MultiQubit multiQubit, int *controlQubits, int numControlQubits) {
	pure_multiControlledPhaseGate(multiQubit, controlQubits, numControlQubits);
}

void controlledPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2) {
	pure_controlledPhaseGate (multiQubit, idQubit1, idQubit2);
}

void sGate(MultiQubit multiQubit, const int targetQubit) {
	pure_sGate(multiQubit, targetQubit);
}

void tGate(MultiQubit multiQubit, const int targetQubit) {
	pure_tGate(multiQubit, targetQubit);
}

void compactUnitary(MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta) {
	pure_compactUnitary(multiQubit, targetQubit, alpha, beta);
}

void unitary(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u) {
	pure_unitary(multiQubit, targetQubit, u);
}

void rotateX(MultiQubit multiQubit, const int rotQubit, REAL angle) {
	pure_rotateX(multiQubit, rotQubit, angle);
}

void rotateY(MultiQubit multiQubit, const int rotQubit, REAL angle) {
	pure_rotateY(multiQubit, rotQubit, angle);
}

void rotateZ(MultiQubit multiQubit, const int rotQubit, REAL angle) {
	pure_rotateZ(multiQubit, rotQubit, angle);
}

void rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector axis) {
	pure_rotateAroundAxis(multiQubit, rotQubit, angle, axis);
}

void controlledRotateX(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateX(multiQubit, controlQubit, targetQubit, angle);
}

void controlledRotateY(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateY(multiQubit, controlQubit, targetQubit, angle);
}

void controlledRotateZ(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle) {
	pure_controlledRotateZ(multiQubit, controlQubit, targetQubit, angle);
}

void controlledRotateAroundAxis(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle, Vector axis) {
	pure_controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, axis);
}

void controlledCompactUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) {
	pure_controlledCompactUnitary(multiQubit, controlQubit, targetQubit, alpha, beta);
}

void controlledUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u) {
	pure_controlledUnitary(multiQubit, controlQubit, targetQubit, u);
}

void multiControlledUnitary(MultiQubit multiQubit, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) {
	pure_multiControlledUnitary(multiQubit, controlQubits, numControlQubits, targetQubit, u);
}

void sigmaX(MultiQubit multiQubit, const int targetQubit) {
	pure_sigmaX(multiQubit, targetQubit);
}

void sigmaY(MultiQubit multiQubit, const int targetQubit) {
	pure_sigmaY(multiQubit, targetQubit);
}

void sigmaZ(MultiQubit multiQubit, const int targetQubit) {
	pure_sigmaZ(multiQubit, targetQubit);
}

void hadamard(MultiQubit multiQubit, const int targetQubit) {
	pure_hadamard(multiQubit, targetQubit);
}

void controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit) {
	pure_controlledNot(multiQubit, controlQubit, targetQubit);
}




#ifdef __cplusplus
}
#endif