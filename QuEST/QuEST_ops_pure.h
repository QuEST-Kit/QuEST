
# include "QuEST.h"
# include "QuEST_precision.h"

#ifdef __cplusplus
extern "C" {
#endif

void pure_initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env);
	
void pure_reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank);

int pure_compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision);

void pure_initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome);

void pure_phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);

void pure_createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env);

void pure_destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env);

int pure_getNumQubits(MultiQubit multiQubit);

int pure_getNumAmps(MultiQubit multiQubit);

void pure_initStateZero(MultiQubit multiQubit);

void pure_initStatePlus(MultiQubit multiQubit);

void pure_initStateDebug(MultiQubit multiQubit);

void pure_initClassicalState(MultiQubit multiQubit, long long int stateInd);

void pure_multiControlledPhaseGate(MultiQubit multiQubit, int *controlQubits, int numControlQubits);

void pure_controlledPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2);

void pure_sGate(MultiQubit multiQubit, const int targetQubit);

void pure_tGate(MultiQubit multiQubit, const int targetQubit);

REAL pure_getRealAmpEl(MultiQubit multiQubit, long long int index);

REAL pure_getImagAmpEl(MultiQubit multiQubit, long long int index);

REAL pure_getProbEl(MultiQubit multiQubit, long long int index);

REAL pure_calcTotalProbability(MultiQubit multiQubit);

void pure_compactUnitary(MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta);

void pure_unitary(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u);

void pure_rotateX(MultiQubit multiQubit, const int rotQubit, REAL angle);

void pure_rotateY(MultiQubit multiQubit, const int rotQubit, REAL angle);

void pure_rotateZ(MultiQubit multiQubit, const int rotQubit, REAL angle);

void pure_rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector axis);

void pure_controlledRotateX(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle);

void pure_controlledRotateY(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle);

void pure_controlledRotateZ(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle);

void pure_controlledRotateAroundAxis(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle, Vector axis);

void pure_controlledCompactUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta);

void pure_controlledUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void pure_multiControlledUnitary(MultiQubit multiQubit, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u);

void pure_sigmaX(MultiQubit multiQubit, const int targetQubit);

void pure_sigmaY(MultiQubit multiQubit, const int targetQubit);

void pure_sigmaZ(MultiQubit multiQubit, const int targetQubit);

void pure_hadamard(MultiQubit multiQubit, const int targetQubit);

void pure_controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit);

REAL pure_findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome);

REAL pure_collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome);

int pure_measure(MultiQubit multiQubit, int measureQubit);

int pure_measureWithStats(MultiQubit multiQubit, int measureQubit, REAL *stateProb);


#ifdef __cplusplus
}
#endif