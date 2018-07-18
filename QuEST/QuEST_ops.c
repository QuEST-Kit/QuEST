
# include "QuEST.h"
# include "QuEST_precision.h"

#ifdef __cplusplus
extern "C" {
#endif

	
int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision){
	return 0;
}
void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome) {}
void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env){}
void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank) {}






void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type){}

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env){}

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env){}

int getNumQubits(MultiQubit multiQubit){ return -1; }

int getNumAmps(MultiQubit multiQubit){ return -1; }

void initStateZero(MultiQubit multiQubit){}

void initStatePlus(MultiQubit multiQubit){}

void initStateDebug(MultiQubit multiQubit){}

void initClassicalState(MultiQubit multiQubit, long long int stateInd){}

void multiControlledPhaseGate(MultiQubit multiQubit, int *controlQubits, int numControlQubits){}

void controlledPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2){}

void sGate(MultiQubit multiQubit, const int targetQubit){}

void tGate(MultiQubit multiQubit, const int targetQubit){}

REAL getRealAmpEl(MultiQubit multiQubit, long long int index){ return -1; }

REAL getImagAmpEl(MultiQubit multiQubit, long long int index){ return -1; }

REAL getProbEl(MultiQubit multiQubit, long long int index){ return -1; }

REAL calcTotalProbability(MultiQubit multiQubit){ return -1; }

void compactUnitary(MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta){}

void unitary(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u){}

void rotateX(MultiQubit multiQubit, const int rotQubit, REAL angle){}

void rotateY(MultiQubit multiQubit, const int rotQubit, REAL angle){}

void rotateZ(MultiQubit multiQubit, const int rotQubit, REAL angle){}

void rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector axis){}

void controlledRotateX(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){}

void controlledRotateY(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){}

void controlledRotateZ(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){}

void controlledRotateAroundAxis(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle, Vector axis){}

void controlledCompactUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta){}

void controlledUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u){}

void multiControlledUnitary(MultiQubit multiQubit, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u){}

void sigmaX(MultiQubit multiQubit, const int targetQubit){}

void sigmaY(MultiQubit multiQubit, const int targetQubit){}

void sigmaZ(MultiQubit multiQubit, const int targetQubit){}

void hadamard(MultiQubit multiQubit, const int targetQubit){}

void controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit){}

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome){ return -1; }

REAL collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome){ return -1; }

int measure(MultiQubit multiQubit, int measureQubit){ return -1; }

int measureWithStats(MultiQubit multiQubit, int measureQubit, REAL *stateProb){ return -1; }

#ifdef __cplusplus
}
#endif