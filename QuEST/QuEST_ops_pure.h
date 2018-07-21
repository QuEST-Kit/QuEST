// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Functions operating upon pure states which are provided by a hardware-specific backend
 */

# ifndef QuEST_OPS_PURE
# define QuEST_OPS_PURE

# include "QuEST.h"
# include "QuEST_precision.h"

# ifdef __cplusplus
extern "C" {
# endif
	
void pure_reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank);

int pure_compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision);

void pure_initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env);

void pure_initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome);

void pure_createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env);

void pure_destroyQubitRegister(QubitRegister qureg, QuESTEnv env);

int pure_getNumQubits(QubitRegister qureg);

int pure_getNumAmps(QubitRegister qureg);

void pure_initStateZero(QubitRegister qureg);

void pure_initStatePlus(QubitRegister qureg);

void pure_initStateDebug(QubitRegister qureg);

void pure_initClassicalState(QubitRegister qureg, long long int stateInd);

void pure_initPureState(QubitRegister targetQureg, QubitRegister copyQureg);

void pure_multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits);

void pure_controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2);

void pure_phaseShift(QubitRegister qureg, const int targetQubit, REAL angle);

void pure_controlledPhaseShift(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle); 

void pure_sGate(QubitRegister qureg, const int targetQubit);

void pure_tGate(QubitRegister qureg, const int targetQubit);

void pure_sGateConj(QubitRegister qureg, const int targetQubit);

void pure_tGateConj(QubitRegister qureg, const int targetQubit);

void pure_sigmaX(QubitRegister qureg, const int targetQubit);

void pure_sigmaY(QubitRegister qureg, const int targetQubit);

void pure_sigmaZ(QubitRegister qureg, const int targetQubit);

REAL pure_getRealAmpEl(QubitRegister qureg, long long int index);

REAL pure_getImagAmpEl(QubitRegister qureg, long long int index);

REAL pure_getProbEl(QubitRegister qureg, long long int index);

REAL pure_calcTotalProbability(QubitRegister qureg);

void pure_compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta);

void pure_unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u);

void pure_rotateX(QubitRegister qureg, const int rotQubit, REAL angle);

void pure_rotateY(QubitRegister qureg, const int rotQubit, REAL angle);

void pure_rotateZ(QubitRegister qureg, const int rotQubit, REAL angle);

void pure_rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis);

void pure_controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle);

void pure_controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle);

void pure_controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle);

void pure_controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis);

void pure_controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta);

void pure_controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void pure_multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u);

void pure_hadamard(QubitRegister qureg, const int targetQubit);

void pure_controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit);

REAL pure_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome);

REAL pure_collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome);

int pure_measure(QubitRegister qureg, int measureQubit);

int pure_measureWithStats(QubitRegister qureg, int measureQubit, REAL *stateProb);

# ifdef __cplusplus
}
# endif

# endif //QuEST_OPS_PURE