// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Functions operating upon pure states which are provided by a hardware-specific backend
 */

# ifndef QuEST_OPS
# define QuEST_OPS

# include "QuEST.h"
# include "QuEST_precision.h"

# ifdef __cplusplus
extern "C" {
# endif

/* operations upon state-vectors */

void densmatr_initStatePlus(QubitRegister targetQureg);

void densmatr_initClassicalState(QubitRegister qureg, long long int stateInd);

void densmatr_initPureState(QubitRegister targetQureg, QubitRegister copyQureg);

REAL densmatr_calcTotalProbability(QubitRegister qureg);

REAL densmatr_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome);
	
/* operations upon density matrices */
	
void statevec_reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank);

int statevec_compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision);

void statevec_initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env);

void statevec_initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome);

void statevec_createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env);

void statevec_destroyQubitRegister(QubitRegister qureg, QuESTEnv env);

void statevec_initStateZero(QubitRegister qureg);

void statevec_initStatePlus(QubitRegister qureg);

void statevec_initStateDebug(QubitRegister qureg);

void statevec_initClassicalState(QubitRegister qureg, long long int stateInd);

void statevec_initPureState(QubitRegister targetQureg, QubitRegister copyQureg);

void statevec_multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits);

void statevec_controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2);

void statevec_phaseShift(QubitRegister qureg, const int targetQubit, REAL angle);

void statevec_phaseShiftByTerm(QubitRegister qureg, const int targetQubit, Complex term);

void statevec_controlledPhaseShift(QubitRegister qureg, const int idQubit1, const int idQubit2, REAL angle); 

void statevec_multiControlledPhaseShift(QubitRegister qureg, int *controlQubits, int numControlQubits, REAL angle);

void statevec_sGate(QubitRegister qureg, const int targetQubit);

void statevec_tGate(QubitRegister qureg, const int targetQubit);

void statevec_sGateConj(QubitRegister qureg, const int targetQubit);

void statevec_tGateConj(QubitRegister qureg, const int targetQubit);

void statevec_sigmaX(QubitRegister qureg, const int targetQubit);

void statevec_sigmaY(QubitRegister qureg, const int targetQubit);

void statevec_sigmaYConj(QubitRegister qureg, const int targetQubit);

void statevec_sigmaZ(QubitRegister qureg, const int targetQubit);

void statevec_controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit);

void statevec_controlledSigmaYConj(QubitRegister qureg, const int controlQubit, const int targetQubit);

REAL statevec_getRealAmpEl(QubitRegister qureg, long long int index);

REAL statevec_getImagAmpEl(QubitRegister qureg, long long int index);

REAL statevec_getProbEl(QubitRegister qureg, long long int index);

REAL statevec_calcTotalProbability(QubitRegister qureg);

void statevec_compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta);

void statevec_unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u);

void statevec_rotateX(QubitRegister qureg, const int rotQubit, REAL angle);

void statevec_rotateY(QubitRegister qureg, const int rotQubit, REAL angle);

void statevec_rotateZ(QubitRegister qureg, const int rotQubit, REAL angle);

void statevec_rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis);

void statevec_rotateAroundAxisConj(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis);

void statevec_controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle);

void statevec_controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle);

void statevec_controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle);

void statevec_controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis);

void statevec_controlledRotateAroundAxisConj(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis);

void statevec_controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta);

void statevec_controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void statevec_multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u);

void statevec_hadamard(QubitRegister qureg, const int targetQubit);

void statevec_controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit);

REAL statevec_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome);

REAL statevec_collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome);

int statevec_measure(QubitRegister qureg, int measureQubit);

int statevec_measureWithStats(QubitRegister qureg, int measureQubit, REAL *stateProb);

# ifdef __cplusplus
}
# endif

# endif //QuEST_OPS