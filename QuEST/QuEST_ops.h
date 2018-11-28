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

/* operations upon density matrices */

void densmatr_initPlusState(QubitRegister targetQureg);

void densmatr_initClassicalState(QubitRegister qureg, long long int stateInd);

void densmatr_initPureState(QubitRegister targetQureg, QubitRegister copyQureg);

REAL densmatr_calcTotalProb(QubitRegister qureg);

REAL densmatr_calcPurity(QubitRegister qureg);

REAL densmatr_calcFidelity(QubitRegister qureg, QubitRegister pureState);

REAL densmatr_calcProbOfOutcome(QubitRegister qureg, const int measureQubit, int outcome);

void densmatr_collapseToKnownProbOutcome(QubitRegister qureg, const int measureQubit, int outcome, REAL outcomeProb);
    
int densmatr_measureWithStats(QubitRegister qureg, int measureQubit, REAL *outcomeProb);

void densmatr_oneQubitDephase(QubitRegister qureg, const int targetQubit, REAL dephase);

void densmatr_twoQubitDephase(QubitRegister qureg, const int qubit1, const int qubit2, REAL dephase);

void densmatr_oneQubitDepolarise(QubitRegister qureg, const int targetQubit, REAL depolLevel);

void densmatr_twoQubitDepolarise(QubitRegister qureg, int qubit1, int qubit2, REAL depolLevel);

void densmatr_addDensityMatrix(QubitRegister combineQureg, REAL otherProb, QubitRegister otherQureg);
    
/* operations upon state vectors */
    
void statevec_reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank);

int statevec_compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision);

int statevec_initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env);

void statevec_initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome);

void statevec_createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env);

void statevec_destroyQubitRegister(QubitRegister qureg, QuESTEnv env);

void statevec_initZeroState(QubitRegister qureg);

void statevec_initPlusState(QubitRegister qureg);

void statevec_initStateDebug(QubitRegister qureg);

void statevec_initClassicalState(QubitRegister qureg, long long int stateInd);

void statevec_initStateFromAmps(QubitRegister qureg, long long int startInd, REAL* reals, REAL* imags, long long int numAmps);

void statevec_cloneQubitRegister(QubitRegister targetQureg, QubitRegister copyQureg);

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

REAL statevec_calcTotalProb(QubitRegister qureg);

REAL statevec_calcFidelity(QubitRegister qureg, QubitRegister pureState);

Complex statevec_calcInnerProduct(QubitRegister bra, QubitRegister ket);

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

REAL statevec_calcProbOfOutcome(QubitRegister qureg, const int measureQubit, int outcome);

void statevec_collapseToKnownProbOutcome(QubitRegister qureg, const int measureQubit, int outcome, REAL outcomeProb);

int statevec_measureWithStats(QubitRegister qureg, int measureQubit, REAL *outcomeProb);



# ifdef __cplusplus
}
# endif

# endif //QuEST_OPS