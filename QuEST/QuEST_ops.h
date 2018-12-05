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

void densmatr_initPlusState(Qureg targetQureg);

void densmatr_initClassicalState(Qureg qureg, long long int stateInd);

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg);

REAL densmatr_calcTotalProb(Qureg qureg);

REAL densmatr_calcPurity(Qureg qureg);

REAL densmatr_calcFidelity(Qureg qureg, Qureg pureState);

REAL densmatr_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome);

void densmatr_collapseToKnownProbOutcome(Qureg qureg, const int measureQubit, int outcome, REAL outcomeProb);
    
int densmatr_measureWithStats(Qureg qureg, int measureQubit, REAL *outcomeProb);

void densmatr_oneQubitDephase(Qureg qureg, const int targetQubit, REAL dephase);

void densmatr_twoQubitDephase(Qureg qureg, const int qubit1, const int qubit2, REAL dephase);

void densmatr_oneQubitDepolarise(Qureg qureg, const int targetQubit, REAL depolLevel);

void densmatr_twoQubitDepolarise(Qureg qureg, int qubit1, int qubit2, REAL depolLevel);

void densmatr_addDensityMatrix(Qureg combineQureg, REAL otherProb, Qureg otherQureg);
    
/* operations upon state vectors */
    
void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank);

int statevec_compareStates(Qureg mq1, Qureg mq2, REAL precision);

int statevec_initStateFromSingleFile(Qureg *qureg, char filename[200], QuESTEnv env);

void statevec_initStateOfSingleQubit(Qureg *qureg, int qubitId, int outcome);

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env);

void statevec_destroyQureg(Qureg qureg, QuESTEnv env);

void statevec_initZeroState(Qureg qureg);

void statevec_initPlusState(Qureg qureg);

void statevec_initStateDebug(Qureg qureg);

void statevec_initClassicalState(Qureg qureg, long long int stateInd);

void statevec_setAmps(Qureg qureg, long long int startInd, REAL* reals, REAL* imags, long long int numAmps);

void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg);

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits);

void statevec_controlledPhaseFlip(Qureg qureg, const int idQubit1, const int idQubit2);

void statevec_phaseShift(Qureg qureg, const int targetQubit, REAL angle);

void statevec_phaseShiftByTerm(Qureg qureg, const int targetQubit, Complex term);

void statevec_controlledPhaseShift(Qureg qureg, const int idQubit1, const int idQubit2, REAL angle); 

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, REAL angle);

void statevec_sGate(Qureg qureg, const int targetQubit);

void statevec_tGate(Qureg qureg, const int targetQubit);

void statevec_sGateConj(Qureg qureg, const int targetQubit);

void statevec_tGateConj(Qureg qureg, const int targetQubit);

void statevec_pauliX(Qureg qureg, const int targetQubit);

void statevec_pauliY(Qureg qureg, const int targetQubit);

void statevec_pauliYConj(Qureg qureg, const int targetQubit);

void statevec_pauliZ(Qureg qureg, const int targetQubit);

void statevec_controlledPauliY(Qureg qureg, const int controlQubit, const int targetQubit);

void statevec_controlledPauliYConj(Qureg qureg, const int controlQubit, const int targetQubit);

REAL statevec_getRealAmp(Qureg qureg, long long int index);

REAL statevec_getImagAmp(Qureg qureg, long long int index);

REAL statevec_getProbAmp(Qureg qureg, long long int index);

REAL statevec_calcTotalProb(Qureg qureg);

REAL statevec_calcFidelity(Qureg qureg, Qureg pureState);

Complex statevec_calcInnerProduct(Qureg bra, Qureg ket);

void statevec_compactUnitary(Qureg qureg, const int targetQubit, Complex alpha, Complex beta);

void statevec_unitary(Qureg qureg, const int targetQubit, ComplexMatrix2 u);

void statevec_rotateX(Qureg qureg, const int rotQubit, REAL angle);

void statevec_rotateY(Qureg qureg, const int rotQubit, REAL angle);

void statevec_rotateZ(Qureg qureg, const int rotQubit, REAL angle);

void statevec_rotateAroundAxis(Qureg qureg, const int rotQubit, REAL angle, Vector axis);

void statevec_rotateAroundAxisConj(Qureg qureg, const int rotQubit, REAL angle, Vector axis);

void statevec_controlledRotateX(Qureg qureg, const int controlQubit, const int targetQubit, REAL angle);

void statevec_controlledRotateY(Qureg qureg, const int controlQubit, const int targetQubit, REAL angle);

void statevec_controlledRotateZ(Qureg qureg, const int controlQubit, const int targetQubit, REAL angle);

void statevec_controlledRotateAroundAxis(Qureg qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis);

void statevec_controlledRotateAroundAxisConj(Qureg qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis);

void statevec_controlledCompactUnitary(Qureg qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta);

void statevec_controlledUnitary(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void statevec_multiControlledUnitary(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u);

void statevec_hadamard(Qureg qureg, const int targetQubit);

void statevec_controlledNot(Qureg qureg, const int controlQubit, const int targetQubit);

REAL statevec_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome);

void statevec_collapseToKnownProbOutcome(Qureg qureg, const int measureQubit, int outcome, REAL outcomeProb);

int statevec_measureWithStats(Qureg qureg, int measureQubit, REAL *outcomeProb);



# ifdef __cplusplus
}
# endif

# endif //QuEST_OPS