// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * General functions used internally, supplied by QuEST_common or by hardware-specific backends.
 * Note that some bespoke functions used only internally exist in QuEST_qasm.h and QuEST_validation.h
 *
 * @author Ania Brown (statevecs, original architecture)
 * @author Tyson Jones (re-architecture, statevecs, density matrices)
 */

# ifndef QUEST_INTERNAL_H
# define QUEST_INTERNAL_H

# include "QuEST.h"
# include "QuEST_precision.h"

# ifdef __cplusplus
extern "C" {
# endif

    
/*
 * general functions
 */
  
long long int getQubitBitMask(int* controlQubits, int numControlQubits);

long long int getControlFlipMask(int* controlQubits, int* controlState, int numControlQubits);

unsigned long int hashString(char *str);

qreal getVectorMagnitude(Vector vec);

Complex getConjugateScalar(Complex scalar);

ComplexMatrix2 getConjugateMatrix2(ComplexMatrix2 src);

ComplexMatrix4 getConjugateMatrix4(ComplexMatrix4 src);

void setConjugateMatrixN(ComplexMatrixN m);

void ensureIndsIncrease(int* ind1, int* ind2);

void getComplexPairFromRotation(qreal angle, Vector axis, Complex* alpha, Complex* beta);

void getZYZRotAnglesFromComplexPair(Complex alpha, Complex beta, qreal* rz2, qreal* ry, qreal* rz1);

void getComplexPairAndPhaseFromUnitary(ComplexMatrix2 u, Complex* alpha, Complex* beta, qreal* globalPhase);

void shiftIndices(int* indices, int numIndices, int shift);

void shiftSubregIndices(int* allInds, int* numIndsPerReg, int numRegs, int shift);

void conjugateMatrixN(ComplexMatrixN u);

void getQuESTDefaultSeedKey(unsigned long int *key);


/*
 * operations upon density matrices 
 */

void densmatr_initPlusState(Qureg targetQureg);

void densmatr_initClassicalState(Qureg qureg, long long int stateInd);

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg);

qreal densmatr_calcTotalProb(Qureg qureg);

qreal densmatr_calcPurity(Qureg qureg);

qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState);

qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b);

qreal densmatr_calcInnerProduct(Qureg a, Qureg b);

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome);

void densmatr_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits);

void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb);
    
int densmatr_measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb);

void densmatr_mixDephasing(Qureg qureg, int targetQubit, qreal dephase);

void densmatr_mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal dephase);

void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depolLevel);

void densmatr_mixDamping(Qureg qureg, int targetQubit, qreal damping);

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel);

void densmatr_mixPauli(Qureg qureg, int qubit, qreal pX, qreal pY, qreal pZ);

void densmatr_mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg);

void densmatr_mixKrausMap(Qureg qureg, int target, ComplexMatrix2 *ops, int numOps);

void densmatr_mixTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps);

void densmatr_mixMultiQubitKrausMap(Qureg qureg, int* targets, int numTargets, ComplexMatrixN* ops, int numOps);

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op);

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op);

void densmatr_applyMultiQubitKrausSuperoperator(Qureg qureg, int *targets, int numTargets, ComplexMatrixN superOp); // exposed for QuESTlink Monte Carlo

void densmatr_setQuregToPauliHamil(Qureg qureg, PauliHamil hamil);


/* 
 * operations upon state vectors
 */
    
void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank);

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env);

void statevec_destroyQureg(Qureg qureg, QuESTEnv env);

void statevec_initBlankState(Qureg qureg);

void statevec_initZeroState(Qureg qureg);

void statevec_initPlusState(Qureg qureg);

void statevec_initDebugState(Qureg qureg);

void statevec_initClassicalState(Qureg qureg, long long int stateInd);

void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps);

void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg);

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits);

void statevec_controlledPhaseFlip(Qureg qureg, int idQubit1, int idQubit2);

void statevec_phaseShift(Qureg qureg, int targetQubit, qreal angle);

void statevec_phaseShiftByTerm(Qureg qureg, int targetQubit, Complex term);

void statevec_controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle); 

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle);

void statevec_sGate(Qureg qureg, int targetQubit);

void statevec_tGate(Qureg qureg, int targetQubit);

void statevec_sGateConj(Qureg qureg, int targetQubit);

void statevec_tGateConj(Qureg qureg, int targetQubit);

void statevec_pauliX(Qureg qureg, int targetQubit);

void statevec_pauliY(Qureg qureg, int targetQubit);

void statevec_pauliYConj(Qureg qureg, int targetQubit);

void statevec_pauliZ(Qureg qureg, int targetQubit);

void statevec_controlledPauliY(Qureg qureg, int controlQubit, int targetQubit);

void statevec_controlledPauliYConj(Qureg qureg, int controlQubit, int targetQubit);

qreal statevec_getRealAmp(Qureg qureg, long long int index);

qreal statevec_getImagAmp(Qureg qureg, long long int index);

qreal statevec_getProbAmp(Qureg qureg, long long int index);

qreal statevec_calcTotalProb(Qureg qureg);

qreal statevec_calcFidelity(Qureg qureg, Qureg pureState);

Complex statevec_calcInnerProduct(Qureg bra, Qureg ket);

qreal statevec_calcExpecPauliProd(Qureg qureg, int* targetQubits, enum pauliOpType* pauliCodes, int numTargets, Qureg workspace);

qreal statevec_calcExpecPauliSum(Qureg qureg, enum pauliOpType* allCodes, qreal* termCoeffs, int numSumTerms, Qureg workspace);

void statevec_compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta);

void statevec_unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u);

void statevec_twoQubitUnitary(Qureg qureg, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

void statevec_controlledTwoQubitUnitary(Qureg qureg, int controlQubit, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

void statevec_multiQubitUnitary(Qureg qureg, int* targets, int numTargets, ComplexMatrixN u);

void statevec_controlledMultiQubitUnitary(Qureg qureg, int ctrl, int* targets, int numTargets, ComplexMatrixN u);

void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u);

void statevec_rotateX(Qureg qureg, int rotQubit, qreal angle);

void statevec_rotateY(Qureg qureg, int rotQubit, qreal angle);

void statevec_rotateZ(Qureg qureg, int rotQubit, qreal angle);

void statevec_rotateAroundAxis(Qureg qureg, int rotQubit, qreal angle, Vector axis);

void statevec_rotateAroundAxisConj(Qureg qureg, int rotQubit, qreal angle, Vector axis);

void statevec_controlledRotateX(Qureg qureg, int controlQubit, int targetQubit, qreal angle);

void statevec_controlledRotateY(Qureg qureg, int controlQubit, int targetQubit, qreal angle);

void statevec_controlledRotateZ(Qureg qureg, int controlQubit, int targetQubit, qreal angle);

void statevec_controlledRotateAroundAxis(Qureg qureg, int controlQubit, int targetQubit, qreal angle, Vector axis);

void statevec_controlledRotateAroundAxisConj(Qureg qureg, int controlQubit, int targetQubit, qreal angle, Vector axis);

void statevec_controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta);

void statevec_controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u);

void statevec_multiControlledUnitary(Qureg qureg, long long int ctrlQubitsMask, long long int ctrlFlipMask, int targetQubit, ComplexMatrix2 u);

void statevec_hadamard(Qureg qureg, int targetQubit);

void statevec_controlledNot(Qureg qureg, int controlQubit, int targetQubit);

void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask);

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome);

void statevec_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits);

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb);

int statevec_measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb);

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2);

void statevec_sqrtSwapGate(Qureg qureg, int qb1, int qb2);

void statevec_sqrtSwapGateConj(Qureg qureg, int qb1, int qb2);

void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle);

void statevec_multiControlledMultiRotateZ(Qureg qureg, long long int ctrlMask, long long int targMask, qreal angle);

void statevec_multiRotatePauli(Qureg qureg, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle, int applyConj);

void statevec_multiControlledMultiRotatePauli(Qureg qureg, long long int ctrlMask, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle, int applyConj);

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out);

void statevec_applyPauliSum(Qureg inQureg, enum pauliOpType* allCodes, qreal* termCoeffs, int numSumTerms, Qureg outQureg);

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op);

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op);

void statevec_applyPhaseFuncOverrides(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, qreal* overridePhases, int numOverrides, int conj);

void statevec_applyMultiVarPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg, long long int* overrideInds, qreal* overridePhases, int numOverrides, int conj);

void statevec_applyParamNamedPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams, long long int* overrideInds, qreal* overridePhases, int numOverrides, int conj);

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps);

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps);

void statevec_applySubDiagonalOp(Qureg qureg, int* targets, SubDiagonalOp op, int conj);


/* 
 * operations which differentiate between state-vectors and density matrices internally 
 * (or do not operate upon Quregs)
 */
 
void agnostic_applyTrotterCircuit(Qureg qureg, PauliHamil hamil, qreal time, int order, int reps);

void agnostic_applyQFT(Qureg qureg, int* qubits, int numQubits);

DiagonalOp agnostic_createDiagonalOp(int numQubits, QuESTEnv env);

void agnostic_destroyDiagonalOp(DiagonalOp op);

void agnostic_syncDiagonalOp(DiagonalOp op);

void agnostic_setDiagonalOpElems(DiagonalOp op, long long int startInd, qreal* real, qreal* imag, long long int numElems);

void agnostic_initDiagonalOpFromPauliHamil(DiagonalOp op, PauliHamil hamil);

# ifdef __cplusplus
}
# endif

# endif // QUEST_INTERNAL_H
