// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * An implementation of QuEST's backend (../QuEST_internal.h) using NVIDIA's cuQuantum library
 *
 * @author Tyson Jones
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include <custatevec.h>



#ifdef __cplusplus
extern "C" {
#endif



void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps)
{
}


void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg)
{
}

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg)
{
}

void densmatr_initPlusState(Qureg qureg)
{
}

void densmatr_initClassicalState(Qureg qureg, long long int stateInd)
{
}

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env)
{   
}

void statevec_destroyQureg(Qureg qureg, QuESTEnv env)
{
}

void copyStateToGPU(Qureg qureg)
{
}

void copyStateFromGPU(Qureg qureg)
{
}

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
}

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
}

void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank)
{
}

qreal statevec_getRealAmp(Qureg qureg, long long int index)
{
    return -1;
}

qreal statevec_getImagAmp(Qureg qureg, long long int index)
{
    return -1;
}

void statevec_initBlankState(Qureg qureg)
{
}

void statevec_initZeroState(Qureg qureg)
{
}

void statevec_initPlusState(Qureg qureg)
{
}

void statevec_initClassicalState(Qureg qureg, long long int stateInd)
{
}

void statevec_initDebugState(Qureg qureg)
{
}

void statevec_initStateOfSingleQubit(Qureg *qureg, int qubitId, int outcome)
{
}

int statevec_initStateFromSingleFile(Qureg *qureg, char filename[200], QuESTEnv env)
{
    return -1;
}

int statevec_compareStates(Qureg mq1, Qureg mq2, qreal precision)
{
    return -1;
}

void statevec_compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta) 
{
}

void statevec_controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta) 
{
}

void statevec_unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u)
{
}

void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u)
{
}

void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int q1, int q2, ComplexMatrix4 u)
{
}

void statevec_controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u)
{
}

void statevec_multiControlledUnitary(
    Qureg qureg, 
    long long int ctrlQubitsMask, long long int ctrlFlipMask, 
    int targetQubit, ComplexMatrix2 u)
{
}

void statevec_pauliX(Qureg qureg, int targetQubit) 
{
}

void statevec_pauliY(Qureg qureg, int targetQubit) 
{
}

void statevec_pauliYConj(Qureg qureg, int targetQubit) 
{
}

void statevec_controlledPauliY(Qureg qureg, int controlQubit, int targetQubit)
{
}

void statevec_controlledPauliYConj(Qureg qureg, int controlQubit, int targetQubit)
{
}

void statevec_phaseShiftByTerm(Qureg qureg, int targetQubit, Complex term)
{   
}

void statevec_controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle)
{
}

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle)
{   
}

void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle)
{   
}

void statevec_multiControlledMultiRotateZ(Qureg qureg, long long int ctrlMask, long long int targMask, qreal angle)
{   
}

qreal densmatr_calcTotalProb(Qureg qureg)
{
    return -1;
}

qreal statevec_calcTotalProb(Qureg qureg)
{
    return -1;
}

void statevec_controlledPhaseFlip(Qureg qureg, int idQubit1, int idQubit2)
{
}

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits)
{
}

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) 
{
}

void statevec_hadamard(Qureg qureg, int targetQubit) 
{
}

void statevec_controlledNot(Qureg qureg, int controlQubit, int targetQubit)
{
}

void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask)
{
}

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    return -1;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    return -1;
}

void statevec_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits)
{
}

void densmatr_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits)
{
}

qreal densmatr_calcInnerProduct(Qureg a, Qureg b)
{
    return -1;
}

Complex statevec_calcInnerProduct(Qureg bra, Qureg ket)
{
    return (Complex) {.real=-1, .imag=-1};
}

qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState)
{
    return -1;
}

qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b)
{
    return -1;
}

qreal densmatr_calcPurity(Qureg qureg)
{
    return -1;
}

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{        
}

void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{
}

void densmatr_mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg)
{
}

void densmatr_mixDephasing(Qureg qureg, int targetQubit, qreal dephase) 
{
}

void densmatr_mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal dephase)
{
}

void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depolLevel)
{
}

void densmatr_mixDamping(Qureg qureg, int targetQubit, qreal damping)
{
}

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel)
{
}

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out)
{
}

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) 
{
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op) {
}

void statevec_applySubDiagonalOp(Qureg qureg, int* targets, SubDiagonalOp op, int conj)
{
}

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    return (Complex) {.real=-1, .imag=-1};
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    return (Complex) {.real=-1, .imag=-1};
}

void statevec_applyPhaseFuncOverrides(
    Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int numTerms, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
}

void statevec_applyMultiVarPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int* numTermsPerReg, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
}

void statevec_applyParamNamedPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    enum phaseFunc phaseFuncName, qreal* params, int numParams,
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
}

void densmatr_setQuregToPauliHamil(Qureg qureg, PauliHamil hamil)
{
}



#ifdef __cplusplus
}
#endif
