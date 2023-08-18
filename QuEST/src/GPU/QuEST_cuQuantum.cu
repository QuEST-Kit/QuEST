// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * An implementation of QuEST's backend (../QuEST_internal.h) using NVIDIA's cuQuantum library.
 * This makes no use of the ComplexArray qureg.deviceStateVec, used by the bespoke GPU kernels,
 * which is not malloc'd in this deployment. Instead, this cuQuantum backend mallocs and uses
 * two dedicated arrays of 'cuAmp' complex primitives; qureg.cuStateVec (CPU memory) and
 * qureg.deviceCuStateVec (GPU memory)
 *
 * @author Tyson Jones
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include <custatevec.h>



// precision-overloaded macros for creating a cuAmp
# if QuEST_PREC==1
    # define toCuAmp(re, im) make_cuFloatComplex(re, im)
    # define cuAmpReal(amp) cuCrealf(amp)
    # define cuAmpImag(amp) cuCimagf(amp)
# elif QuEST_PREC==2
    # define toCuAmp(re, im) make_cuDoubleComplex(re, im)
    # define cuAmpReal(amp) cuCreal(amp)
    # define cuAmpImag(amp) cuCimag(amp)
# elif QuEST_PREC==4
    # define toCuAmp(re, im) -1 // invalid precision config
    # define cuAmpReal(amp) -1
    # define cuAmpImag(amp) -1
#endif



#ifdef __cplusplus
extern "C" {
#endif



/* 
 * QUREG CREATION AND AMP SET/GET
 */

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env)
{   
    // set fields
    long long int numAmps = 1LL << numQubits;
    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsPerChunk = numAmps;
    qureg->numAmpsTotal = numAmps;
    qureg->chunkId = 0;
    qureg->numChunks = 1;
    qureg->isDensityMatrix = 0;

    // allocate user-facing CPU memory
    qureg->stateVec.real = (qreal*) malloc(numAmps * sizeof(qureg->stateVec.real));
    qureg->stateVec.imag = (qreal*) malloc(numAmps * sizeof(qureg->stateVec.imag));
    validateQuregAllocation(qureg, env, __func__);

    // allocate cuQuantum GPU memory (unvalidated)
    cudaMalloc( &(qureg->deviceCuStateVec), numAmps * sizeof(*(qureg->deviceCuStateVec)) );

    // allocate private cuQuantum CPU memory (for exchanging with GPU memory)
    qureg->cuStateVec = (cuAmp*) malloc(numAmps * sizeof(qureg->cuStateVec));
}

void statevec_destroyQureg(Qureg qureg, QuESTEnv env)
{
    // free user-facing CPU memory 
    free(qureg.stateVec.real);
    free(qureg.stateVec.imag);

    // free private cuQuantum CPU memory
    free(qureg.cuStateVec);

    // free cuQuantum GPU memory
    cudaFree(qureg.deviceCuStateVec);
}

void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps)
{
    // slowly manually overwrite subset of private cuQuantum CPU memory
    for (long long int i=0; i<numAmps; i++)
        qureg.cuStateVec[i+startInd] = toCuAmp(reals[i], imags[i]);

    // cuda-copy subset to GPU memory subset
    cudaMemcpy(
        &(qureg.deviceCuStateVec[startInd]), 
        &(qureg.cuStateVec[startInd]), 
        numAmps * sizeof(cuAmp), cudaMemcpyHostToDevice);
}

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    statevec_setAmps(qureg, startInd, &(qureg.stateVec.real[startInd]), &(qureg.stateVec.imag[startInd]), numAmps);
}

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    // cuda-copy subset of GPU memory to private cuQuantum CPU memory
    cudaMemcpy(
        &(qureg.cuStateVec[startInd]), 
        &(qureg.deviceCuStateVec[startInd]), 
        numAmps * sizeof(cuAmp), cudaMemcpyDeviceToHost);

    // slowly manually overwrite public CPU memory from private
    for (long long int i=startInd; i<(startInd+numAmps); i++) {
        qureg.stateVec.real[i] = cuAmpReal(qureg.cuStateVec[i]);
        qureg.stateVec.imag[i] = cuAmpImag(qureg.cuStateVec[i]);
    }
}

void copyStateToGPU(Qureg qureg)
{
    statevec_copySubstateToGPU(qureg, 0, qureg.numAmpsTotal);
}

void copyStateFromGPU(Qureg qureg)
{
    statevec_copySubstateFromGPU(qureg, 0, qureg.numAmpsTotal);
}

void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg)
{
    // directly cuda-copy the GPU memory 
    cudaMemcpy(
        targetQureg.deviceCuStateVec,
        copyQureg.deviceCuStateVec,
        copyQureg.numAmpsTotal * sizeof(cuAmp),
        cudaMemcpyDeviceToDevice);
}

qreal statevec_getRealAmp(Qureg qureg, long long int index)
{
    cuAmp amp;
    cudaMemcpy(&amp, &(qureg.deviceCuStateVec[index]), sizeof(cuAmp), cudaMemcpyDeviceToHost);
    return cuAmpReal(amp);
}

qreal statevec_getImagAmp(Qureg qureg, long long int index)
{
    cuAmp amp;
    cudaMemcpy(&amp, &(qureg.deviceCuStateVec[index]), sizeof(cuAmp), cudaMemcpyDeviceToHost);
    return cuAmpImag(amp);
}



/*
 * STATE INITIALISATION
 */

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg)
{
}

void densmatr_initPlusState(Qureg qureg)
{
}

void densmatr_initClassicalState(Qureg qureg, long long int stateInd)
{
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

void densmatr_setQuregToPauliHamil(Qureg qureg, PauliHamil hamil)
{
}

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out)
{
}



/*
 * DEBUG
 */

void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank)
{
}

int statevec_compareStates(Qureg mq1, Qureg mq2, qreal precision)
{
    return -1;
}



/*
 * OPERATORS
 */

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

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) 
{
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op)
{
}

void statevec_applySubDiagonalOp(Qureg qureg, int* targets, SubDiagonalOp op, int conj)
{
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



/*
 * DECOHERENCE
 */

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



/*
 * CALCULATIONS
 */

qreal densmatr_calcTotalProb(Qureg qureg)
{
    return -1;
}

qreal statevec_calcTotalProb(Qureg qureg)
{
    return -1;
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

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    return (Complex) {.real=-1, .imag=-1};
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    return (Complex) {.real=-1, .imag=-1};
}



/*
 * REDUCTIONS
 */

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{        
}

void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{
}



#ifdef __cplusplus
}
#endif
