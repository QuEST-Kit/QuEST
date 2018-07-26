// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the pure backend in ../QuEST_ops_pure.h for a local (non-MPI, non-GPU) environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 */

# include "../QuEST.h"
# include "../QuEST_internal.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"

# include "QuEST_cpu_internal.h"

// debug: remove this after all validation is removed
# include "../QuEST_validation.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <sys/types.h>

# ifdef _OPENMP
# include <omp.h>
# endif




// @TODO
void densmatr_initPureState(QubitRegister targetQureg, QubitRegister copyQureg) {
	
	QuESTAssert(targetQureg.isDensityMatrix, 14, __func__);
	QuESTAssert( !copyQureg.isDensityMatrix, 12, __func__);
	QuESTAssert(targetQureg.numQubitsRepresented==copyQureg.numQubitsInStateVec, 13, __func__);
	
	densmatr_initPureStateLocal(targetQureg, copyQureg);
}

REAL densmatr_calcTotalProbability(QubitRegister qureg) {
	
	// computes the trace using Kahan summation
	REAL pTotal=0;
	REAL y, t, c;
	c = 0;
	
	long long int numCols = 1LL << qureg.numQubitsRepresented;
	long long diagIndex;
	
	for (int col=0; col< numCols; col++) {
		diagIndex = col*(numCols + 1);
		y = qureg.stateVec.real[diagIndex] - c;
		t = pTotal + y;
		c = ( t - pTotal ) - y; // brackets are important
		pTotal = t;
	}
	
	// @TODO should maybe do a cheap test that imaginary components are ~0
	
	return pTotal;
}

REAL statevec_calcTotalProbability(QubitRegister qureg){
    // implemented using Kahan summation for greater accuracy at a slight floating
    // point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    REAL pTotal=0; 
    REAL y, t, c;
    long long int index;
    long long int numAmpsPerRank = qureg.numAmpsPerChunk;
    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){ 
        // Perform pTotal+=qureg.stateVec.real[index]*qureg.stateVec.real[index]; by Kahan

        y = qureg.stateVec.real[index]*qureg.stateVec.real[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;

        // Perform pTotal+=qureg.stateVec.imag[index]*qureg.stateVec.imag[index]; by Kahan

        y = qureg.stateVec.imag[index]*qureg.stateVec.imag[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;


    } 
    return pTotal;
}


void initQuESTEnv(QuESTEnv *env){
    // init MPI environment
    env->rank=0;
    env->numRanks=1;
	
	seedQuESTDefault();
}

void syncQuESTEnv(QuESTEnv env){
    // MPI Barrier goes here in MPI version. 
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void closeQuESTEnv(QuESTEnv env){
    // MPI finalize goes here in MPI version. Call this function anyway for consistency
}

void reportQuESTEnv(QuESTEnv env){
    printf("EXECUTION ENVIRONMENT:\n");
    printf("Running locally on one node\n");
    printf("Number of ranks is %d\n", env.numRanks);
# ifdef _OPENMP
    printf("OpenMP enabled\n");
    printf("Number of threads available is %d\n", omp_get_max_threads());
# else
    printf("OpenMP disabled\n");
# endif
    printf("Precision: size of REAL is %ld bytes\n", sizeof(REAL));
}

void reportNodeList(QuESTEnv env){
    printf("Hostname unknown: running locally\n");
}

REAL statevec_getRealAmpEl(QubitRegister qureg, long long int index){
    return qureg.stateVec.real[index];
}

REAL statevec_getImagAmpEl(QubitRegister qureg, long long int index){
    return qureg.stateVec.imag[index];
}

void statevec_compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(isComplexPairUnitary(alpha, beta), 6, __func__);

    // all values required to update state vector lie in this rank
    statevec_compactUnitaryLocal(qureg, targetQubit, alpha, beta);
}

void statevec_unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(isMatrixUnitary(u), 5, __func__);

    // all values required to update state vector lie in this rank
    statevec_unitaryLocal(qureg, targetQubit, u);
}

void statevec_controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(isComplexPairUnitary(alpha, beta), 6, __func__);


    statevec_controlledCompactUnitaryLocal(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(isMatrixUnitary(u), 5, __func__);

    statevec_controlledUnitaryLocal(qureg, controlQubit, targetQubit, u);
}

void statevec_multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(numControlQubits > 0 && numControlQubits <= qureg.numQubitsInStateVec, 4, __func__);
    QuESTAssert(isMatrixUnitary(u), 5, __func__);

    long long int mask=0; 
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<qureg.numQubitsInStateVec)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<targetQubit)) != (1LL<<targetQubit), 3, __func__);

    statevec_multiControlledUnitaryLocal(qureg, targetQubit, mask, u);
}

void statevec_sigmaX(QubitRegister qureg, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    statevec_sigmaXLocal(qureg, targetQubit);
}

void statevec_sigmaY(QubitRegister qureg, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
	int conjFac = 1;
    statevec_sigmaYLocal(qureg, targetQubit, conjFac);
}

void statevec_sigmaYConj(QubitRegister qureg, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
	int conjFac = -1;
    statevec_sigmaYLocal(qureg, targetQubit, conjFac);
}

void statevec_controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
	int conjFac = 1;
	statevec_controlledSigmaYLocal(qureg, controlQubit, targetQubit, conjFac);
}

void statevec_controlledSigmaYConj(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
	int conjFac = -1;
	statevec_controlledSigmaYLocal(qureg, controlQubit, targetQubit, conjFac);
}

void statevec_hadamard(QubitRegister qureg, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    statevec_hadamardLocal(qureg, targetQubit);
}

void statevec_controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    statevec_controlledNotLocal(qureg, controlQubit, targetQubit);
}

REAL statevec_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubitsInStateVec, 2, __func__);
    REAL stateProb=0;
    stateProb = statevec_findProbabilityOfZeroLocal(qureg, measureQubit);
    if (outcome==1) stateProb = 1.0 - stateProb;
    return stateProb;
}

REAL densmatr_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	
	REAL outcomeProb = densmatr_findProbabilityOfZeroLocal(qureg, measureQubit);
	if (outcome == 1)
		outcomeProb = 1.0 - outcomeProb;
	return outcomeProb;
}

REAL statevec_collapseToKnownProbOutcome(QubitRegister qureg, const int measureQubit, int outcome, REAL stateProb)
{
    statevec_collapseToOutcomeLocal(qureg, measureQubit, outcome, stateProb);
    return stateProb;
}
