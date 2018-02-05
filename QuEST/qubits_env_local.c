/** @file
An implementation of the API in qubits.h for a local (non-MPI) environment.
*/

# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "precision.h"
# include "qubits.h"
# include "qubits_internal.h"

void initQuESTEnv(QuESTEnv *env){
        // init MPI environment
	env->rank=0;
	env->numRanks=1;
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

REAL calcTotalProbability(MultiQubit multiQubit){
  /* IJB - implemented using Kahan summation for greater accuracy at a slight floating
     point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm */
  /* Don't change the bracketing in this routine! */
  REAL pTotal=0; 
  REAL y, t, c;
  long long int index;
  long long int numAmpsPerRank = multiQubit.numAmps;
  c = 0.0;
  for (index=0; index<numAmpsPerRank; index++){ 
    /* Perform pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index]; by Kahan */
   // pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];
    
    y = multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index] - c;
    t = pTotal + y;
    c = ( t - pTotal ) - y;
    pTotal = t;
    
    /* Perform pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index]; by Kahan */
    //pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];
    
    
    y = multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index] - c;
    t = pTotal + y;
    c = ( t - pTotal ) - y;
    pTotal = t;
    
    
  } 
  return pTotal;
}

REAL getRealAmpEl(MultiQubit multiQubit, long long int index){
	return multiQubit.stateVec.real[index];
}

REAL getImagAmpEl(MultiQubit multiQubit, long long int index){
	return multiQubit.stateVec.imag[index];
}

void compactUnitary(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

	// all values required to update state vector lie in this rank
	compactUnitaryLocal(multiQubit, rotQubit, alpha, beta);
}

void unitary(MultiQubit multiQubit, const int rotQubit, ComplexMatrix2 u) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

	// all values required to update state vector lie in this rank
	unitaryLocal(multiQubit, rotQubit, u);
}

void controlledCompactUnitary(MultiQubit multiQubit, const int rotQubit, const int controlQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != rotQubit, 3, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);
    

	controlledCompactUnitaryLocal(multiQubit, rotQubit, controlQubit, alpha, beta);
}

void controlledUnitary(MultiQubit multiQubit, const int rotQubit, const int controlQubit, ComplexMatrix2 u) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != rotQubit, 3, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);
   
	controlledUnitaryLocal(multiQubit, rotQubit, controlQubit, u);
}

void multiControlledUnitary(MultiQubit multiQubit, int* controlQubits, const int numControlQubits, const int rotQubit, ComplexMatrix2 u) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(numControlQubits >= 0 && numControlQubits < multiQubit.numQubits, 4, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    long long int mask=0; 
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<multiQubit.numQubits)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<rotQubit)) != (1LL<<rotQubit), 3, __func__);
	
    multiControlledUnitaryLocal(multiQubit, rotQubit, mask, u);
}

void sigmaX(MultiQubit multiQubit, const int rotQubit) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
	sigmaXLocal(multiQubit, rotQubit);
}

void sigmaY(MultiQubit multiQubit, const int rotQubit) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
	sigmaYLocal(multiQubit, rotQubit);
}

void phaseGate(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type)
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
	phaseGateLocal(multiQubit, rotQubit, type);
}

void hadamard(MultiQubit multiQubit, const int rotQubit) 
{
    QuESTAssert(rotQubit >= 0 && rotQubit < multiQubit.numQubits, 1, __func__);
	hadamardLocal(multiQubit, rotQubit);
}

void controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
	controlledNotLocal(multiQubit, controlQubit, targetQubit);
}

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
	REAL stateProb=0;
	stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	if (outcome==1) stateProb = 1.0 - stateProb;
	return stateProb;
}

REAL collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    REAL stateProb;
	stateProb = findProbabilityOfOutcome(multiQubit, measureQubit, outcome);
    QuESTAssert(fabs(stateProb>10e-13), 8, __func__);
    collapseToOutcomeLocal(multiQubit, measureQubit, stateProb, outcome);
    return stateProb;
}

void exitWithError(int errorCode, const char* func){
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, errorCodes[errorCode]);
    printf("!!!\n");
    printf("exiting..\n");
    exit(errorCode);
}

void QuESTAssert(int isValid, int errorCode, const char* func){
    if (!isValid) exitWithError(errorCode, func);
}
