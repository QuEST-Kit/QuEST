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

int syncQuESTSuccess(QuESTEnv env, int successCode){
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

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta) 
{
	// all values required to update state vector lie in this rank
	rotateQubitLocal(multiQubit, rotQubit, alpha, beta);
}

void controlRotateQubit(MultiQubit multiQubit, const int rotQubit, const int controlQubit, Complex alpha, Complex beta) 
{
	controlRotateQubitLocal(multiQubit, rotQubit, controlQubit, alpha, beta);
}

void sigmaX(MultiQubit multiQubit, const int rotQubit) 
{
	sigmaXLocal(multiQubit, rotQubit);
}

void sigmaY(MultiQubit multiQubit, const int rotQubit) 
{
	sigmaYLocal(multiQubit, rotQubit);
}

void phaseGate(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type)
{
	phaseGateLocal(multiQubit, rotQubit, type);
}

void hadamard(MultiQubit multiQubit, const int rotQubit) 
{
	hadamardLocal(multiQubit, rotQubit);
}

void controlNot(MultiQubit multiQubit, const int targetQubit, const int controlQubit) 
{
	controlNotLocal(multiQubit, targetQubit, controlQubit);
}

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
	REAL stateProb=0;
	stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	if (outcome==1) stateProb = 1.0 - stateProb;
	return stateProb;
}

REAL measureInState(MultiQubit multiQubit, const int measureQubit, int outcome)
{
        REAL stateProb;
	stateProb = findProbabilityOfOutcome(multiQubit, measureQubit, outcome);
        if (stateProb!=0) measureInStateLocal(multiQubit, measureQubit, stateProb, outcome);
        return stateProb;
}

REAL filterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
        REAL stateProb=0;
        stateProb = probOfFilterOut111(multiQubit, idQubit1, idQubit2, idQubit3);
        if (stateProb!=0) filterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3, stateProb);
        return stateProb;
}

REAL probOfFilterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
        REAL stateProb=0;
        stateProb = probOfFilterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3);
        return stateProb;
}



