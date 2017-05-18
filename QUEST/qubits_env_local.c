/** @file
An implementation of the API in qubits.h for a local (non-MPI) environment.
*/

# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "precision.h"
# include "qubits.h"
# include "qubits_internal.h"

void initQUESTEnv(QUESTEnv *env){
        // init MPI environment
	env->rank=0;
	env->numRanks=1;
}

void syncQUESTEnv(QUESTEnv env){
	// MPI Barrier goes here in MPI version. 
} 

void closeQUESTEnv(QUESTEnv env){
	// MPI finalize goes here in MPI version. Call this function anyway for consistency
}

void reportQUESTEnv(QUESTEnv env){
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

REAL calcTotalProbability(MultiQubit multiQubit){
        REAL pTotal=0; 
	long long int index;
	long long int numAmpsPerRank = multiQubit.numAmps;
        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];      
                pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];      
        } 
	return pTotal;
}

REAL getRealAmpEl(MultiQubit multiQubit, long long int index){
	return multiQubit.stateVec.real[index];
}

REAL getImagAmpEl(MultiQubit multiQubit, long long int index){
	return multiQubit.stateVec.imag[index];
}

REAL getProbEl(MultiQubit multiQubit, long long int index){
	REAL real;
	REAL imag;
	real = multiQubit.stateVec.real[index];
	imag = multiQubit.stateVec.imag[index];
	return real*real + imag*imag;
}

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta) 
{
	// all values required to update state vector lie in this rank
	rotateQubitLocal(multiQubit, rotQubit, alpha, beta);
}

REAL findProbabilityOfZero(MultiQubit multiQubit,
                const int measureQubit)
{
	REAL stateProb=0;
	stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	return stateProb;
}

REAL measureInZero(MultiQubit multiQubit, const int measureQubit)
{
        REAL stateProb;
	stateProb = findProbabilityOfZero(multiQubit, measureQubit);
        measureInZeroLocal(multiQubit, measureQubit, stateProb);
        return stateProb;
}

REAL filterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
        REAL stateProb=0;
        stateProb = probOfFilterOut111(multiQubit, idQubit1, idQubit2, idQubit3);
        filterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3, stateProb);
        return stateProb;
}

REAL probOfFilterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
        REAL stateProb=0;
        stateProb = probOfFilterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3);
        return stateProb;
}



