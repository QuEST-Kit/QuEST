/** @file
An implementation of the API in qubits_env_wrapper.h for a local (non-MPI) environment.
*/

# include <stdlib.h>
# include "qubits.h"
# include "qubits_env_wrapper.h"

void initQUESTEnv(QUESTEnv *env){
        // init MPI environment
        int rank, numRanks;
	env->rank=0;
	env->numRanks=1;
}

void syncQUESTEnv(QUESTEnv env){
	// MPI Barrier goes here in MPI version. 
} 

void closeQUESTEnv(QUESTEnv env){
	// MPI finalize goes here in MPI version. Call this function anyway for consistency
}

double calcTotalProbability(MultiQubit multiQubit){
        double pTotal=0; 
	long long int index;
	long long int numAmpsPerRank = multiQubit.numAmps;
        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];      
                pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];      
        } 
	return pTotal;
}

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta)

{
	// all values required to update state vector lie in this rank
	rotateQubitLocal(multiQubit, rotQubit, alpha, beta);
}

double findProbabilityOfZero(MultiQubit multiQubit,
                const int measureQubit)
{
	double stateProb=0;
	stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	return stateProb;
}


