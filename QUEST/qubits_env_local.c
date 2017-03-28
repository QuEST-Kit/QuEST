
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

// ==================================================================== //
//                                                                      //
//     calcTotalProbability -- calculate total probability over all 	//
//     			ranks by taking the norm of the state vector.	//
//     			Should be equal to one.				// 
//                                                                      //
//     input:                                                           //
//                    rank -- mpi proc id				// 
//                    numAmpsPerRank  -- number of amps on one rank	//
//                    numRanks -- num mpi processes			//
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector                 //
//                                                                      //
//     output:                                                          //
//     		      double -- total probability			//
//                                                                      //
// ==================================================================== //


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

// ==================================================================== 
// rotateQubit: 
// inputs:
// 	numAmpsPerRank -- number of amps in mpi proc
//      numQubits -- number of qubits                     
// 	rotQubit -- qubit being rotated 
//      aRe,    -- real/imag part of                
//      aIm        rotation angle alpha             
//      bRe,     -- real/imag part of                
//      bIm         rotation angle beta              
//
// temp:
//      stateVecRealPair, -- real/imag parts of the state vector offset by half a block
//      stateVecImagPair     from the state vector updated on this rank, used to update
//      		     this rank's state vector
//
// outputs:
//      stateVecReal, -- real/imag parts of               
//      stateVecImag     the state vector updated on this rank                 
// ==================================================================== 

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta)

{
	// all values required to update state vector lie in this rank
	rotateQubitLocal(multiQubit, rotQubit, alpha, beta);
}

// ==================================================================== 
// findProbabilityOfZero: 
// inputs:
//      rank -- mpi proc id				 
// 	numAmpsPerRank -- number of amps in mpi proc
//      numQubits -- number of qubits                     
// 	measureQubit -- qubit being measured
//      stateVecReal, -- real/imag parts of               
//      stateVecImag     the state vector updated on this rank                 
// ==================================================================== 


double findProbabilityOfZero(MultiQubit multiQubit,
                const int measureQubit)
{
	double stateProb=0;
	stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	return stateProb;
}


