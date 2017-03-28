
# include <mpi.h>
# include <stdlib.h>
# include "qubits.h"
# include "qubits_mpi.h"

void initQUESTEnv(QUESTEnv *env){
        // init MPI environment
        int rank, numRanks;
        #if USE_MPI==1
                int provided;
                MPI_Init_thread(&narg, &varg, MPI_THREAD_FUNNELED, &provided);
                MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);

                if (DEBUG) {
                        char hostName[256];
                        int hostNameLen;
                        MPI_Get_processor_name(hostName, &hostNameLen);
                        printf("rank %d on host %s\n", rank, hostName);
                }
        #else
                rank=0; numRanks=1;
        #endif
	env->rank=rank;
	env->numRanks=numRanks;
}

void closeQUESTEnv(QUESTEnv env){
	MPI_Finalize();
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
        double allRankTotals=0;
	long long int index;
	long long int numAmpsPerRank = multiQubit.numAmps;
        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];      
                pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];      
        } 
	if (multiQubit.numChunks>1) MPI_Reduce(&pTotal, &allRankTotals, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return allRankTotals;
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
	// flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
	int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmps, rotQubit);
	Complex rot1, rot2;

	// rank's chunk is in upper half of block 
	int rankIsUpper;
	int pairRank; // rank of corresponding chunk

        // MPI send/receive vars
	int TAG=100;
        MPI_Status status;

	double *stateVecReal, stateVecImag, stateVecRealPair, stateVecImagPair;
	

	if (useLocalDataOnly){
		// all values required to update state vector lie in this rank
		rotateQubitLocal(multiQubit, rotQubit, alpha, beta);
	} else {
		// need to get corresponding chunk of state vector from other rank
		rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmps, rotQubit);
		getRotAngle(rankIsUpper, &rot1, &rot2, alpha, beta);
		pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmps, rotQubit);
		//printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
		// get corresponding values from my pair
		MPI_Sendrecv(multiQubit.stateVec.real, multiQubit.numAmps, MPI_DOUBLE, pairRank, TAG,
				 multiQubit.pairStateVec.real, multiQubit.numAmps, MPI_DOUBLE, pairRank, TAG,
				 MPI_COMM_WORLD, &status);
		//printf("rank: %d err: %d\n", multiQubit.rank, err);
		MPI_Sendrecv(multiQubit.stateVec.imag, multiQubit.numAmps, MPI_DOUBLE, pairRank, TAG,
				multiQubit.pairStateVec.imag, multiQubit.numAmps, MPI_DOUBLE, pairRank, TAG,
				MPI_COMM_WORLD, &status);
		// this rank's values are either in the upper of lower half of the block. send values to rotateQubitDistributed
		// in the correct order
		if (rankIsUpper){
			rotateQubitDistributed(multiQubit,rotQubit,rot1,rot2,
				multiQubit.stateVec.real,multiQubit.stateVec.imag,
				multiQubit.pairStateVec.real,multiQubit.pairStateVec.imag,
				multiQubit.stateVec.real,multiQubit.stateVec.imag);
		} else {
			rotateQubitDistributed(multiQubit,rotQubit,rot1,rot2,
				multiQubit.pairStateVec.real,multiQubit.pairStateVec.imag,
				multiQubit.stateVec.real,multiQubit.stateVec.imag,
				multiQubit.stateVec.real,multiQubit.stateVec.imag);
		}
	}
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
	double stateProb=0, totalStateProb=0;
	int skipValuesWithinRank = halfMatrixBlockFitsInChunk(multiQubit.numAmps, measureQubit);
	if (skipValuesWithinRank) {
		stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	} else {
		if (!isChunkToSkipInFindPZero(multiQubit.chunkId, multiQubit.numAmps, measureQubit)){
			stateProb = findProbabilityOfZeroDistributed(multiQubit, measureQubit);
		} else stateProb = 0;
	}
	MPI_Reduce(&stateProb, &totalStateProb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	return totalStateProb;
}


