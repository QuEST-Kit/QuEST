
# include <mpi.h>
# include "qubits_mpi.h"
# include "qubits.h"

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


double calcTotalProbability(int rank, long int numAmpsPerRank, int numRanks, double *stateVecReal, double *stateVecImag){
        double pTotal=0; 
        double *allRankTotals=0;
	int index;
        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=stateVecReal[index]*stateVecReal[index];      
                pTotal+=stateVecImag[index]*stateVecImag[index];      
        } 

       	if (rank==0) allRankTotals = malloc(numRanks*sizeof(allRankTotals)); 
        MPI_Gather(&pTotal, 1, MPI_DOUBLE, allRankTotals, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        if (rank==0){ 
                pTotal=0; 
                for (index=0; index<numRanks; index++){ 
                        pTotal+=allRankTotals[index]; 
                } 
        	free(allRankTotals); 
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

void rotateQubit(const long int numAmpsPerRank, const int numQubits, const int rotQubit,
		double aRe, double aIm, double bRe,  double bIm,
                double *restrict stateVecReal, double *restrict stateVecImag,
                double *restrict stateVecRealPair, double *restrict stateVecImagPair, int rank)

{
	// flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
	int useLocalDataOnly = halfMatrixBlockFitsInChunk(numAmpsPerRank, rotQubit);
	double rot1Re,rot1Im;
	double rot2Re,rot2Im;

	// rank's chunk is in upper half of block 
	int rankIsUpper;
	int pairRank; // rank of corresponding chunk

        // MPI send/receive vars
	int TAG=100;
        MPI_Status status;
	

	if (useLocalDataOnly){
		// all values required to update state vector lie in this rank
		rotateQubitLocal(numAmpsPerRank>>1,numQubits,rotQubit,aRe,aIm,bRe,bIm,stateVecReal,stateVecImag);
	} else {
		// need to get corresponding chunk of state vector from other rank
		rankIsUpper = chunkIsUpper(rank, numAmpsPerRank, rotQubit);
		getAlphaBeta(rankIsUpper, &rot1Re, &rot1Im, &rot2Re, &rot2Im, aRe, aIm, bRe, bIm);
		pairRank = getChunkPairId(rankIsUpper, rank, numAmpsPerRank, rotQubit);
		//printf("%d rank has pair rank: %d\n", rank, pairRank);
		// get corresponding values from my pair
		MPI_Sendrecv(stateVecReal, numAmpsPerRank, MPI_DOUBLE, pairRank, TAG,
				 stateVecRealPair, numAmpsPerRank, MPI_DOUBLE, pairRank, TAG,
				 MPI_COMM_WORLD, &status);
		//printf("rank: %d err: %d\n", rank, err);
		MPI_Sendrecv(stateVecImag, numAmpsPerRank, MPI_DOUBLE, pairRank, TAG,
				stateVecImagPair, numAmpsPerRank, MPI_DOUBLE, pairRank, TAG,
				MPI_COMM_WORLD, &status);
		// this rank's values are either in the upper of lower half of the block. send values to rotateQubitDistributed
		// in the correct order
		if (rankIsUpper){
			rotateQubitDistributed(numAmpsPerRank,numQubits,rotQubit,rot1Re,rot1Im,rot2Re,rot2Im,
				stateVecReal,stateVecImag,stateVecRealPair,stateVecImagPair,stateVecReal,stateVecImag);
		} else {
			rotateQubitDistributed(numAmpsPerRank,numQubits,rotQubit,rot1Re,rot1Im,rot2Re,rot2Im,
				stateVecRealPair,stateVecImagPair,stateVecReal,stateVecImag,stateVecReal,stateVecImag);
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


double findProbabilityOfZero(int rank, const long int numAmpsPerRank, const int numQubits,
                const int measureQubit,
                double *restrict stateVecReal,
                double *restrict stateVecImag)
{
	double stateProb=0, totalStateProb=0;
	int skipValuesWithinRank = halfMatrixBlockFitsInChunk(numAmpsPerRank, measureQubit);
	if (skipValuesWithinRank) {
		stateProb = findProbabilityOfZeroLocal(numAmpsPerRank>>1, numQubits, measureQubit, stateVecReal,stateVecImag);
	} else {
		if (!isChunkToSkipInFindPZero(rank, numAmpsPerRank, measureQubit)){
			stateProb = findProbabilityOfZeroDistributed(numAmpsPerRank, numQubits, measureQubit, stateVecReal,stateVecImag);
//			printf("measure %d, don't skip: %d. p %.8f\n", measureQubit, rank, stateProb);
		} else stateProb = 0;
	}
	MPI_Reduce(&stateProb, &totalStateProb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	return totalStateProb;
}


