
# include <mpi.h>
# include <stdlib.h>
# include "qubits.h"
# include "qubits_mpi.h"

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


double calcTotalProbability(Circuit circuit){
        double pTotal=0; 
        double allRankTotals=0;
	long long int index;
	long long int numAmpsPerRank = circuit.numAmps;
        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=circuit.stateVec.real[index]*circuit.stateVec.real[index];      
                pTotal+=circuit.stateVec.imag[index]*circuit.stateVec.imag[index];      
        } 
	if (circuit.numChunks>1) MPI_Reduce(&pTotal, &allRankTotals, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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

void rotateQubit(const int rotQubit, Complex alpha, Complex beta,
		Circuit *circuit)

{
	// flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
	int useLocalDataOnly = halfMatrixBlockFitsInChunk(circuit->numAmps, rotQubit);
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
		rotateQubitLocal(circuit, rotQubit, alpha, beta);
	} else {
		// need to get corresponding chunk of state vector from other rank
		rankIsUpper = chunkIsUpper(circuit->chunkId, circuit->numAmps, rotQubit);
		getRotAngle(rankIsUpper, &rot1, &rot2, alpha, beta);
		pairRank = getChunkPairId(rankIsUpper, circuit->chunkId, circuit->numAmps, rotQubit);
		//printf("%d rank has pair rank: %d\n", circuit->rank, pairRank);
		// get corresponding values from my pair
		MPI_Sendrecv(circuit->stateVec.real, circuit->numAmps, MPI_DOUBLE, pairRank, TAG,
				 circuit->pairStateVec.real, circuit->numAmps, MPI_DOUBLE, pairRank, TAG,
				 MPI_COMM_WORLD, &status);
		//printf("rank: %d err: %d\n", circuit->rank, err);
		MPI_Sendrecv(circuit->stateVec.imag, circuit->numAmps, MPI_DOUBLE, pairRank, TAG,
				circuit->pairStateVec.imag, circuit->numAmps, MPI_DOUBLE, pairRank, TAG,
				MPI_COMM_WORLD, &status);
		// this rank's values are either in the upper of lower half of the block. send values to rotateQubitDistributed
		// in the correct order
		if (rankIsUpper){
			rotateQubitDistributed(circuit,rotQubit,rot1,rot2,
				circuit->stateVec.real,circuit->stateVec.imag,
				circuit->pairStateVec.real,circuit->pairStateVec.imag,
				circuit->stateVec.real,circuit->stateVec.imag);
		} else {
			rotateQubitDistributed(circuit,rotQubit,rot1,rot2,
				circuit->pairStateVec.real,circuit->pairStateVec.imag,
				circuit->stateVec.real,circuit->stateVec.imag,
				circuit->stateVec.real,circuit->stateVec.imag);
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


double findProbabilityOfZero(Circuit *circuit,
                const int measureQubit)
{
	double stateProb=0, totalStateProb=0;
	int skipValuesWithinRank = halfMatrixBlockFitsInChunk(circuit->numAmps, measureQubit);
	if (skipValuesWithinRank) {
		stateProb = findProbabilityOfZeroLocal(circuit, measureQubit);
	} else {
		if (!isChunkToSkipInFindPZero(circuit->chunkId, circuit->numAmps, measureQubit)){
			stateProb = findProbabilityOfZeroDistributed(circuit, measureQubit);
		} else stateProb = 0;
	}
	MPI_Reduce(&stateProb, &totalStateProb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	return totalStateProb;
}


