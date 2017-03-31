/** @file
An implementation of the API in qubits_env_wrapper.h for an MPI environment.
*/
# include <mpi.h>
# include <stdlib.h>
# include "qubits.h"
# include "qubits_env_wrapper.h"

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

void syncQUESTEnv(QUESTEnv env){
	MPI_Barrier(MPI_COMM_WORLD);
}

void closeQUESTEnv(QUESTEnv env){
	MPI_Finalize();
}

/** Find chunks to skip when calculating probability of qubit being zero.
When calculating probability of a bit q being zero,
sum up 2^q values, then skip 2^q values, etc. This function finds if an entire chunk
is in the range of values to be skipped

@param[in] chunkId id of chunk in state vector
@param[in] chunkSize number of amps in chunk
@param[in] measureQubi qubit being measured
@return int -- 1: skip, 0: don't skip
*/

int isChunkToSkipInFindPZero(int chunkId, int chunkSize, int measureQubit){
        long long int sizeHalfBlock = 1LL << (measureQubit);
        int numChunksToSkip = sizeHalfBlock/chunkSize;
        // calculate probability by summing over numChunksToSkip, then skipping numChunksToSkip, etc
        int bitToCheck = chunkId & numChunksToSkip;
        return bitToCheck;
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

/** Returns whether a given chunk in position chunkId is in the upper or lower half of
a block.

@param[in] chunkId id of chunk in state vector
@param[in] chunkSize number of amps in chunk
@param[in] rotQubit qubit being rotated 
@return 1: chunk is in upper half of block, 0: chunk is in lower half of block 
*/

int chunkIsUpper(int chunkId, int chunkSize, int rotQubit)
{       
        long long int sizeHalfBlock = 1LL << (rotQubit);
        long long int sizeBlock = sizeHalfBlock*2;
        long long int posInBlock = (chunkId*chunkSize) % sizeBlock;
        return posInBlock<sizeHalfBlock;
}

/** Get rotation values for a given chunk
@param[in] chunkIsUpper 1: chunk is in upper half of block, 0: chunk is in lower half

@param[out] rot1, rot2 rotation values to use, allocated for upper/lower such that
@verbatim
stateUpper = rot1 * stateUpper + conj(rot2)  * stateLower
@endverbatim
or
@verbatim
stateLower = rot1 * stateUpper + conj(rot2)  * stateLower
@endverbatim
@param[in] alpha, beta initial rotation values 
*/
void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta)
{
        if (chunkIsUpper){
                *rot1=alpha;
                rot2->real=-beta.real;
                rot2->imag=-beta.imag;
        } else {
                *rot1=beta;
                *rot2=alpha;
        }
}
/** get position of corresponding chunk, holding values required to
update values in my chunk (with chunkId) when rotating rotQubit.

@param[in] chunkIsUpper 1: chunk is in upper half of block, 0: chunk is in lower half
@param[in] chunkId id of chunk in state vector
@param[in] chunkSize number of amps in chunk
@param[in] rotQubit qubit being rotated 
@return chunkId of chunk required to rotate rotQubit 
*/

int getChunkPairId(int chunkIsUpper, int chunkId, int chunkSize, int rotQubit)
{
        long long int sizeHalfBlock = 1LL << (rotQubit);
        int chunksPerHalfBlock = sizeHalfBlock/chunkSize;
        if (chunkIsUpper){
                return chunkId + chunksPerHalfBlock;
        } else {
                return chunkId - chunksPerHalfBlock;
        }
}

/** return whether the current qubit rotation will use
blocks that fit within a single chunk.

@param[in] chunkSize number of amps in chunk
@param[in] rotQubit qubit being rotated 
@return 1: one chunk fits in one block 0: chunk is larger than block
*/

int halfMatrixBlockFitsInChunk(int chunkSize, int rotQubit)
{
        long long int sizeHalfBlock = 1LL << (rotQubit);
        if (chunkSize > sizeHalfBlock) return 1;
        else return 0;
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
				multiQubit.stateVec, //upper
				multiQubit.pairStateVec, //lower
				multiQubit.stateVec); //output
		} else {
			rotateQubitDistributed(multiQubit,rotQubit,rot1,rot2,
				multiQubit.pairStateVec, //upper
				multiQubit.stateVec, //lower
				multiQubit.stateVec); //output
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


