/** @file
An implementation of the API in qubits.h for an MPI environment.
*/
# include <mpi.h>
# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "qubits.h"
# include "qubits_internal.h"

# define DEBUG 0
static int isChunkToSkipInFindPZero(int chunkId, int chunkSize, int measureQubit);
static int chunkIsUpper(int chunkId, int chunkSize, int rotQubit);
static void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta);
static int getChunkPairId(int chunkIsUpper, int chunkId, int chunkSize, int rotQubit);
static int halfMatrixBlockFitsInChunk(int chunkSize, int rotQubit);

void initQUESTEnv(QUESTEnv *env){
        // init MPI environment
        int rank, numRanks, initialized;
	MPI_Initialized(&initialized);
	if (!initialized){
		MPI_Init(NULL, NULL);
		MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		if (DEBUG) {
			char hostName[256];
			int hostNameLen;
			MPI_Get_processor_name(hostName, &hostNameLen);
			printf("rank %d on host %s\n", rank, hostName);
		}
		env->rank=rank;
		env->numRanks=numRanks;
	} else printf("ERROR: Trying to initialize QUESTEnv multiple times. Ignoring\n");
}

void syncQUESTEnv(QUESTEnv env){
	MPI_Barrier(MPI_COMM_WORLD);
}

void closeQUESTEnv(QUESTEnv env){
	int finalized;
	MPI_Finalized(&finalized);
	if (!finalized) MPI_Finalize();
	else printf("ERROR: Trying to close QUESTEnv multiple times. Ignoring\n");
}

void reportQUESTEnv(QUESTEnv env){
	if (env.rank==0){
		printf("EXECUTION ENVIRONMENT:\n"); 
		printf("Running distributed (MPI) version\n");
		printf("Number of ranks is %d\n", env.numRanks);
# ifdef _OPENMP
		printf("OpenMP enabled\n");
		printf("Number of threads available is %d\n", omp_get_max_threads());
# else
		printf("OpenMP disabled\n");
# endif 
	}
}



double calcTotalProbability(MultiQubit multiQubit){
        double pTotal=0; 
        double allRankTotals=0;
	long long int index;
	long long int numAmpsPerRank = multiQubit.numAmps;
        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];      
                pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];      
        } 
	if (DEBUG) printf("before calc prob. %d\n", multiQubit.numChunks);
	if (multiQubit.numChunks>1) MPI_Allreduce(&pTotal, &allRankTotals, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	else allRankTotals=pTotal;

	return allRankTotals;
}

/** Returns whether a given chunk in position chunkId is in the upper or lower half of
a block.

@param[in] chunkId id of chunk in state vector
@param[in] chunkSize number of amps in chunk
@param[in] rotQubit qubit being rotated 
@return 1: chunk is in upper half of block, 0: chunk is in lower half of block 
*/

static int chunkIsUpper(int chunkId, int chunkSize, int rotQubit)
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
static void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta)
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

static int getChunkPairId(int chunkIsUpper, int chunkId, int chunkSize, int rotQubit)
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

static int halfMatrixBlockFitsInChunk(int chunkSize, int rotQubit)
{
        long long int sizeHalfBlock = 1LL << (rotQubit);
        if (chunkSize > sizeHalfBlock) return 1;
        else return 0;
}

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

/** Find chunks to skip when calculating probability of qubit being zero.
When calculating probability of a bit q being zero,
sum up 2^q values, then skip 2^q values, etc. This function finds if an entire chunk
is in the range of values to be skipped

@param[in] chunkId id of chunk in state vector
@param[in] chunkSize number of amps in chunk
@param[in] measureQubi qubit being measured
@return int -- 1: skip, 0: don't skip
*/

static int isChunkToSkipInFindPZero(int chunkId, int chunkSize, int measureQubit)
{
        long long int sizeHalfBlock = 1LL << (measureQubit);
        int numChunksToSkip = sizeHalfBlock/chunkSize;
        // calculate probability by summing over numChunksToSkip, then skipping numChunksToSkip, etc
        int bitToCheck = chunkId & numChunksToSkip;
        return bitToCheck;
}

double findProbabilityOfZero(MultiQubit multiQubit, const int measureQubit)
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
	MPI_Allreduce(&stateProb, &totalStateProb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return totalStateProb;
}


double measureInZero(MultiQubit multiQubit, const int measureQubit)
{
	double totalStateProb=findProbabilityOfZero(multiQubit, measureQubit);
	int skipValuesWithinRank = halfMatrixBlockFitsInChunk(multiQubit.numAmps, measureQubit);
	if (skipValuesWithinRank) {
		measureInZeroLocal(multiQubit, measureQubit, totalStateProb);
	} else {
		if (!isChunkToSkipInFindPZero(multiQubit.chunkId, multiQubit.numAmps, measureQubit)){
			measureInZeroDistributedRenorm(multiQubit, measureQubit, totalStateProb);
		} else {
			measureInZeroDistributedSetZero(multiQubit, measureQubit);
		}
	}
	return totalStateProb;
}

double filterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
	double stateProb=0;
	stateProb = probOfFilterOut111(multiQubit, idQubit1, idQubit2, idQubit3);
	filterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3, stateProb);
	return stateProb;
}

double probOfFilterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
	double stateProb=0, totalStateProb=0;
	stateProb = probOfFilterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3);
	MPI_Allreduce(&stateProb, &totalStateProb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return totalStateProb;
}

