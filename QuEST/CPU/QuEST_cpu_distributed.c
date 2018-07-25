// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the pure backend in ../QuEST_ops_pure.h for an MPI environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 */

# include "../QuEST.h"
# include "../QuEST_internal.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"

# include "QuEST_cpu_internal.h"

#define _BSD_SOURCE
# include <unistd.h>
# include <mpi.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <sys/types.h>

# ifdef _OPENMP
# include <omp.h>
# endif


// @TODO
void densmatr_initPureState(QubitRegister targetQureg, QubitRegister copyQureg) {
	densmatr_initPureStateDistributed(targetQureg, copyQureg);
}

REAL densmatr_calcTotalProbability(QubitRegister qureg) {
	
	// computes the trace by summing every element ("diag") with global index (2^n + 1)i for i in [0, 2^n-1]
	
	// computes first local index containing a diagonal element
	long long int diagSpacing = 1LL + (1LL << qureg.numQubitsRepresented);
	long long int numPrevDiags = (qureg.chunkId * qureg.numAmpsPerChunk) / diagSpacing;
	long long int globalIndNextDiag = diagSpacing * numPrevDiags;
	long long int localIndNextDiag = globalIndNextDiag % qureg.numAmpsPerChunk;
	long long int index;
	
	REAL rankTotal = 0;
	REAL y, t, c;
	c = 0;
	
	// iterates every local diagonal
	for (index=localIndNextDiag; index < qureg.numAmpsPerChunk; index += diagSpacing) {
		
		// Kahan summation - brackets are important
		y = qureg.stateVec.real[index] - c;
		t = rankTotal + y;
		c = ( t - rankTotal ) - y;
		rankTotal = t;
	}
	
	// combine each node's sum of diagonals
	REAL globalTotal;
	if (qureg.numChunks > 1)
		MPI_Allreduce(&rankTotal, &globalTotal, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
	else
		globalTotal = rankTotal;
	
	return globalTotal;
}

REAL statevec_calcTotalProbability(QubitRegister qureg){
    // Implemented using Kahan summation for greater accuracy at a slight floating
    //   point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    REAL pTotal=0; 
    REAL y, t, c;
    REAL allRankTotals=0;
    long long int index;
    long long int numAmpsPerRank = qureg.numAmpsPerChunk;
    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){ 
        // Perform pTotal+=qureg.stateVec.real[index]*qureg.stateVec.real[index]; by Kahan
        y = qureg.stateVec.real[index]*qureg.stateVec.real[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;
        // Perform pTotal+=qureg.stateVec.imag[index]*qureg.stateVec.imag[index]; by Kahan
        y = qureg.stateVec.imag[index]*qureg.stateVec.imag[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;
    } 
    if (qureg.numChunks>1)
		MPI_Allreduce(&pTotal, &allRankTotals, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    else 
		allRankTotals=pTotal;

    return allRankTotals;
}


static int isChunkToSkipInFindPZero(int chunkId, long long int chunkSize, int measureQubit);
static int chunkIsUpper(int chunkId, long long int chunkSize, int targetQubit);
static void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta);
static int getChunkPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit);
static int halfMatrixBlockFitsInChunk(long long int chunkSize, int targetQubit);
static int getChunkIdFromIndex(QubitRegister qureg, long long int index);

void initQuESTEnv(QuESTEnv *env){
    // init MPI environment
    int rank, numRanks, initialized;
    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        env->rank=rank;
        env->numRanks=numRanks;

    } else printf("ERROR: Trying to initialize QuESTEnv multiple times. Ignoring\n");
	
	seedQuESTDefault();
}

void syncQuESTEnv(QuESTEnv env){
    MPI_Barrier(MPI_COMM_WORLD);
}

int syncQuESTSuccess(int successCode){
    int totalSuccess;
    MPI_Allreduce(&successCode, &totalSuccess, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return totalSuccess;
}

void closeQuESTEnv(QuESTEnv env){
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Finalize();
    else printf("ERROR: Trying to close QuESTEnv multiple times. Ignoring\n");
}

void reportQuESTEnv(QuESTEnv env){
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
        printf("Precision: size of REAL is %ld bytes\n", sizeof(REAL));
    }
}

void reportNodeList(QuESTEnv env){
    char hostName[256];
    gethostname(hostName, 255);
    printf("hostname on rank %d: %s\n", env.rank, hostName);
}

int getChunkIdFromIndex(QubitRegister qureg, long long int index){
    return index/qureg.numAmpsPerChunk; // this is numAmpsPerChunk
}

REAL statevec_getRealAmpEl(QubitRegister qureg, long long int index){
    int chunkId = getChunkIdFromIndex(qureg, index);
    REAL el; 
    if (qureg.chunkId==chunkId){
        el = qureg.stateVec.real[index-chunkId*qureg.numAmpsPerChunk];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
} 

REAL statevec_getImagAmpEl(QubitRegister qureg, long long int index){
    int chunkId = getChunkIdFromIndex(qureg, index);
    REAL el; 
    if (qureg.chunkId==chunkId){
        el = qureg.stateVec.imag[index-chunkId*qureg.numAmpsPerChunk];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
}

/** Returns whether a given chunk in position chunkId is in the upper or lower half of
  a block.
 * 
 * @param[in] chunkId id of chunk in state vector
 * @param[in] chunkSize number of amps in chunk
 * @param[in] targetQubit qubit being rotated 
 * @return 1: chunk is in upper half of block, 0: chunk is in lower half of block 
 */
//! fix -- is this the same as isChunkToSkip?
static int chunkIsUpper(int chunkId, long long int chunkSize, int targetQubit)
{       
    long long int sizeHalfBlock = 1LL << (targetQubit);
    long long int sizeBlock = sizeHalfBlock*2;
    long long int posInBlock = (chunkId*chunkSize) % sizeBlock;
    return posInBlock<sizeHalfBlock;
}

/** Get rotation values for a given chunk
 * @param[in] chunkIsUpper 1: chunk is in upper half of block, 0: chunk is in lower half
 * 
 * @param[out] rot1, rot2 rotation values to use, allocated for upper/lower such that
 * @verbatim
 stateUpper = rot1 * stateUpper + conj(rot2)  * stateLower
 @endverbatim
 * or
 * @verbatim
 stateLower = rot1 * stateUpper + conj(rot2)  * stateLower
 @endverbatim
 *
 * @param[in] alpha, beta initial rotation values 
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

/** Get rotation values for a given chunk given a unitary matrix
 * @param[in] chunkIsUpper 1: chunk is in upper half of block, 0: chunk is in lower half
 * 
 * @param[out] rot1, rot2 rotation values to use, allocated for upper/lower such that
 * @verbatim
 stateUpper = rot1 * stateUpper + conj(rot2)  * stateLower
 @endverbatim
 * or
 * @verbatim
 stateLower = rot1 * stateUpper + conj(rot2)  * stateLower
 @endverbatim
 * @param[in] u unitary matrix operation
 */
static void getRotAngleFromUnitaryMatrix(int chunkIsUpper, Complex *rot1, Complex *rot2, ComplexMatrix2 u)
{
    if (chunkIsUpper){
        *rot1=u.r0c0;
        *rot2=u.r0c1;
    } else {
        *rot1=u.r1c0;
        *rot2=u.r1c1;
    }
}

/** get position of corresponding chunk, holding values required to
 * update values in my chunk (with chunkId) when rotating targetQubit.
 * 
 * @param[in] chunkIsUpper 1: chunk is in upper half of block, 0: chunk is in lower half
 * @param[in] chunkId id of chunk in state vector
 * @param[in] chunkSize number of amps in chunk
 * @param[in] targetQubit qubit being rotated 
 * @return chunkId of chunk required to rotate targetQubit 
 */
static int getChunkPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit)
{
    long long int sizeHalfBlock = 1LL << (targetQubit);
    int chunksPerHalfBlock = sizeHalfBlock/chunkSize;
    if (chunkIsUpper){
        return chunkId + chunksPerHalfBlock;
    } else {
        return chunkId - chunksPerHalfBlock;
    }
}

/** return whether the current qubit rotation will use
 * blocks that fit within a single chunk.
 * 
 * @param[in] chunkSize number of amps in chunk
 * @param[in] targetQubit qubit being rotated 
 * @return 1: one chunk fits in one block 0: chunk is larger than block
 */
static int halfMatrixBlockFitsInChunk(long long int chunkSize, int targetQubit)
{
    long long int sizeHalfBlock = 1LL << (targetQubit);
    if (chunkSize > sizeHalfBlock) return 1;
    else return 0;
}

void exchangeStateVectors(QubitRegister qureg, int pairRank){
    // MPI send/receive vars
    int TAG=100;
    MPI_Status status;

    // Multiple messages are required as MPI uses int rather than long long int for count
    // For openmpi, messages are further restricted to 2GB in size -- do this for all cases
    // to be safe
    long long int maxMessageCount = 1LL<<29;
    if (sizeof(REAL)==8) maxMessageCount = (1LL<<28);
    else if (sizeof(REAL)==16) maxMessageCount = (1LL<<27);

    if (qureg.numAmpsPerChunk<maxMessageCount) maxMessageCount = qureg.numAmpsPerChunk;
    int numMessages = qureg.numAmpsPerChunk/maxMessageCount;
    int i;
    long long int offset;
    // send my state vector to pairRank's qureg.pairStateVec
    // receive pairRank's state vector into qureg.pairStateVec
    for (i=0; i<numMessages; i++){
        offset = i*maxMessageCount;
        MPI_Sendrecv(&qureg.stateVec.real[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &qureg.pairStateVec.real[offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
        //printf("rank: %d err: %d\n", qureg.rank, err);
        MPI_Sendrecv(&qureg.stateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &qureg.pairStateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
    }
}

void statevec_compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_compactUnitaryLocal(qureg, targetQubit, alpha, beta);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        getRotAngle(rankIsUpper, &rot1, &rot2, alpha, beta);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);

        // this rank's values are either in the upper of lower half of the block. 
        // send values to compactUnitaryDistributed in the correct order
        if (rankIsUpper){
            statevec_compactUnitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_compactUnitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void statevec_unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_unitaryLocal(qureg, targetQubit, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);

        // this rank's values are either in the upper of lower half of the block. 
        // send values to compactUnitaryDistributed in the correct order
        if (rankIsUpper){
            statevec_unitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_unitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }


}

void statevec_controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledCompactUnitaryLocal(qureg, controlQubit, targetQubit, alpha, beta);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        getRotAngle(rankIsUpper, &rot1, &rot2, alpha, beta);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to controlledCompactUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            statevec_controlledCompactUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_controlledCompactUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void statevec_controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, 
        ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledUnitaryLocal(qureg, controlQubit, targetQubit, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to controlledUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            statevec_controlledUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_controlledUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void statevec_multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(numControlQubits > 0 && numControlQubits <= qureg.numQubitsInStateVec, 4, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<qureg.numQubitsInStateVec)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<targetQubit)) != (1LL<<targetQubit), 3, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_multiControlledUnitaryLocal(qureg, targetQubit, mask, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to multiControlledUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            statevec_multiControlledUnitaryDistributed(qureg,targetQubit,mask,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_multiControlledUnitaryDistributed(qureg,targetQubit,mask,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}
void statevec_sigmaX(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_sigmaXLocal(qureg, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block. sigmaX just replaces
        // this rank's values with pair values
        statevec_sigmaXDistributed(qureg, targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec); // out
    }
}

void statevec_controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper; 	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledNotLocal(qureg, controlQubit, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        if (rankIsUpper){
            statevec_controlledNotDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec); //out
        } else {
            statevec_controlledNotDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec); //out
        }
    }
}

void statevec_sigmaY(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
	
	int conjFac = 1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper;	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        statevec_sigmaYLocal(qureg, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        statevec_sigmaYDistributed(qureg,targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper, conjFac);
    }
}

void statevec_sigmaYConj(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
	
	int conjFac = -1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper;	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        statevec_sigmaYLocal(qureg, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        statevec_sigmaYDistributed(qureg,targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper, conjFac);
    }
}

void statevec_controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
	
	int conjFac = 1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper; 	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledSigmaYLocal(qureg, controlQubit, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        if (rankIsUpper){
            statevec_controlledSigmaYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        } else {
            statevec_controlledSigmaYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        }
    }
}

void statevec_controlledSigmaYConj(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
	
	int conjFac = -1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper; 	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledSigmaYLocal(qureg, controlQubit, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        if (rankIsUpper){
            statevec_controlledSigmaYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        } else {
            statevec_controlledSigmaYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        }
    }
}

void statevec_hadamard(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubitsInStateVec, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_hadamardLocal(qureg, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block. send values to hadamardDistributed
        // in the correct order
        if (rankIsUpper){
            statevec_hadamardDistributed(qureg,targetQubit,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec, rankIsUpper); //output
        } else {
            statevec_hadamardDistributed(qureg,targetQubit,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec, rankIsUpper); //output
        }
    }
}


/** Find chunks to skip when calculating probability of qubit being zero.
 * When calculating probability of a bit q being zero,
 * sum up 2^q values, then skip 2^q values, etc. This function finds if an entire chunk
 * is in the range of values to be skipped
 * 
 * @param[in] chunkId id of chunk in state vector
 * @param[in] chunkSize number of amps in chunk
 * @param[in] measureQubi qubit being measured
 * @return int -- 1: skip, 0: don't skip
 */
static int isChunkToSkipInFindPZero(int chunkId, long long int chunkSize, int measureQubit)
{
    long long int sizeHalfBlock = 1LL << (measureQubit);
    int numChunksToSkip = sizeHalfBlock/chunkSize;
    // calculate probability by summing over numChunksToSkip, then skipping numChunksToSkip, etc
    int bitToCheck = chunkId & numChunksToSkip;
    return bitToCheck;
}

REAL statevec_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubitsInStateVec, 2, __func__);

    REAL stateProb=0, totalStateProb=0;
    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        stateProb = statevec_findProbabilityOfZeroLocal(qureg, measureQubit);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            stateProb = statevec_findProbabilityOfZeroDistributed(qureg, measureQubit);
        } else stateProb = 0;
    }
    MPI_Allreduce(&stateProb, &totalStateProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    if (outcome==1) totalStateProb = 1.0 - totalStateProb;
    return totalStateProb;
}

REAL densmatr_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome) {
	
	REAL zeroProb = densmatr_findProbabilityOfZeroLocal(qureg, measureQubit);
	
	REAL outcomeProb;
	MPI_Allreduce(&zeroProb, &outcomeProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
	if (outcome == 1)
		outcomeProb = 1.0 - outcomeProb;
	
	return outcomeProb;
}

REAL statevec_collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubitsInStateVec, 2, __func__);
    QuESTAssert((outcome==0 || outcome==1), 10, __func__);

    REAL totalStateProb = statevec_findProbabilityOfOutcome(qureg, measureQubit, outcome);
    QuESTAssert(absReal(totalStateProb)>REAL_EPS, 8, __func__);

    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        statevec_collapseToOutcomeLocal(qureg, measureQubit, totalStateProb, outcome);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) statevec_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    totalStateProb);
            else statevec_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        } else {
            // chunk has amps for q=1
            if (outcome==1) statevec_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    totalStateProb);
            else statevec_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        }
    }
    return totalStateProb;
}

int statevec_measureWithStats(QubitRegister qureg, int measureQubit, REAL *stateProb){
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubitsInStateVec, 2, __func__);

    int outcome;
    // find probability of qubit being in state 1
    REAL stateProbInternal = statevec_findProbabilityOfOutcome(qureg, measureQubit, 1);

    // we can't collapse to a state that has a probability too close to zero
    if (stateProbInternal<REAL_EPS) outcome=0;
    else if (1-stateProbInternal<REAL_EPS) outcome=1;
    else {
        // ok. both P(0) and P(1) are large enough to resolve
        // generate random float on [0,1]
        float randNum = genrand_real1();
        if (randNum<=stateProbInternal) outcome = 1;
        else outcome = 0;
    } 
    if (outcome==0) stateProbInternal = 1-stateProbInternal;

    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        statevec_collapseToOutcomeLocal(qureg, measureQubit, stateProbInternal, outcome);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) statevec_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    stateProbInternal);
            else statevec_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        } else {
            // chunk has amps for q=1
            if (outcome==1) statevec_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    stateProbInternal);
            else statevec_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        }
    }

    *stateProb = stateProbInternal;
    return outcome;
}

int statevec_measure(QubitRegister qureg, int measureQubit){
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubitsInStateVec, 2, __func__);
    REAL stateProb; 
    return statevec_measureWithStats(qureg, measureQubit, &stateProb); 
}
