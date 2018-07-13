// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt 
// for details 

/** @file
  An implementation of the API in qubits.h for an MPI environment.
  */

# include "../QuEST.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"
# include "QuEST_internal.h"

#define _BSD_SOURCE // defining prototype for gethostname 
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

static int isChunkToSkipInFindPZero(int chunkId, long long int chunkSize, int measureQubit);
static int chunkIsUpper(int chunkId, long long int chunkSize, int targetQubit);
static void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta);
static int getChunkPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit);
static int halfMatrixBlockFitsInChunk(long long int chunkSize, int targetQubit);
static int getChunkIdFromIndex(MultiQubit multiQubit, long long int index);

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
	
	QuESTSeedRandomDefault();
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

int getChunkIdFromIndex(MultiQubit multiQubit, long long int index){
    return index/multiQubit.numAmpsDividedByNumChunks; // this is numAmpsPerChunk
}

REAL getRealAmpEl(MultiQubit multiQubit, long long int index){
    int chunkId = getChunkIdFromIndex(multiQubit, index);
    REAL el; 
    if (multiQubit.chunkId==chunkId){
        el = multiQubit.stateVec.real[index-chunkId*multiQubit.numAmpsDividedByNumChunks];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
} 

REAL getImagAmpEl(MultiQubit multiQubit, long long int index){
    int chunkId = getChunkIdFromIndex(multiQubit, index);
    REAL el; 
    if (multiQubit.chunkId==chunkId){
        el = multiQubit.stateVec.imag[index-chunkId*multiQubit.numAmpsDividedByNumChunks];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
}

REAL calcTotalProbability(MultiQubit multiQubit){
    // Implemented using Kahan summation for greater accuracy at a slight floating
    //   point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    REAL pTotal=0; 
    REAL y, t, c;
    REAL allRankTotals=0;
    long long int index;
    long long int numAmpsPerRank = multiQubit.numAmpsDividedByNumChunks;
    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){ 
        // Perform pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index]; by Kahan
        y = multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;
        // Perform pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index]; by Kahan
        y = multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;
    } 
    if (multiQubit.numChunks>1) MPI_Allreduce(&pTotal, &allRankTotals, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    else allRankTotals=pTotal;

    return allRankTotals;
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

void exchangeStateVectors(MultiQubit multiQubit, int pairRank){
    // MPI send/receive vars
    int TAG=100;
    MPI_Status status;

    // Multiple messages are required as MPI uses int rather than long long int for count
    // For openmpi, messages are further restricted to 2GB in size -- do this for all cases
    // to be safe
    long long int maxMessageCount = 1LL<<29;
    if (sizeof(REAL)==8) maxMessageCount = (1LL<<28);
    else if (sizeof(REAL)==16) maxMessageCount = (1LL<<27);

    if (multiQubit.numAmpsDividedByNumChunks<maxMessageCount) maxMessageCount = multiQubit.numAmpsDividedByNumChunks;
    int numMessages = multiQubit.numAmpsDividedByNumChunks/maxMessageCount;
    int i;
    long long int offset;
    // send my state vector to pairRank's multiQubit.pairStateVec
    // receive pairRank's state vector into multiQubit.pairStateVec
    for (i=0; i<numMessages; i++){
        offset = i*maxMessageCount;
        MPI_Sendrecv(&multiQubit.stateVec.real[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &multiQubit.pairStateVec.real[offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
        //printf("rank: %d err: %d\n", multiQubit.rank, err);
        MPI_Sendrecv(&multiQubit.stateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, TAG,
                &multiQubit.pairStateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
    }
}

void compactUnitary(MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        compactUnitaryLocal(multiQubit, targetQubit, alpha, beta);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        getRotAngle(rankIsUpper, &rot1, &rot2, alpha, beta);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);

        // this rank's values are either in the upper of lower half of the block. 
        // send values to compactUnitaryDistributed in the correct order
        if (rankIsUpper){
            compactUnitaryDistributed(multiQubit,targetQubit,rot1,rot2,
                    multiQubit.stateVec, //upper
                    multiQubit.pairStateVec, //lower
                    multiQubit.stateVec); //output
        } else {
            compactUnitaryDistributed(multiQubit,targetQubit,rot1,rot2,
                    multiQubit.pairStateVec, //upper
                    multiQubit.stateVec, //lower
                    multiQubit.stateVec); //output
        }
    }
}

void unitary(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        unitaryLocal(multiQubit, targetQubit, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);

        // this rank's values are either in the upper of lower half of the block. 
        // send values to compactUnitaryDistributed in the correct order
        if (rankIsUpper){
            unitaryDistributed(multiQubit,targetQubit,rot1,rot2,
                    multiQubit.stateVec, //upper
                    multiQubit.pairStateVec, //lower
                    multiQubit.stateVec); //output
        } else {
            unitaryDistributed(multiQubit,targetQubit,rot1,rot2,
                    multiQubit.pairStateVec, //upper
                    multiQubit.stateVec, //lower
                    multiQubit.stateVec); //output
        }
    }


}

void controlledCompactUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        controlledCompactUnitaryLocal(multiQubit, controlQubit, targetQubit, alpha, beta);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        getRotAngle(rankIsUpper, &rot1, &rot2, alpha, beta);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to controlledCompactUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            controlledCompactUnitaryDistributed(multiQubit,controlQubit,targetQubit,rot1,rot2,
                    multiQubit.stateVec, //upper
                    multiQubit.pairStateVec, //lower
                    multiQubit.stateVec); //output
        } else {
            controlledCompactUnitaryDistributed(multiQubit,controlQubit,targetQubit,rot1,rot2,
                    multiQubit.pairStateVec, //upper
                    multiQubit.stateVec, //lower
                    multiQubit.stateVec); //output
        }
    }
}

void controlledUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, 
        ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        controlledUnitaryLocal(multiQubit, controlQubit, targetQubit, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to controlledUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            controlledUnitaryDistributed(multiQubit,controlQubit,targetQubit,rot1,rot2,
                    multiQubit.stateVec, //upper
                    multiQubit.pairStateVec, //lower
                    multiQubit.stateVec); //output
        } else {
            controlledUnitaryDistributed(multiQubit,controlQubit,targetQubit,rot1,rot2,
                    multiQubit.pairStateVec, //upper
                    multiQubit.stateVec, //lower
                    multiQubit.stateVec); //output
        }
    }
}

void multiControlledUnitary(MultiQubit multiQubit, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(numControlQubits > 0 && numControlQubits <= multiQubit.numQubits, 4, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<multiQubit.numQubits)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<targetQubit)) != (1LL<<targetQubit), 3, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        multiControlledUnitaryLocal(multiQubit, targetQubit, mask, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to multiControlledUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            multiControlledUnitaryDistributed(multiQubit,targetQubit,mask,rot1,rot2,
                    multiQubit.stateVec, //upper
                    multiQubit.pairStateVec, //lower
                    multiQubit.stateVec); //output
        } else {
            multiControlledUnitaryDistributed(multiQubit,targetQubit,mask,rot1,rot2,
                    multiQubit.pairStateVec, //upper
                    multiQubit.stateVec, //lower
                    multiQubit.stateVec); //output
        }
    }
}
void sigmaX(MultiQubit multiQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        sigmaXLocal(multiQubit, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);
        // this rank's values are either in the upper of lower half of the block. sigmaX just replaces
        // this rank's values with pair values
        sigmaXDistributed(multiQubit, targetQubit,
                multiQubit.pairStateVec, // in
                multiQubit.stateVec); // out
    }
}

void controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        controlledNotLocal(multiQubit, controlQubit, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);
        // this rank's values are either in the upper of lower half of the block. send values to controlledNot
        // in the correct order
        if (rankIsUpper){
            controlledNotDistributed(multiQubit,controlQubit,targetQubit,
                    multiQubit.pairStateVec, //in
                    multiQubit.stateVec); //out
        } else {
            controlledNotDistributed(multiQubit,controlQubit,targetQubit,
                    multiQubit.pairStateVec, //in
                    multiQubit.stateVec); //out
        }
    }
}

void sigmaY(MultiQubit multiQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        sigmaYLocal(multiQubit, targetQubit);
    } else {
        //! fix -- put duplicate code (sigmaX, sigmaY) in seperate function
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);
        // this rank's values are either in the upper of lower half of the block. sigmaX just replaces
        // this rank's values with pair values
        sigmaYDistributed(multiQubit,targetQubit,
                multiQubit.pairStateVec, // in
                multiQubit.stateVec, // out
                rankIsUpper);
    }
}

void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;

    if (useLocalDataOnly){
        phaseGateLocal(multiQubit, targetQubit, type);
    } else {
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        if (!rankIsUpper) phaseGateDistributed(multiQubit, targetQubit, type);
    }
}

void hadamard(MultiQubit multiQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        hadamardLocal(multiQubit, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, targetQubit);
        //printf("%d rank has pair rank: %d\n", multiQubit.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(multiQubit, pairRank);
        // this rank's values are either in the upper of lower half of the block. send values to hadamardDistributed
        // in the correct order
        if (rankIsUpper){
            hadamardDistributed(multiQubit,targetQubit,
                    multiQubit.stateVec, //upper
                    multiQubit.pairStateVec, //lower
                    multiQubit.stateVec, rankIsUpper); //output
        } else {
            hadamardDistributed(multiQubit,targetQubit,
                    multiQubit.pairStateVec, //upper
                    multiQubit.stateVec, //lower
                    multiQubit.stateVec, rankIsUpper); //output
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

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);

    REAL stateProb=0, totalStateProb=0;
    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, measureQubit);
    if (skipValuesWithinRank) {
        stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
    } else {
        if (!isChunkToSkipInFindPZero(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, measureQubit)){
            stateProb = findProbabilityOfZeroDistributed(multiQubit, measureQubit);
        } else stateProb = 0;
    }
    MPI_Allreduce(&stateProb, &totalStateProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    if (outcome==1) totalStateProb = 1.0 - totalStateProb;
    return totalStateProb;
}


REAL collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert((outcome==0 || outcome==1), 10, __func__);

    REAL totalStateProb=findProbabilityOfOutcome(multiQubit, measureQubit, outcome);
    QuESTAssert(fabs(totalStateProb)>REAL_EPS, 8, __func__);

    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, measureQubit);
    if (skipValuesWithinRank) {
        collapseToOutcomeLocal(multiQubit, measureQubit, totalStateProb, outcome);
    } else {
        if (!isChunkToSkipInFindPZero(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) collapseToOutcomeDistributedRenorm(multiQubit, measureQubit, 
                    totalStateProb);
            else collapseToOutcomeDistributedSetZero(multiQubit, measureQubit);
        } else {
            // chunk has amps for q=1
            if (outcome==1) collapseToOutcomeDistributedRenorm(multiQubit, measureQubit, 
                    totalStateProb);
            else collapseToOutcomeDistributedSetZero(multiQubit, measureQubit);
        }
    }
    return totalStateProb;
}


int measure(MultiQubit multiQubit, int measureQubit){
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    REAL stateProb; 
    return measureWithStats(multiQubit, measureQubit, &stateProb); 
}

int measureWithStats(MultiQubit multiQubit, int measureQubit, REAL *stateProb){
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);

    int outcome;
    // find probability of qubit being in state 1
    REAL stateProbInternal = findProbabilityOfOutcome(multiQubit, measureQubit, 1);

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

    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(multiQubit.numAmpsDividedByNumChunks, measureQubit);
    if (skipValuesWithinRank) {
        collapseToOutcomeLocal(multiQubit, measureQubit, stateProbInternal, outcome);
    } else {
        if (!isChunkToSkipInFindPZero(multiQubit.chunkId, multiQubit.numAmpsDividedByNumChunks, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) collapseToOutcomeDistributedRenorm(multiQubit, measureQubit, 
                    stateProbInternal);
            else collapseToOutcomeDistributedSetZero(multiQubit, measureQubit);
        } else {
            // chunk has amps for q=1
            if (outcome==1) collapseToOutcomeDistributedRenorm(multiQubit, measureQubit, 
                    stateProbInternal);
            else collapseToOutcomeDistributedSetZero(multiQubit, measureQubit);
        }
    }

    *stateProb = stateProbInternal;
    return outcome;
}

void exitWithError(int errorCode, const char* func){
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, errorCodes[errorCode]);
    printf("!!!\n");
    printf("exiting..\n");
    MPI_Abort(MPI_COMM_WORLD, errorCode);
}

void QuESTAssert(int isValid, int errorCode, const char* func){
    if (!isValid) exitWithError(errorCode, func);
}
