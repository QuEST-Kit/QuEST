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
void mixed_initPureState(QubitRegister targetQureg, QubitRegister copyQureg) {
	mixed_initPureStateDistributed(targetQureg, copyQureg);
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

REAL pure_getRealAmpEl(QubitRegister qureg, long long int index){
    int chunkId = getChunkIdFromIndex(qureg, index);
    REAL el; 
    if (qureg.chunkId==chunkId){
        el = qureg.stateVec.real[index-chunkId*qureg.numAmpsPerChunk];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
} 

REAL pure_getImagAmpEl(QubitRegister qureg, long long int index){
    int chunkId = getChunkIdFromIndex(qureg, index);
    REAL el; 
    if (qureg.chunkId==chunkId){
        el = qureg.stateVec.imag[index-chunkId*qureg.numAmpsPerChunk];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
}

REAL pure_calcTotalProbability(QubitRegister qureg){
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
    if (qureg.numChunks>1) MPI_Allreduce(&pTotal, &allRankTotals, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
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

void pure_compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_compactUnitaryLocal(qureg, targetQubit, alpha, beta);
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
            pure_compactUnitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            pure_compactUnitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void pure_unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_unitaryLocal(qureg, targetQubit, u);
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
            pure_unitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            pure_unitaryDistributed(qureg,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }


}

void pure_controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubits, 2, __func__);
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
        pure_controlledCompactUnitaryLocal(qureg, controlQubit, targetQubit, alpha, beta);
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
            pure_controlledCompactUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            pure_controlledCompactUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void pure_controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, 
        ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubits, 2, __func__);
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
        pure_controlledUnitaryLocal(qureg, controlQubit, targetQubit, u);
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
            pure_controlledUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            pure_controlledUnitaryDistributed(qureg,controlQubit,targetQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void pure_multiControlledUnitary(QubitRegister qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);
    QuESTAssert(numControlQubits > 0 && numControlQubits <= qureg.numQubits, 4, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<qureg.numQubits)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<targetQubit)) != (1LL<<targetQubit), 3, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_multiControlledUnitaryLocal(qureg, targetQubit, mask, u);
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
            pure_multiControlledUnitaryDistributed(qureg,targetQubit,mask,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            pure_multiControlledUnitaryDistributed(qureg,targetQubit,mask,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}
void pure_sigmaX(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_sigmaXLocal(qureg, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block. sigmaX just replaces
        // this rank's values with pair values
        pure_sigmaXDistributed(qureg, targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec); // out
    }
}

void pure_controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < qureg.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_controlledNotLocal(qureg, controlQubit, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block. send values to controlledNot
        // in the correct order
        if (rankIsUpper){
            pure_controlledNotDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec); //out
        } else {
            pure_controlledNotDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec); //out
        }
    }
}

void pure_sigmaY(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_sigmaYLocal(qureg, targetQubit);
    } else {
        //! fix -- put duplicate code (sigmaX, sigmaY) in seperate function
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block. sigmaX just replaces
        // this rank's values with pair values
        pure_sigmaYDistributed(qureg,targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper);
    }
}

// @ probably wrong
void pure_phaseShift(QubitRegister qureg, const int targetQubit, REAL angle)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;

    if (useLocalDataOnly){
        pure_phaseShiftLocal(qureg, targetQubit, angle);
    } else {
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        if (!rankIsUpper) pure_phaseShiftDistributed(qureg, targetQubit, angle);
    }
}

void pure_hadamard(QubitRegister qureg, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < qureg.numQubits, 1, __func__);

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        pure_hadamardLocal(qureg, targetQubit);
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
            pure_hadamardDistributed(qureg,targetQubit,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec, rankIsUpper); //output
        } else {
            pure_hadamardDistributed(qureg,targetQubit,
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

REAL pure_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubits, 2, __func__);

    REAL stateProb=0, totalStateProb=0;
    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        stateProb = pure_findProbabilityOfZeroLocal(qureg, measureQubit);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            stateProb = pure_findProbabilityOfZeroDistributed(qureg, measureQubit);
        } else stateProb = 0;
    }
    MPI_Allreduce(&stateProb, &totalStateProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    if (outcome==1) totalStateProb = 1.0 - totalStateProb;
    return totalStateProb;
}


REAL pure_collapseToOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubits, 2, __func__);
    QuESTAssert((outcome==0 || outcome==1), 10, __func__);

    REAL totalStateProb = pure_findProbabilityOfOutcome(qureg, measureQubit, outcome);
    QuESTAssert(absReal(totalStateProb)>REAL_EPS, 8, __func__);

    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        pure_collapseToOutcomeLocal(qureg, measureQubit, totalStateProb, outcome);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) pure_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    totalStateProb);
            else pure_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        } else {
            // chunk has amps for q=1
            if (outcome==1) pure_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    totalStateProb);
            else pure_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        }
    }
    return totalStateProb;
}

int pure_measureWithStats(QubitRegister qureg, int measureQubit, REAL *stateProb){
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubits, 2, __func__);

    int outcome;
    // find probability of qubit being in state 1
    REAL stateProbInternal = pure_findProbabilityOfOutcome(qureg, measureQubit, 1);

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
        pure_collapseToOutcomeLocal(qureg, measureQubit, stateProbInternal, outcome);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) pure_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    stateProbInternal);
            else pure_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        } else {
            // chunk has amps for q=1
            if (outcome==1) pure_collapseToOutcomeDistributedRenorm(qureg, measureQubit, 
                    stateProbInternal);
            else pure_collapseToOutcomeDistributedSetZero(qureg, measureQubit);
        }
    }

    *stateProb = stateProbInternal;
    return outcome;
}

int pure_measure(QubitRegister qureg, int measureQubit){
    QuESTAssert(measureQubit >= 0 && measureQubit < qureg.numQubits, 2, __func__);
    REAL stateProb; 
    return pure_measureWithStats(qureg, measureQubit, &stateProb); 
}
