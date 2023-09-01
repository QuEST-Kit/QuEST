// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the backend in ../QuEST_ops.h for an MPI environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 *
 * @author Ania Brown
 * @author Tyson Jones
 * @author Balint Koczor
 */

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "mt19937ar.h"

# include "QuEST_cpu_internal.h"

# define _BSD_SOURCE
# include <unistd.h>
# include <mpi.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>    // for memcpy
# include <math.h>
# include <time.h>
# include <sys/types.h>

# ifdef _OPENMP
# include <omp.h>
# endif


Complex statevec_calcInnerProduct(Qureg bra, Qureg ket) {
    
    Complex localInnerProd = statevec_calcInnerProductLocal(bra, ket);
    if (bra.numChunks == 1)
        return localInnerProd;
    
    qreal localReal = localInnerProd.real;
    qreal localImag = localInnerProd.imag;
    qreal globalReal, globalImag;
    MPI_Allreduce(&localReal, &globalReal, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localImag, &globalImag, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    Complex globalInnerProd;
    globalInnerProd.real = globalReal;
    globalInnerProd.imag = globalImag;
    return globalInnerProd;
}

qreal densmatr_calcTotalProb(Qureg qureg) {
	
	// computes the trace by summing every element ("diag") with global index (2^n + 1)i for i in [0, 2^n-1]
    
    long long int firstGlobalInd = qureg.chunkId * qureg.numAmpsPerChunk;
    long long int diagSpacing = 1LL + (1LL << qureg.numQubitsRepresented);
    long long int distFromPrevDiag = (firstGlobalInd % diagSpacing);
    long long int distToNextDiag = (distFromPrevDiag > 0)? diagSpacing - distFromPrevDiag : 0;
	
    // Kahan parameters
	qreal rankTotal = 0;
	qreal y, t, c;
	c = 0;
    
    long long int index = distToNextDiag;
    
    while (index < qureg.numAmpsPerChunk) {
        
        // Kahan summation - brackets are important
        y = qureg.stateVec.real[index] - c;
        t = rankTotal + y;
        c = ( t - rankTotal ) - y;
        rankTotal = t;
        
        index += diagSpacing;
    }

	// combine each node's sum of diagonals
	qreal globalTotal;
	if (qureg.numChunks > 1)
		MPI_Allreduce(&rankTotal, &globalTotal, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
	else
		globalTotal = rankTotal;
	
	return globalTotal;
}

qreal statevec_calcTotalProb(Qureg qureg){
    // Implemented using Kahan summation for greater accuracy at a slight floating
    //   point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    qreal pTotal=0; 
    qreal y, t, c;
    qreal allRankTotals=0;
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
static int chunkIsUpperInOuterBlock(int chunkId, long long int chunkSize, int targetQubit, int numQubits);
static void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta);
static int getChunkPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit);
static int getChunkOuterBlockPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit, int numQubits);
static int halfMatrixBlockFitsInChunk(long long int chunkSize, int targetQubit);
static int getChunkIdFromIndex(Qureg qureg, long long int index);

QuESTEnv createQuESTEnv(void) {
    
    QuESTEnv env;
    
    // init MPI environment
    int rank, numRanks, initialized;
    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        env.rank=rank;
        env.numRanks=numRanks;

    } else {
        
        printf("ERROR: Trying to initialize QuESTEnv multiple times. Ignoring...\n");
        
        // ensure env is initialised anyway, so the compiler is happy
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        env.rank=rank;
        env.numRanks=numRanks;
	}
    
    validateNumRanks(env.numRanks, __func__);
    
    env.seeds = NULL;
    env.numSeeds = 0;
	seedQuESTDefault(&env);
    
    return env;
}

void syncQuESTEnv(QuESTEnv env){
    MPI_Barrier(MPI_COMM_WORLD);
}

int syncQuESTSuccess(int successCode){
    int totalSuccess;
    MPI_Allreduce(&successCode, &totalSuccess, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return totalSuccess;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);
    
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
        printf("Precision: size of qreal is %ld bytes\n", sizeof(qreal) );
    }
}

void getEnvironmentString(QuESTEnv env, char str[200]){
    int ompStatus=0;
    int numThreads=1;
# ifdef _OPENMP
    ompStatus=1;
    numThreads=omp_get_max_threads(); 
# endif
    sprintf(str, "CUDA=0 OpenMP=%d MPI=1 threads=%d ranks=%d", ompStatus, numThreads, env.numRanks);
}

int getChunkIdFromIndex(Qureg qureg, long long int index){
    return index/qureg.numAmpsPerChunk; // this is numAmpsPerChunk
}

qreal statevec_getRealAmp(Qureg qureg, long long int index){
    int chunkId = getChunkIdFromIndex(qureg, index);
    qreal el; 
    if (qureg.chunkId==chunkId){
        el = qureg.stateVec.real[index-chunkId*qureg.numAmpsPerChunk];
    }
    MPI_Bcast(&el, 1, MPI_QuEST_REAL, chunkId, MPI_COMM_WORLD);
    return el; 
} 

qreal statevec_getImagAmp(Qureg qureg, long long int index){
    int chunkId = getChunkIdFromIndex(qureg, index);
    qreal el; 
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

//! fix -- do with masking instead
static int chunkIsUpperInOuterBlock(int chunkId, long long int chunkSize, int targetQubit, int numQubits)
{       
    long long int sizeOuterHalfBlock = 1LL << (targetQubit+numQubits);
    long long int sizeOuterBlock = sizeOuterHalfBlock*2;
    long long int posInBlock = (chunkId*chunkSize) % sizeOuterBlock;
    return posInBlock<sizeOuterHalfBlock;
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
        *rot1=(Complex) {.real=u.real[0][0], .imag=u.imag[0][0]};
        *rot2=(Complex) {.real=u.real[0][1], .imag=u.imag[0][1]};
    } else {
        *rot1=(Complex) {.real=u.real[1][0], .imag=u.imag[1][0]};
        *rot2=(Complex) {.real=u.real[1][1], .imag=u.imag[1][1]};
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

static int getChunkOuterBlockPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit, int numQubits)
{
    long long int sizeOuterHalfBlock = 1LL << (targetQubit+numQubits);
    int chunksPerOuterHalfBlock = sizeOuterHalfBlock/chunkSize;
    if (chunkIsUpper){
        return chunkId + chunksPerOuterHalfBlock;
    } else {
        return chunkId - chunksPerOuterHalfBlock;
    }
}

static int getChunkOuterBlockPairIdForPart3(int chunkIsUpperSmallerQubit, int chunkIsUpperBiggerQubit, int chunkId, 
        long long int chunkSize, int smallerQubit, int biggerQubit, int numQubits)
{
    long long int sizeOuterHalfBlockBiggerQubit = 1LL << (biggerQubit+numQubits);
    long long int sizeOuterHalfBlockSmallerQubit = 1LL << (smallerQubit+numQubits);
    int chunksPerOuterHalfBlockSmallerQubit = sizeOuterHalfBlockSmallerQubit/chunkSize;
    int chunksPerOuterHalfBlockBiggerQubit = sizeOuterHalfBlockBiggerQubit/chunkSize;
    int rank;
    if (chunkIsUpperBiggerQubit){
        rank = chunkId + chunksPerOuterHalfBlockBiggerQubit;
    } else {
        rank = chunkId - chunksPerOuterHalfBlockBiggerQubit;
    }

    if (chunkIsUpperSmallerQubit){
        rank = rank + chunksPerOuterHalfBlockSmallerQubit;
    } else {
        rank = rank - chunksPerOuterHalfBlockSmallerQubit;
    }

    return rank;
}

/** return whether the current qubit rotation will use
 * blocks that fit within a single chunk.
 * 
 * @param[in] chunkSize number of amps in chunk
 * @param[in] targetQubit qubit being rotated 
 * @return 1: one chunk fits in one block 0: chunk is larger than block
 */
//! fix -- this should be renamed to matrixBlockFitsInChunk
static int halfMatrixBlockFitsInChunk(long long int chunkSize, int targetQubit)
{
    long long int sizeHalfBlock = 1LL << (targetQubit);
    if (chunkSize > sizeHalfBlock) return 1;
    else return 0;
}

static int densityMatrixBlockFitsInChunk(long long int chunkSize, int numQubits, int targetQubit) {
    long long int sizeOuterHalfBlock = 1LL << (targetQubit+numQubits);
    if (chunkSize > sizeOuterHalfBlock) return 1;
    else return 0;
}

/** This copies/clones vec (a statevector) into every node's matr pairState.
 * (where matr is a density matrix or equal number of qubits as vec) */
void copyVecIntoMatrixPairState(Qureg matr, Qureg vec) {
    
    // Remember that for every amplitude that `vec` stores on the node,
    // `matr` stores an entire column. Ergo there are always an integer
    // number (in fact, a power of 2) number of  `matr`s columns on each node.
    // Since the total size of `vec` (between all nodes) is one column
    // and each node stores (possibly) multiple columns (vec.numAmpsPerChunk as many), 
    // `vec` can be fit entirely inside a single node's matr.pairStateVec (with excess!)
    
    // copy this node's vec segment into this node's matr pairState (in the right spot)
    long long int numLocalAmps = vec.numAmpsPerChunk;
    long long int myOffset = vec.chunkId * numLocalAmps;
    memcpy(&matr.pairStateVec.real[myOffset], vec.stateVec.real, numLocalAmps * sizeof(qreal));
    memcpy(&matr.pairStateVec.imag[myOffset], vec.stateVec.imag, numLocalAmps * sizeof(qreal));
    
    // we now want to share this node's vec segment with other node, so that 
    // vec is cloned in every node's matr.pairStateVec 

    // work out how many messages needed to send vec chunks (2GB limit)
    long long int maxMsgSize = MPI_MAX_AMPS_IN_MSG;
    if (numLocalAmps < maxMsgSize) 
        maxMsgSize = numLocalAmps;
    // safely assume MPI_MAX... = 2^n, so division always exact:
    int numMsgs = numLocalAmps / maxMsgSize;
    
    // every node gets a turn at being the broadcaster
    for (int broadcaster=0; broadcaster < vec.numChunks; broadcaster++) {
        
        long long int otherOffset = broadcaster * numLocalAmps;
    
        // every node sends a slice of qureg's pairState to every other
        for (int i=0; i< numMsgs; i++) {
    
            // by sending that slice in further slices (due to bandwidth limit)
            MPI_Bcast(
                &matr.pairStateVec.real[otherOffset + i*maxMsgSize], 
                maxMsgSize,  MPI_QuEST_REAL, broadcaster, MPI_COMM_WORLD);
            MPI_Bcast(
                &matr.pairStateVec.imag[otherOffset + i*maxMsgSize], 
                maxMsgSize,  MPI_QuEST_REAL, broadcaster, MPI_COMM_WORLD);
        }
    }
}

qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState) {
    
    // set qureg's pairState is to be the full pureState (on every node)
    copyVecIntoMatrixPairState(qureg, pureState);
 
    // collect calcFidelityLocal by every machine
    qreal localSum = densmatr_calcFidelityLocal(qureg, pureState);
    
    // sum each localSum
    qreal globalSum;
    MPI_Allreduce(&localSum, &globalSum, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    return globalSum;
}

qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b) {
    
    qreal localSum = densmatr_calcHilbertSchmidtDistanceSquaredLocal(a, b);
    
    qreal globalSum;
    MPI_Allreduce(&localSum, &globalSum, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    qreal dist = sqrt(globalSum);
    return dist;
}

qreal densmatr_calcInnerProduct(Qureg a, Qureg b) {
    
    qreal localSum = densmatr_calcInnerProductLocal(a, b);
    
    qreal globalSum;
    MPI_Allreduce(&localSum, &globalSum, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    qreal dist = globalSum;
    return dist;
}

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg) {

    if (targetQureg.numChunks==1){
        // local version
        // save pointers to qureg's pair state
        qreal* quregPairRePtr = targetQureg.pairStateVec.real;
        qreal* quregPairImPtr = targetQureg.pairStateVec.imag;

        // populate qureg pair state with pure state (by repointing)
        targetQureg.pairStateVec.real = copyQureg.stateVec.real;
        targetQureg.pairStateVec.imag = copyQureg.stateVec.imag;

        // populate density matrix via it's pairState
        densmatr_initPureStateLocal(targetQureg, copyQureg);

        // restore pointers
        targetQureg.pairStateVec.real = quregPairRePtr;
        targetQureg.pairStateVec.imag = quregPairImPtr;
    } else {
        // set qureg's pairState is to be the full pure state (on every node)
        copyVecIntoMatrixPairState(targetQureg, copyQureg);
        
        // update every density matrix chunk using pairState
        densmatr_initPureStateLocal(targetQureg, copyQureg);
    }
}

void exchangeStateVectors(Qureg qureg, int pairRank){
    // MPI send/receive vars
    int sendTag;
    int recvTag;
    int rank;
    MPI_Request * requests;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Multiple messages are required as MPI uses int rather than long long int for count
    // For openmpi, messages are further restricted to 2GB in size -- do this for all cases
    // to be safe
    long long int maxMessageCount = MPI_MAX_AMPS_IN_MSG;
    if (qureg.numAmpsPerChunk < maxMessageCount) 
        maxMessageCount = qureg.numAmpsPerChunk;
    
    // safely assume MPI_MAX... = 2^n, so division always exact
    int numMessages = qureg.numAmpsPerChunk/maxMessageCount;
    requests = (MPI_Request*) malloc(4 * numMessages * sizeof(MPI_Request));
    int i;
    long long int offset;
    
    // send my state vector to pairRank's qureg.pairStateVec
    // receive pairRank's state vector into qureg.pairStateVec
    for (i=0; i<numMessages; i++){
        offset = i*maxMessageCount;
        sendTag = rank * numMessages + i;
        recvTag = pairRank * numMessages + i;

        MPI_Isend(&qureg.stateVec.real[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, sendTag, MPI_COMM_WORLD, requests + 4*i);
        MPI_Isend(&qureg.stateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, sendTag, MPI_COMM_WORLD, requests + 4*i + 1);
        MPI_Irecv(&qureg.pairStateVec.real[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, recvTag, MPI_COMM_WORLD, requests + 4*i + 2);
        MPI_Irecv(&qureg.pairStateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL, pairRank, recvTag, MPI_COMM_WORLD, requests + 4*i + 3);
    }

    MPI_Waitall(4 * numMessages, requests, MPI_STATUSES_IGNORE);

    free(requests);
}

void exchangePairStateVectorHalves(Qureg qureg, int pairRank){
    // MPI send/receive vars
    int TAG=100;
    MPI_Status status;
    long long int numAmpsToSend = qureg.numAmpsPerChunk >> 1;

    // Multiple messages are required as MPI uses int rather than long long int for count
    // For openmpi, messages are further restricted to 2GB in size -- do this for all cases
    // to be safe
    long long int maxMessageCount = MPI_MAX_AMPS_IN_MSG;
    if (numAmpsToSend < maxMessageCount) 
        maxMessageCount = numAmpsToSend;
    
    // safely assume MPI_MAX... = 2^n, so division always exact
    int numMessages = numAmpsToSend/maxMessageCount;
    int i;
    long long int offset;
    // send the bottom half of my state vector to the top half of pairRank's qureg.pairStateVec
    // receive pairRank's state vector into the top of qureg.pairStateVec
    for (i=0; i<numMessages; i++){
        offset = i*maxMessageCount;
        MPI_Sendrecv(&qureg.pairStateVec.real[offset+numAmpsToSend], maxMessageCount, 
                MPI_QuEST_REAL, pairRank, TAG,
                &qureg.pairStateVec.real[offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
        //printf("rank: %d err: %d\n", qureg.rank, err);
        MPI_Sendrecv(&qureg.pairStateVec.imag[offset+numAmpsToSend], maxMessageCount, 
                MPI_QuEST_REAL, pairRank, TAG,
                &qureg.pairStateVec.imag[offset], maxMessageCount, MPI_QuEST_REAL,
                pairRank, TAG, MPI_COMM_WORLD, &status);
    }
}

// @todo decide where this function should go. It is a preparation for MPI data transfer function
void compressPairVectorForSingleQubitDepolarise(Qureg qureg, int targetQubit){
    long long int sizeInnerBlock, sizeInnerHalfBlock;
    long long int sizeOuterColumn, sizeOuterHalfColumn;
    long long int thisInnerBlock, // current block
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInOuterColumn,
         thisIndexInInnerBlock;
         
    int outerBit;

    long long int thisTask;
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeInnerHalfBlock = 1LL << targetQubit;
    sizeInnerBlock     = 2LL * sizeInnerHalfBlock;
    sizeOuterHalfColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterColumn     = 2LL * sizeOuterHalfColumn;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlock,sizeInnerHalfBlock,sizeOuterColumn,sizeOuterHalfColumn, \
                qureg,numTasks,targetQubit) \
    private  (thisTask,thisInnerBlock,thisOuterColumn,thisIndex,thisIndexInOuterColumn, \
                thisIndexInInnerBlock,outerBit) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // thisTask iterates over half the elements in this process' chunk of the density matrix
        // treat this as iterating over all columns, then iterating over half the values
        // within one column.
        // If this function has been called, this process' chunk contains half an 
        // outer block or less
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // we want to process all columns in the density matrix,
            // updating the values for half of each column (one half of each inner block)
            thisOuterColumn = thisTask / sizeOuterHalfColumn;
            thisIndexInOuterColumn = thisTask&(sizeOuterHalfColumn-1); // thisTask % sizeOuterHalfColumn
            thisInnerBlock = thisIndexInOuterColumn/sizeInnerHalfBlock;
            // get index in state vector corresponding to upper inner block
            thisIndexInInnerBlock = thisTask&(sizeInnerHalfBlock-1); // thisTask % sizeInnerHalfBlock
            thisIndex = thisOuterColumn*sizeOuterColumn + thisInnerBlock*sizeInnerBlock
                + thisIndexInInnerBlock;
            // check if we are in the upper or lower half of an outer block
            outerBit = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBit*(sizeInnerHalfBlock);

            // NOTE: at this point thisIndex should be the index of the element we want to 
            // dephase in the chunk of the state vector on this process, in the 
            // density matrix representation. 
            // thisTask is the index of the pair element in pairStateVec
            // we will populate the second half of pairStateVec with this process'
            // data to send

            qureg.pairStateVec.real[thisTask+numTasks] = qureg.stateVec.real[thisIndex];
            qureg.pairStateVec.imag[thisTask+numTasks] = qureg.stateVec.imag[thisIndex];

        }
    }
}

void compressPairVectorForTwoQubitDepolarise(Qureg qureg, int targetQubit,
        int qubit2) {

    long long int sizeInnerBlockQ1, sizeInnerHalfBlockQ1;
    long long int sizeInnerBlockQ2, sizeInnerHalfBlockQ2, sizeInnerQuarterBlockQ2;
    long long int sizeOuterColumn, sizeOuterQuarterColumn;
    long long int 
         thisInnerBlockQ2,
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInOuterColumn,
         thisIndexInInnerBlockQ1,
         thisIndexInInnerBlockQ2,
         thisInnerBlockQ1InInnerBlockQ2;
    int outerBitQ1, outerBitQ2;

    long long int thisTask;
    long long int numTasks=qureg.numAmpsPerChunk>>2;

    // set dimensions
    sizeInnerHalfBlockQ1 = 1LL << targetQubit;
    sizeInnerHalfBlockQ2 = 1LL << qubit2;
    sizeInnerQuarterBlockQ2 = sizeInnerHalfBlockQ2 >> 1;
    sizeInnerBlockQ2 = sizeInnerHalfBlockQ2 << 1;
    sizeInnerBlockQ1 = 2LL * sizeInnerHalfBlockQ1;
    sizeOuterColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterQuarterColumn = sizeOuterColumn >> 2;
 
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlockQ1,sizeInnerHalfBlockQ1,sizeInnerQuarterBlockQ2,sizeInnerHalfBlockQ2,sizeInnerBlockQ2, \
                sizeOuterColumn,sizeOuterQuarterColumn,qureg,numTasks,targetQubit,qubit2) \
    private  (thisTask,thisInnerBlockQ2,thisOuterColumn,thisIndex,thisIndexInOuterColumn, \
                thisIndexInInnerBlockQ1,thisIndexInInnerBlockQ2,thisInnerBlockQ1InInnerBlockQ2,outerBitQ1,outerBitQ2) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // thisTask iterates over half the elements in this process' chunk of the density matrix
        // treat this as iterating over all columns, then iterating over half the values
        // within one column.
        // If this function has been called, this process' chunk contains half an 
        // outer block or less
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // we want to process all columns in the density matrix,
            // updating the values for half of each column (one half of each inner block)
            thisOuterColumn = thisTask / sizeOuterQuarterColumn;
            // thisTask % sizeOuterQuarterColumn
            thisIndexInOuterColumn = thisTask&(sizeOuterQuarterColumn-1);
            thisInnerBlockQ2 = thisIndexInOuterColumn / sizeInnerQuarterBlockQ2;
            // thisTask % sizeInnerQuarterBlockQ2;
            thisIndexInInnerBlockQ2 = thisTask&(sizeInnerQuarterBlockQ2-1);
            thisInnerBlockQ1InInnerBlockQ2 = thisIndexInInnerBlockQ2 / sizeInnerHalfBlockQ1;
            // thisTask % sizeInnerHalfBlockQ1;
            thisIndexInInnerBlockQ1 = thisTask&(sizeInnerHalfBlockQ1-1);

            // get index in state vector corresponding to upper inner block
            thisIndex = thisOuterColumn*sizeOuterColumn + thisInnerBlockQ2*sizeInnerBlockQ2
                + thisInnerBlockQ1InInnerBlockQ2*sizeInnerBlockQ1 + thisIndexInInnerBlockQ1;

            // check if we are in the upper or lower half of an outer block for Q1
            outerBitQ1 = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBitQ1*(sizeInnerHalfBlockQ1);

            // check if we are in the upper or lower half of an outer block for Q2
            outerBitQ2 = extractBit(qubit2, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBitQ2*(sizeInnerQuarterBlockQ2<<1);

            // NOTE: at this point thisIndex should be the index of the element we want to 
            // dephase in the chunk of the state vector on this process, in the 
            // density matrix representation. 
            // thisTask is the index of the pair element in pairStateVec

            // state[thisIndex] = (1-depolLevel)*state[thisIndex] + depolLevel*(state[thisIndex]
            //      + pair[thisTask])/2
            qureg.pairStateVec.real[thisTask+numTasks*2] = qureg.stateVec.real[thisIndex];
            qureg.pairStateVec.imag[thisTask+numTasks*2] = qureg.stateVec.imag[thisIndex];
        }
    }
}


void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depolLevel) {
    if (depolLevel == 0)
        return;
    
    int rankIsUpper; // rank is in the upper half of an outer block
    int pairRank; // rank of corresponding chunk

    int useLocalDataOnly = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
            qureg.numQubitsRepresented, targetQubit);

    if (useLocalDataOnly){
        densmatr_mixDepolarisingLocal(qureg, targetQubit, depolLevel);
    } else {
        // pack data to send to my pair process into the first half of pairStateVec
        compressPairVectorForSingleQubitDepolarise(qureg, targetQubit);

        rankIsUpper = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit, 
                qureg.numQubitsRepresented);
        pairRank = getChunkOuterBlockPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, 
                targetQubit, qureg.numQubitsRepresented);

        exchangePairStateVectorHalves(qureg, pairRank);
        densmatr_mixDepolarisingDistributed(qureg, targetQubit, depolLevel);
    }

}

void densmatr_mixDamping(Qureg qureg, int targetQubit, qreal damping) {
    if (damping == 0)
        return;
    
    int rankIsUpper; // rank is in the upper half of an outer block
    int pairRank; // rank of corresponding chunk

    int useLocalDataOnly = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
            qureg.numQubitsRepresented, targetQubit);

    if (useLocalDataOnly){
        densmatr_mixDampingLocal(qureg, targetQubit, damping);
    } else {
        // pack data to send to my pair process into the first half of pairStateVec
        compressPairVectorForSingleQubitDepolarise(qureg, targetQubit);

        rankIsUpper = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit, 
                qureg.numQubitsRepresented);
        pairRank = getChunkOuterBlockPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, 
                targetQubit, qureg.numQubitsRepresented);

        exchangePairStateVectorHalves(qureg, pairRank);
        densmatr_mixDampingDistributed(qureg, targetQubit, damping);
    }

}

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel){
    if (depolLevel == 0)
        return;
    int rankIsUpperBiggerQubit, rankIsUpperSmallerQubit;
    int pairRank; // rank of corresponding chunk
    int biggerQubit, smallerQubit;

    densmatr_mixTwoQubitDephasing(qureg, qubit1, qubit2, depolLevel);
    
    qreal eta = 2/depolLevel;
    qreal delta = eta - 1 - sqrt( (eta-1)*(eta-1) - 1 ); 
    qreal gamma = 1+delta;
    gamma = 1/(gamma*gamma*gamma);
    qreal GAMMA_PARTS_1_OR_2 = 1.0;
    // TODO -- test delta too small
    /*   
    if (fabs(4*delta*(1+delta)*gamma-depolLevel)>1e-5){
        printf("Numerical error in delta; for small error rates try Taylor expansion.\n");
        exit(1);
    }
    */
    
    biggerQubit = qubit1 > qubit2 ? qubit1 : qubit2;
    smallerQubit = qubit1 < qubit2 ? qubit1 : qubit2;
    int useLocalDataOnlyBigQubit, useLocalDataOnlySmallQubit;

    useLocalDataOnlyBigQubit = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
        qureg.numQubitsRepresented, biggerQubit);
    if (useLocalDataOnlyBigQubit){
        // does parts 1, 2 and 3 locally in one go
        densmatr_mixTwoQubitDepolarisingLocal(qureg, qubit1, qubit2, delta, gamma);
    } else {
        useLocalDataOnlySmallQubit = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
            qureg.numQubitsRepresented, smallerQubit);
        if (useLocalDataOnlySmallQubit){
            // do part 1 locally
            densmatr_mixTwoQubitDepolarisingLocalPart1(qureg, smallerQubit, biggerQubit, delta);
            
            // do parts 2 and 3 distributed (if part 2 is distributed part 3 is also distributed)
            // part 2 will be distributed and the value of the small qubit won't matter
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperBiggerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, biggerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperBiggerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    biggerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_mixTwoQubitDepolarisingDistributed(qureg, smallerQubit, biggerQubit, delta, GAMMA_PARTS_1_OR_2);
            
            // part 3 will be distributed but involve rearranging for the smaller qubit
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperBiggerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, biggerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperBiggerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    biggerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_mixTwoQubitDepolarisingQ1LocalQ2DistributedPart3(qureg, smallerQubit, biggerQubit, delta, gamma);
        } else {
            // do part 1, 2 and 3 distributed
            // part 1
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperSmallerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, smallerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperSmallerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    smallerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_mixTwoQubitDepolarisingDistributed(qureg, smallerQubit, biggerQubit, delta, GAMMA_PARTS_1_OR_2);

            // part 2
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperBiggerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, biggerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperBiggerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    biggerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_mixTwoQubitDepolarisingDistributed(qureg, smallerQubit, biggerQubit, delta, GAMMA_PARTS_1_OR_2);

            // part 3
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            pairRank = getChunkOuterBlockPairIdForPart3(rankIsUpperSmallerQubit, rankIsUpperBiggerQubit, 
                    qureg.chunkId, qureg.numAmpsPerChunk, smallerQubit, biggerQubit, qureg.numQubitsRepresented);
            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_mixTwoQubitDepolarisingDistributed(qureg, smallerQubit, biggerQubit, delta, gamma);

        }
    }

}

void statevec_compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta)
{
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
            statevec_compactUnitaryDistributed(qureg,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_compactUnitaryDistributed(qureg,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void statevec_unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u)
{
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
            statevec_unitaryDistributed(qureg,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_unitaryDistributed(qureg,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }


}

void statevec_controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta)
{
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
            statevec_controlledCompactUnitaryDistributed(qureg,controlQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_controlledCompactUnitaryDistributed(qureg,controlQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void statevec_controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, 
        ComplexMatrix2 u)
{
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
            statevec_controlledUnitaryDistributed(qureg,controlQubit,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_controlledUnitaryDistributed(qureg,controlQubit,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}

void statevec_multiControlledUnitary(Qureg qureg, long long int ctrlQubitsMask, long long int ctrlFlipMask, int targetQubit, ComplexMatrix2 u)
{
    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    Complex rot1, rot2;

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_multiControlledUnitaryLocal(qureg, targetQubit, ctrlQubitsMask, ctrlFlipMask, u);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        getRotAngleFromUnitaryMatrix(rankIsUpper, &rot1, &rot2, u);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);

        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);

        // this rank's values are either in the upper of lower half of the block. send values to multiControlledUnitaryDistributed
        // in the correct order
        if (rankIsUpper){
            statevec_multiControlledUnitaryDistributed(qureg,targetQubit,ctrlQubitsMask,ctrlFlipMask,rot1,rot2,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec); //output
        } else {
            statevec_multiControlledUnitaryDistributed(qureg,targetQubit,ctrlQubitsMask,ctrlFlipMask,rot1,rot2,
                    qureg.pairStateVec, //upper
                    qureg.stateVec, //lower
                    qureg.stateVec); //output
        }
    }
}
void statevec_pauliX(Qureg qureg, int targetQubit)
{
    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);

    // rank's chunk is in upper half of block 
    int rankIsUpper;
    int pairRank; // rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_pauliXLocal(qureg, targetQubit);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        //printf("%d rank has pair rank: %d\n", qureg.rank, pairRank);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block. pauliX just replaces
        // this rank's values with pair values
        statevec_pauliXDistributed(qureg,
                qureg.pairStateVec, // in
                qureg.stateVec); // out
    }
}

void statevec_controlledNot(Qureg qureg, int controlQubit, int targetQubit)
{
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
        statevec_controlledNotDistributed(qureg,controlQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec); //out
    }
}

void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask)
{   
    /* operation is the same regardless of control and target ordering, hence 
     * we accept only bitMasks (for convenience of caller, when shifting qubits
     * for density matrices)
     */
    
    // global index of the first basis state in this node
    long long int firstInd = qureg.chunkId * qureg.numAmpsPerChunk;
    
    /* optimisation: if this node doesn't contain any amplitudes for which {ctrls}={1}
     * (and hence, neither does the pair node), then these pair nodes have nothing to modify 
     * nor any need to communicate, and can halt. No ctrls are contained in the node 
     * if the distance from the first index, to the next index where {ctrls}=1, is 
     * greater than the total contained amplitudes. This is a worthwhile optimisation, 
     * since although we must still wait for the slowest node, we have potentially reduced 
     * the network traffic and might avoid saturation. 
     */
    if ((firstInd|ctrlMask) - firstInd >= qureg.numAmpsPerChunk)
        return;
        
    /* nodes communicate pairwise, and (ignoring ctrls) swap all their amplitudes with mate.
     * hence we find |pairState> = X_{targs}|firstStateInNode>, determine which node contains  
     * |pairState>, and swap state-vector with it (unless it happens to be this node).
     */
    
    // global index of the corresponding NOT'd first basis state
    long long int pairInd = firstInd ^ targMask;
    int pairRank = pairInd / qureg.numAmpsPerChunk;
    int useLocalDataOnly = (pairRank == qureg.chunkId);
    
    if (useLocalDataOnly) {
        // swaps amplitudes locally, setting |a>=X|b>, and |b>=X|a>
        statevec_multiControlledMultiQubitNotLocal(qureg, ctrlMask, targMask);
    } else {
        // swaps amplitudes with pair node
        exchangeStateVectors(qureg, pairRank);
        //  modifies only |a>=X|b> (pair node handles the converse)
        statevec_multiControlledMultiQubitNotDistributed(
            qureg, ctrlMask, targMask,
            qureg.pairStateVec, // in
            qureg.stateVec);    // out
    }
}

void statevec_pauliY(Qureg qureg, int targetQubit)
{	
	int conjFac = 1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper;	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        statevec_pauliYLocal(qureg, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        statevec_pauliYDistributed(qureg,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper, conjFac);
    }
}

void statevec_pauliYConj(Qureg qureg, int targetQubit)
{	
	int conjFac = -1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper;	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        statevec_pauliYLocal(qureg, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        statevec_pauliYDistributed(qureg,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper, conjFac);
    }
}

void statevec_controlledPauliY(Qureg qureg, int controlQubit, int targetQubit)
{
	int conjFac = 1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper; 	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledPauliYLocal(qureg, controlQubit, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        if (rankIsUpper){
            statevec_controlledPauliYDistributed(qureg,controlQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        } else {
            statevec_controlledPauliYDistributed(qureg,controlQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					- conjFac); //out
        }
    }
}

void statevec_controlledPauliYConj(Qureg qureg, int controlQubit, int targetQubit)
{
	int conjFac = -1;

    // flag to require memory exchange. 1: an entire block fits on one rank, 0: at most half a block fits on one rank
    int useLocalDataOnly = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targetQubit);
    int rankIsUpper; 	// rank's chunk is in upper half of block 
    int pairRank; 		// rank of corresponding chunk

    if (useLocalDataOnly){
        // all values required to update state vector lie in this rank
        statevec_controlledPauliYLocal(qureg, controlQubit, targetQubit, conjFac);
    } else {
        // need to get corresponding chunk of state vector from other rank
        rankIsUpper = chunkIsUpper(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        pairRank = getChunkPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, targetQubit);
        // get corresponding values from my pair
        exchangeStateVectors(qureg, pairRank);
        // this rank's values are either in the upper of lower half of the block
        if (rankIsUpper){
            statevec_controlledPauliYDistributed(qureg,controlQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        } else {
            statevec_controlledPauliYDistributed(qureg,controlQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					- conjFac); //out
        }
    }
}

void statevec_hadamard(Qureg qureg, int targetQubit)
{
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
            statevec_hadamardDistributed(qureg,
                    qureg.stateVec, //upper
                    qureg.pairStateVec, //lower
                    qureg.stateVec, rankIsUpper); //output
        } else {
            statevec_hadamardDistributed(qureg,
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
 * @param[in] measureQubit qubit being measured
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

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    qreal stateProb=0, totalStateProb=0;
    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        stateProb = statevec_findProbabilityOfZeroLocal(qureg, measureQubit);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            stateProb = statevec_findProbabilityOfZeroDistributed(qureg);
        } else stateProb = 0;
    }
    MPI_Allreduce(&stateProb, &totalStateProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    if (outcome==1) totalStateProb = 1.0 - totalStateProb;
    return totalStateProb;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome) {
	
	qreal zeroProb = densmatr_findProbabilityOfZeroLocal(qureg, measureQubit);
	
	qreal outcomeProb;
	MPI_Allreduce(&zeroProb, &outcomeProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
	if (outcome == 1)
		outcomeProb = 1.0 - outcomeProb;
	
	return outcomeProb;
}

void statevec_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits) {
    
    // each node populates retProbs with contributions from the subset of amps in each
    statevec_calcProbOfAllOutcomesLocal(retProbs, qureg, qubits, numQubits);

    // then, retProbs are summed element-wise
    MPI_Allreduce(MPI_IN_PLACE, retProbs, 1LL<<numQubits, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
}

void densmatr_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits) {
    
    // each node populates retProbs with contributions from the subset of amps in each
    densmatr_calcProbOfAllOutcomesLocal(retProbs, qureg, qubits, numQubits);

    // then, retProbs are summed element-wise
    MPI_Allreduce(MPI_IN_PLACE, retProbs, 1LL<<numQubits, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
}

qreal densmatr_calcPurity(Qureg qureg) {
    
    qreal localPurity = densmatr_calcPurityLocal(qureg);
        
    qreal globalPurity;
    MPI_Allreduce(&localPurity, &globalPurity, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    return globalPurity;
}

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal totalStateProb)
{
    int skipValuesWithinRank = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, measureQubit);
    if (skipValuesWithinRank) {
        statevec_collapseToKnownProbOutcomeLocal(qureg, measureQubit, outcome, totalStateProb);
    } else {
        if (!isChunkToSkipInFindPZero(qureg.chunkId, qureg.numAmpsPerChunk, measureQubit)){
            // chunk has amps for q=0
            if (outcome==0) statevec_collapseToKnownProbOutcomeDistributedRenorm(qureg, measureQubit, 
                    totalStateProb);
            else statevec_collapseToOutcomeDistributedSetZero(qureg);
        } else {
            // chunk has amps for q=1
            if (outcome==1) statevec_collapseToKnownProbOutcomeDistributedRenorm(qureg, measureQubit, 
                    totalStateProb);
            else statevec_collapseToOutcomeDistributedSetZero(qureg);
        }
    }
}

void seedQuEST(QuESTEnv *env, unsigned long int* seedArray, int numSeeds) {

    // it is imperative every node agrees on the seed, so that random decisions 
    // agree on every node. Hence we use only the master node keys.
    MPI_Bcast(seedArray, numSeeds, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    
    // free existing seed array, if exists
    if (env->seeds != NULL)
        free(env->seeds);
        
    // record keys in permanent heap
    env->seeds = malloc(numSeeds * sizeof *(env->seeds));
    for (int i=0; i<numSeeds; i++)
        (env->seeds)[i] = seedArray[i];
    env->numSeeds = numSeeds;
    
    // pass keys to Mersenne Twister seeder
    init_by_array(seedArray, numSeeds); 
}

/** returns -1 if this node contains no amplitudes where qb1 and qb2 
 * have opposite parity, otherwise returns the global index of one
 * of such contained amplitudes (not necessarily the first)
 */
long long int getGlobalIndOfOddParityInChunk(Qureg qureg, int qb1, int qb2) {
    long long int chunkStartInd = qureg.numAmpsPerChunk * qureg.chunkId;
    long long int chunkEndInd = chunkStartInd + qureg.numAmpsPerChunk; // exclusive
    long long int oddParityInd;
    
    if (extractBit(qb1, chunkStartInd) != extractBit(qb2, chunkStartInd))
        return chunkStartInd;
    
    oddParityInd = flipBit(chunkStartInd, qb1);
    if (oddParityInd >= chunkStartInd && oddParityInd < chunkEndInd)
        return oddParityInd;
        
    oddParityInd = flipBit(chunkStartInd, qb2);
    if (oddParityInd >= chunkStartInd && oddParityInd < chunkEndInd)
        return oddParityInd;
    
    return -1;
}

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) {
    
    // perform locally if possible 
    int qbBig = (qb1 > qb2)? qb1 : qb2;
    if (halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, qbBig))
        return statevec_swapQubitAmpsLocal(qureg, qb1, qb2);
        
    // do nothing if this node contains no amplitudes to swap
    long long int oddParityGlobalInd = getGlobalIndOfOddParityInChunk(qureg, qb1, qb2);
    if (oddParityGlobalInd == -1)
        return;

    // determine and swap amps with pair node
    int pairRank = flipBit(flipBit(oddParityGlobalInd, qb1), qb2) / qureg.numAmpsPerChunk;
    exchangeStateVectors(qureg, pairRank);
    statevec_swapQubitAmpsDistributed(qureg, pairRank, qb1, qb2);
}

/** This calls swapQubitAmps only when it would involve a distributed communication;
 * if the qubit chunks already fit in the node, it operates the unitary direct.
 * Note the order of q1 and q2 in the call to twoQubitUnitaryLocal is important.
 * 
 * @todo refactor so that the 'swap back' isn't performed; instead the qubit locations 
 * are updated.
 * @todo the double swap (q1,q2 to 0,1) may be possible simultaneously by a bespoke 
 * swap routine.
 */
void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int q1, int q2, ComplexMatrix4 u) {
    int q1FitsInNode = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, q1);
    int q2FitsInNode = halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, q2);
        
    if (q1FitsInNode && q2FitsInNode) {    
        statevec_multiControlledTwoQubitUnitaryLocal(qureg, ctrlMask, q1, q2, u);
        
    } else if (q1FitsInNode) {
        int qSwap = (q1 > 0)? q1-1 : q1+1;

        // ensure ctrl == qSwap, ensure ctrlMask updates under the swap
        if (maskContainsBit(ctrlMask, qSwap))
            ctrlMask = flipBit(flipBit(ctrlMask, q2), qSwap);

        statevec_swapQubitAmps(qureg, q2, qSwap);
        statevec_multiControlledTwoQubitUnitaryLocal(qureg, ctrlMask, q1, qSwap, u);
        statevec_swapQubitAmps(qureg, q2, qSwap);

    } else if (q2FitsInNode) {
        int qSwap = (q2 > 0)? q2-1 : q2+1;
        
        // ensure ctrl == qSwap, ensure ctrlMask updates under the swap
        if (maskContainsBit(ctrlMask, qSwap))
            ctrlMask = flipBit(flipBit(ctrlMask, q1), qSwap);
        
        statevec_swapQubitAmps(qureg, q1, qSwap);
        statevec_multiControlledTwoQubitUnitaryLocal(qureg, ctrlMask, qSwap, q2, u);
        statevec_swapQubitAmps(qureg, q1, qSwap);
        
    } else {
        // we know with certainty that both q1 and q2 >= 2
        int swap1 = 0;
        int swap2 = 1;
        
        // if ctrl == swap1 or swap2, ensure ctrlMask updates under the swap
        if (maskContainsBit(ctrlMask, swap1))
            ctrlMask = flipBit(flipBit(ctrlMask, swap1), q1);
        if (maskContainsBit(ctrlMask, swap2))
            ctrlMask = flipBit(flipBit(ctrlMask, swap2), q2);
        
        statevec_swapQubitAmps(qureg, q1, swap1);
        statevec_swapQubitAmps(qureg, q2, swap2);
        statevec_multiControlledTwoQubitUnitaryLocal(qureg, ctrlMask, swap1, swap2, u);
        statevec_swapQubitAmps(qureg, q1, swap1);
        statevec_swapQubitAmps(qureg, q2, swap2);
    }
}

/** This calls swapQubitAmps only when it would involve a distributed communication;
 * if the qubit chunks already fit in the node, it operates the unitary direct.
 * It is already gauranteed here that all target qubits can fit on each node (this is
 * validated in the front-end)
 * 
 * @todo refactor so that the 'swap back' isn't performed; instead the qubit locations 
 * are updated.
 */
void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u) {

    // bit mask of target qubits (for quick collision checking)
    long long int targMask = getQubitBitMask(targs, numTargs);
    
    // find lowest qubit available for swapping (isn't in targs)
    int freeQb=0;
    while (maskContainsBit(targMask, freeQb))
        freeQb++;
        
    // assign indices of where each target will be swapped to (else itself)
    int swapTargs[numTargs];
    for (int t=0; t<numTargs; t++) {
        if (halfMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, targs[t]))
            swapTargs[t] = targs[t];
        else {
            // mark swap
            swapTargs[t] = freeQb;
            
            // update ctrlMask if swapped-out qubit was a control
            if (maskContainsBit(ctrlMask, swapTargs[t]))
                ctrlMask = flipBit(flipBit(ctrlMask, swapTargs[t]), targs[t]); // swap targ and ctrl
            
            // locate next available on-chunk qubit
            freeQb++;
            while (maskContainsBit(targMask, freeQb))
                freeQb++;
        }
    }
    
    // perform swaps as necessary 
    for (int t=0; t<numTargs; t++)
        if (swapTargs[t] != targs[t])
            statevec_swapQubitAmps(qureg, targs[t], swapTargs[t]);
    
    // all target qubits have now been swapped into local memory
    statevec_multiControlledMultiQubitUnitaryLocal(qureg, ctrlMask, swapTargs, numTargs, u);
    
    // undo swaps 
    for (int t=0; t<numTargs; t++)
        if (swapTargs[t] != targs[t])
            statevec_swapQubitAmps(qureg, targs[t], swapTargs[t]);
}


void copyDiagOpIntoMatrixPairState(Qureg qureg, DiagonalOp op) {
        
    /* since, for every elem in 2^N op, there is a column in 2^N x 2^N qureg, 
     * we know immediately (by each node containing at least 1 element of op)
     * that every node contains at least 1 column. Hence, we know that pairStateVec 
     * of qureg can fit the entirety of op.
     */

    // load up our local contribution
    long long int localOffset = qureg.chunkId * op.numElemsPerChunk;
    memcpy(&qureg.pairStateVec.real[localOffset], op.real, op.numElemsPerChunk * sizeof(qreal));
    memcpy(&qureg.pairStateVec.imag[localOffset], op.imag, op.numElemsPerChunk * sizeof(qreal));
    
    // work out how many messages are needed to send op chunks (2GB limit)
    long long int maxMsgSize = MPI_MAX_AMPS_IN_MSG;
    if (op.numElemsPerChunk < maxMsgSize) 
        maxMsgSize = op.numElemsPerChunk;
    int numMsgs = op.numElemsPerChunk / maxMsgSize; // since MPI_MAX... = 2^n, division is exact
    
    // each node has a turn at broadcasting its contribution of op
    for (int broadcaster=0; broadcaster < qureg.numChunks; broadcaster++) {
        long long int broadOffset = broadcaster * op.numElemsPerChunk;
    
        // (while keeping each message smaller than MPI max)
        for (int i=0; i<numMsgs; i++) {
            MPI_Bcast(
                &qureg.pairStateVec.real[broadOffset + i*maxMsgSize], 
                maxMsgSize,  MPI_QuEST_REAL, broadcaster, MPI_COMM_WORLD);
            MPI_Bcast(
                &qureg.pairStateVec.imag[broadOffset + i*maxMsgSize], 
                maxMsgSize,  MPI_QuEST_REAL, broadcaster, MPI_COMM_WORLD);
        }
    }
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op) {
    
    copyDiagOpIntoMatrixPairState(qureg, op);
    densmatr_applyDiagonalOpLocal(qureg, op);
}

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {

    Complex localExpec = statevec_calcExpecDiagonalOpLocal(qureg, op);
    if (qureg.numChunks == 1)
        return localExpec;
        
    qreal localReal = localExpec.real;
    qreal localImag = localExpec.imag;
    qreal globalReal, globalImag;
    MPI_Allreduce(&localReal, &globalReal, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localImag, &globalImag, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    Complex globalExpec;
    globalExpec.real = globalReal;
    globalExpec.imag = globalImag;
    return globalExpec;
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {
    
    Complex localVal = densmatr_calcExpecDiagonalOpLocal(qureg, op);
    if (qureg.numChunks == 1)
        return localVal;
    
    qreal localRe = localVal.real;
    qreal localIm = localVal.imag;
    qreal globalRe, globalIm;
    
    MPI_Allreduce(&localRe, &globalRe, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localIm, &globalIm, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    Complex globalVal;
    globalVal.real = globalRe;
    globalVal.imag = globalIm;
    return globalVal;
}
