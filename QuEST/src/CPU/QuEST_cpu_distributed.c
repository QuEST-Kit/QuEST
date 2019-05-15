// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the backend in ../QuEST_ops.h for an MPI environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 */

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "mt19937ar.h"

# include "QuEST_cpu_internal.h"

#define _BSD_SOURCE
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

/** Get the value of the bit at a particular index in a number.
  SCB edit: new definition of extractBit is much faster ***
 * @param[in] locationOfBitFromRight location of bit in theEncodedNumber
 * @param[in] theEncodedNumber number to search
 * @return the value of the bit in theEncodedNumber
 */
// TODO -- This should not also be multiply defined in QuEST_cpu.c -- move to a helper functions file
static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber)
{
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

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
	
	// computes first local index containing a diagonal element
	long long int diagSpacing = 1LL + (1LL << qureg.numQubitsRepresented);
    long long int numPrevDiags = (qureg.chunkId>0)? 1+(qureg.chunkId*qureg.numAmpsPerChunk)/diagSpacing : 0;
	long long int globalIndNextDiag = diagSpacing * numPrevDiags;
	long long int localIndNextDiag = globalIndNextDiag % qureg.numAmpsPerChunk;
	long long int index;
	
	qreal rankTotal = 0;
	qreal y, t, c;
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
static int getChunkOuterBlockPairId(int chunkIsUpper, int chunkId, long long int chunkSize, int targetQubit,        int numQubits);
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
    
	seedQuESTDefault();
    
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

void copyVecIntoMatrixPairState(Qureg matr, Qureg vec) {
    
    // copy this state's pure state section into this qureg's pairState
    long long int numLocalAmps = vec.numAmpsPerChunk;
    long long int myOffset = vec.chunkId * numLocalAmps;
    memcpy(&matr.pairStateVec.real[myOffset], vec.stateVec.real, numLocalAmps * sizeof(qreal) );
    memcpy(&matr.pairStateVec.imag[myOffset], vec.stateVec.imag, numLocalAmps * sizeof(qreal) );

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
    int TAG=100;
    MPI_Status status;

    // Multiple messages are required as MPI uses int rather than long long int for count
    // For openmpi, messages are further restricted to 2GB in size -- do this for all cases
    // to be safe
    long long int maxMessageCount = MPI_MAX_AMPS_IN_MSG;
    if (qureg.numAmpsPerChunk < maxMessageCount) 
        maxMessageCount = qureg.numAmpsPerChunk;
    
    // safely assume MPI_MAX... = 2^n, so division always exact
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

//TODO -- decide where this function should go. It is a preparation for MPI data transfer function
void compressPairVectorForSingleQubitDepolarise(Qureg qureg, const int targetQubit){
    long long int sizeInnerBlock, sizeInnerHalfBlock;
    long long int sizeOuterColumn, sizeOuterHalfColumn;
    long long int thisInnerBlock, // current block
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInOuterColumn,
         thisIndexInInnerBlock;
         
    int outerBit;

    long long int thisTask;
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeInnerHalfBlock = 1LL << targetQubit;
    sizeInnerBlock     = 2LL * sizeInnerHalfBlock;
    sizeOuterHalfColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterColumn     = 2LL * sizeOuterHalfColumn;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlock,sizeInnerHalfBlock,sizeOuterColumn,sizeOuterHalfColumn,qureg) \
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

void compressPairVectorForTwoQubitDepolarise(Qureg qureg, const int targetQubit,
        const int qubit2) {

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
    const long long int numTasks=qureg.numAmpsPerChunk>>2;

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
                sizeOuterColumn, \
                sizeOuterQuarterColumn,qureg) \
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


void densmatr_oneQubitDepolarise(Qureg qureg, const int targetQubit, qreal depolLevel) {
    if (depolLevel == 0)
        return;
    
    int rankIsUpper; // rank is in the upper half of an outer block
    int pairRank; // rank of corresponding chunk

    int useLocalDataOnly = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
            qureg.numQubitsRepresented, targetQubit);

    if (useLocalDataOnly){
        densmatr_oneQubitDepolariseLocal(qureg, targetQubit, depolLevel);
    } else {
        // pack data to send to my pair process into the first half of pairStateVec
        compressPairVectorForSingleQubitDepolarise(qureg, targetQubit);

        rankIsUpper = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit, 
                qureg.numQubitsRepresented);
        pairRank = getChunkOuterBlockPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, 
                targetQubit, qureg.numQubitsRepresented);

        exchangePairStateVectorHalves(qureg, pairRank);
        densmatr_oneQubitDepolariseDistributed(qureg, targetQubit, depolLevel);
    }

}

void densmatr_oneQubitDamping(Qureg qureg, const int targetQubit, qreal damping) {
    if (damping == 0)
        return;
    
    int rankIsUpper; // rank is in the upper half of an outer block
    int pairRank; // rank of corresponding chunk

    int useLocalDataOnly = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
            qureg.numQubitsRepresented, targetQubit);

    if (useLocalDataOnly){
        densmatr_oneQubitDampingLocal(qureg, targetQubit, damping);
    } else {
        // pack data to send to my pair process into the first half of pairStateVec
        compressPairVectorForSingleQubitDepolarise(qureg, targetQubit);

        rankIsUpper = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, targetQubit, 
                qureg.numQubitsRepresented);
        pairRank = getChunkOuterBlockPairId(rankIsUpper, qureg.chunkId, qureg.numAmpsPerChunk, 
                targetQubit, qureg.numQubitsRepresented);

        exchangePairStateVectorHalves(qureg, pairRank);
        densmatr_oneQubitDampingDistributed(qureg, targetQubit, damping);
    }

}

void densmatr_twoQubitDepolarise(Qureg qureg, int qubit1, int qubit2, qreal depolLevel){
    if (depolLevel == 0)
        return;
    int rankIsUpperBiggerQubit, rankIsUpperSmallerQubit;
    int pairRank; // rank of corresponding chunk
    int biggerQubit, smallerQubit;

    densmatr_twoQubitDephase(qureg, qubit1, qubit2, depolLevel);
    
    qreal eta = 2/depolLevel;
    qreal delta = eta - 1 - sqrt( (eta-1)*(eta-1) - 1 ); 
    qreal gamma = 1+delta;
    gamma = 1/(gamma*gamma*gamma);
    const qreal GAMMA_PARTS_1_OR_2 = 1.0;
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
        densmatr_twoQubitDepolariseLocal(qureg, qubit1, qubit2, delta, gamma);
    } else {
        useLocalDataOnlySmallQubit = densityMatrixBlockFitsInChunk(qureg.numAmpsPerChunk, 
            qureg.numQubitsRepresented, smallerQubit);
        if (useLocalDataOnlySmallQubit){
            // do part 1 locally
            densmatr_twoQubitDepolariseLocalPart1(qureg, smallerQubit, biggerQubit, delta);
            
            // do parts 2 and 3 distributed (if part 2 is distributed part 3 is also distributed)
            // part 2 will be distributed and the value of the small qubit won't matter
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperBiggerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, biggerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperBiggerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    biggerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_twoQubitDepolariseDistributed(qureg, smallerQubit, biggerQubit, delta, GAMMA_PARTS_1_OR_2);
            
            // part 3 will be distributed but involve rearranging for the smaller qubit
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperBiggerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, biggerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperBiggerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    biggerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_twoQubitDepolariseQ1LocalQ2DistributedPart3(qureg, smallerQubit, biggerQubit, delta, gamma);
        } else {
            // do part 1, 2 and 3 distributed
            // part 1
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperSmallerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, smallerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperSmallerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    smallerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_twoQubitDepolariseDistributed(qureg, smallerQubit, biggerQubit, delta, GAMMA_PARTS_1_OR_2);

            // part 2
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            rankIsUpperBiggerQubit = chunkIsUpperInOuterBlock(qureg.chunkId, qureg.numAmpsPerChunk, biggerQubit, 
                    qureg.numQubitsRepresented);
            pairRank = getChunkOuterBlockPairId(rankIsUpperBiggerQubit, qureg.chunkId, qureg.numAmpsPerChunk, 
                    biggerQubit, qureg.numQubitsRepresented);

            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_twoQubitDepolariseDistributed(qureg, smallerQubit, biggerQubit, delta, GAMMA_PARTS_1_OR_2);

            // part 3
            compressPairVectorForTwoQubitDepolarise(qureg, smallerQubit, biggerQubit);
            pairRank = getChunkOuterBlockPairIdForPart3(rankIsUpperSmallerQubit, rankIsUpperBiggerQubit, 
                    qureg.chunkId, qureg.numAmpsPerChunk, smallerQubit, biggerQubit, qureg.numQubitsRepresented);
            exchangePairStateVectorHalves(qureg, pairRank);
            densmatr_twoQubitDepolariseDistributed(qureg, smallerQubit, biggerQubit, delta, gamma);

        }
    }

}

void statevec_compactUnitary(Qureg qureg, const int targetQubit, Complex alpha, Complex beta)
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

void statevec_unitary(Qureg qureg, const int targetQubit, ComplexMatrix2 u)
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

void statevec_controlledCompactUnitary(Qureg qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta)
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

void statevec_controlledUnitary(Qureg qureg, const int controlQubit, const int targetQubit, 
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

void statevec_multiControlledUnitary(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);

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
void statevec_pauliX(Qureg qureg, const int targetQubit)
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
        statevec_pauliXDistributed(qureg, targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec); // out
    }
}

void statevec_controlledNot(Qureg qureg, const int controlQubit, const int targetQubit)
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

void statevec_pauliY(Qureg qureg, const int targetQubit)
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
        statevec_pauliYDistributed(qureg,targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper, conjFac);
    }
}

void statevec_pauliYConj(Qureg qureg, const int targetQubit)
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
        statevec_pauliYDistributed(qureg,targetQubit,
                qureg.pairStateVec, // in
                qureg.stateVec, // out
                rankIsUpper, conjFac);
    }
}

void statevec_controlledPauliY(Qureg qureg, const int controlQubit, const int targetQubit)
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
            statevec_controlledPauliYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        } else {
            statevec_controlledPauliYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					-conjFac); //out
        }
    }
}

void statevec_controlledPauliYConj(Qureg qureg, const int controlQubit, const int targetQubit)
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
            statevec_controlledPauliYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					conjFac); //out
        } else {
            statevec_controlledPauliYDistributed(qureg,controlQubit,targetQubit,
                    qureg.pairStateVec, //in
                    qureg.stateVec,
					-conjFac); //out
        }
    }
}

void statevec_hadamard(Qureg qureg, const int targetQubit)
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

qreal statevec_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome)
{
    qreal stateProb=0, totalStateProb=0;
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

qreal densmatr_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome) {
	
	qreal zeroProb = densmatr_findProbabilityOfZeroLocal(qureg, measureQubit);
	
	qreal outcomeProb;
	MPI_Allreduce(&zeroProb, &outcomeProb, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
	if (outcome == 1)
		outcomeProb = 1.0 - outcomeProb;
	
	return outcomeProb;
}

qreal densmatr_calcPurity(Qureg qureg) {
    
    qreal localPurity = densmatr_calcPurityLocal(qureg);
        
    qreal globalPurity;
    MPI_Allreduce(&localPurity, &globalPurity, 1, MPI_QuEST_REAL, MPI_SUM, MPI_COMM_WORLD);
    
    return globalPurity;
}

void statevec_collapseToKnownProbOutcome(Qureg qureg, const int measureQubit, int outcome, qreal totalStateProb)
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

void seedQuESTDefault(){
    // init MT random number generator with three keys -- time and pid
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process

    unsigned long int key[2];
    getQuESTDefaultSeedKey(key);
    // this seed will be used to generate the same random number on all procs,
    // therefore we want to make sure all procs receive the same key
    MPI_Bcast(key, 2, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    init_by_array(key, 2);
}
