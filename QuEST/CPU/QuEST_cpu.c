// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file qubits.c
 * The core of the CPU backend functionality. The CPU/MPI implementations of the pure state functions in
 * ../QuEST_ops_pure.h are in QuEST_cpu_local.c and QuEST_cpu_distributed.c which mostly wrap the core
 * functions defined here. Some additional hardware-agnostic functions are defined here
 */

# include "../QuEST.h"
# include "../QuEST_internal.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"

# include "QuEST_cpu_internal.h"

# include <math.h>  
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>

# ifdef _OPENMP
# include <omp.h>
# endif


/** Get the value of the bit at a particular index in a number.
  SCB edit: new definition of extractBit is much faster ***
 * @param[in] locationOfBitFromRight location of bit in theEncodedNumber
 * @param[in] theEncodedNumber number to search
 * @return the value of the bit in theEncodedNumber
 */
static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber)
{
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}


void densmatr_initClassicalState (QubitRegister qureg, long long int stateInd)
{
    // dimension of the state vector
    long long int densityNumElems = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    REAL *densityReal = qureg.stateVec.real;
    REAL *densityImag = qureg.stateVec.imag;

    // initialise the state to all zeros
	long long int index;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (densityNumElems, densityReal, densityImag) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<densityNumElems; index++) {
            densityReal[index] = 0.0;
            densityImag[index] = 0.0;
        }
    }
	
	// index of the single density matrix elem to set non-zero
	long long int densityDim = 1LL << qureg.numDensityQubits;
	long long int densityInd = (densityDim + 1)*stateInd;

	// give the specified classical state prob 1
    if (qureg.chunkId == densityInd / densityDim){
        densityReal[densityInd % densityDim] = 1.0;
        densityImag[densityInd % densityDim] = 0.0;
    }
}


void densmatr_initStatePlus (QubitRegister qureg)
{
    long long int chunkSize, stateVecSize;
    long long int index;

    // dimension of the state vector
    chunkSize = qureg.numAmpsPerChunk;
    stateVecSize = chunkSize*qureg.numChunks;
    REAL probFactor = 1.0/((REAL)stateVecSize);

    // Can't use qureg->stateVec as a private OMP var
    REAL *densityReal = qureg.stateVec.real;
    REAL *densityImag = qureg.stateVec.imag;

    // initialise the state to |+++..+++> = 1/normFactor {1, 1, 1, ...}
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkSize, densityReal, densityImag, probFactor) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<chunkSize; index++) {
            densityReal[index] = probFactor;
            densityImag[index] = 0.0;
        }
    }
}

void densmatr_initPureStateLocal(QubitRegister targetQureg, QubitRegister copyQureg) {
	
	// targetQureg is a density matrix of size (2^N)^2, copy is pure of size 2^N

    // dimension of the pure state vector
    long long int stateVecSize = copyQureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    REAL *targetDensityReal = targetQureg.stateVec.real;
    REAL *targetDensityImag = targetQureg.stateVec.imag;
    REAL *copyStateVecReal = copyQureg.stateVec.real;
    REAL *copyStateVecImag = copyQureg.stateVec.imag;
	
	// iterates density entries
	long long int row, col;
	
	// involved elements
	REAL realRow, imagRow, realCol, imagCol;

    // initialise targetQureg to be 100% likely in the puerstate of copyQureg
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecSize, targetDensityReal, targetDensityImag, copyStateVecReal, copyStateVecImag) \
    private  (row, col, realRow, imagRow, realCol, imagCol) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (row=0; row < stateVecSize; row++) {
			
			realRow = copyStateVecReal[row];
			imagRow = copyStateVecImag[row];
			
			for (col=0; col < stateVecSize; col++) {

				realCol =   copyStateVecReal[col];
				imagCol = - copyStateVecImag[col]; //note minus for conjugation
			
				targetDensityReal[col*stateVecSize + row] = realRow*realCol - imagRow*imagCol;
				targetDensityImag[col*stateVecSize + row] = realRow*imagCol + imagRow*realCol;
			}
        }
    }
}




// @TODO
void densmatr_initPureStateDistributed(QubitRegister targetQureg, QubitRegister copyQureg) {
	QuESTAssert(0, 0, __func__);
}


void statevec_createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env)
{
    QuESTAssert(numQubits>0, 9, __func__);
    long long int numAmps = 1L << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;

    qureg->stateVec.real = malloc(numAmpsPerRank * sizeof(*(qureg->stateVec.real)));
    qureg->stateVec.imag = malloc(numAmpsPerRank * sizeof(*(qureg->stateVec.imag)));
    if (env.numRanks>1){
        qureg->pairStateVec.real = malloc(numAmpsPerRank * sizeof(*(qureg->pairStateVec.real)));
        qureg->pairStateVec.imag = malloc(numAmpsPerRank * sizeof(*(qureg->pairStateVec.imag)));
    }

    if ( (!(qureg->stateVec.real) || !(qureg->stateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!");
        exit (EXIT_FAILURE);
    }

    if ( env.numRanks>1 && (!(qureg->pairStateVec.real) || !(qureg->pairStateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!");
        exit (EXIT_FAILURE);
    }

    qureg->numQubits = numQubits;
	qureg->numAmpsTotal = numAmps;
    qureg->numAmpsPerChunk = numAmpsPerRank;
    qureg->chunkId = env.rank;
    qureg->numChunks = env.numRanks;
	qureg->isDensityMatrix = 0;
}

void statevec_destroyQubitRegister(QubitRegister qureg, QuESTEnv env){
    free(qureg.stateVec.real);
    free(qureg.stateVec.imag);
    if (env.numRanks>1){
        free(qureg.pairStateVec.real);
        free(qureg.pairStateVec.imag);
    }
}

void statevec_reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank){
    long long int index;
    int rank;
    if (qureg.numQubits<=5){
        for (rank=0; rank<qureg.numChunks; rank++){
            if (qureg.chunkId==rank){
                if (reportRank) {
                    printf("Reporting state from rank %d [\n", qureg.chunkId);
                    printf("real, imag\n");
                } else if (rank==0) {
                    printf("Reporting state [\n");
                    printf("real, imag\n");
                }

                for(index=0; index<qureg.numAmpsPerChunk; index++){
                    printf(REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", qureg.stateVec.real[index], qureg.stateVec.imag[index]);
                }
                if (reportRank || rank==qureg.numChunks-1) printf("]\n");
            }
            syncQuESTEnv(env);
        }
    } else printf("Error: reportStateToScreen will not print output for systems of more than 5 qubits.\n");
}

void statevec_getEnvironmentString(QuESTEnv env, QubitRegister qureg, char str[200]){
    int numThreads=1;
# ifdef _OPENMP
    numThreads=omp_get_max_threads(); 
# endif
    sprintf(str, "%dqubits_CPU_%dranksx%dthreads", qureg.numQubits, env.numRanks, numThreads);
}

void statevec_initStateZero (QubitRegister qureg)
{
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

    // initialise the state-vector to all-zeroes
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecSize, stateVecReal, stateVecImag) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            stateVecReal[index] = 0.0;
            stateVecImag[index] = 0.0;
        }
    }

    if (qureg.chunkId==0){
        // zero state |0000..0000> has probability 1
        stateVecReal[0] = 1.0;
        stateVecImag[0] = 0.0;
    }
}

void statevec_initStatePlus (QubitRegister qureg)
{
    long long int chunkSize, stateVecSize;
    long long int index;

    // dimension of the state vector
    chunkSize = qureg.numAmpsPerChunk;
    stateVecSize = chunkSize*qureg.numChunks;
    REAL normFactor = 1.0/sqrt((REAL)stateVecSize);

    // Can't use qureg->stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

    // initialise the state to |+++..+++> = 1/normFactor {1, 1, 1, ...}
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkSize, stateVecReal, stateVecImag, normFactor) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<chunkSize; index++) {
            stateVecReal[index] = normFactor;
            stateVecImag[index] = 0.0;
        }
    }
}

void statevec_initClassicalState (QubitRegister qureg, long long int stateInd)
{
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

    // initialise the state to vector to all zeros
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecSize, stateVecReal, stateVecImag) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            stateVecReal[index] = 0.0;
            stateVecImag[index] = 0.0;
        }
    }

	// give the specified classical state prob 1
    if (qureg.chunkId == stateInd/stateVecSize){
        stateVecReal[stateInd % stateVecSize] = 1.0;
        stateVecImag[stateInd % stateVecSize] = 0.0;
    }
}

void statevec_initPureState(QubitRegister targetQureg, QubitRegister copyQureg) {
	
	// registers are equal sized, so nodes hold the same state-vector partitions
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = targetQureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    REAL *targetStateVecReal = targetQureg.stateVec.real;
    REAL *targetStateVecImag = targetQureg.stateVec.imag;
    REAL *copyStateVecReal = copyQureg.stateVec.real;
    REAL *copyStateVecImag = copyQureg.stateVec.imag;

    // initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecSize, targetStateVecReal, targetStateVecImag, copyStateVecReal, copyStateVecImag) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            targetStateVecReal[index] = copyStateVecReal[index];
            targetStateVecImag[index] = copyStateVecImag[index];
        }
    }
}

/**
 * Initialise the state vector of probability amplitudes such that one qubit is set to 'outcome' and all other qubits are in an equal superposition of zero and one.
 * @param[in,out] qureg object representing the set of qubits to be initialised
 * @param[in] qubitId id of qubit to set to state 'outcome'
 * @param[in] value of qubit 'qubitId'
 */
void statevec_initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome)
{
    long long int chunkSize, stateVecSize;
    long long int index;
    int bit;
    const long long int chunkId=qureg->chunkId;

    // dimension of the state vector
    chunkSize = qureg->numAmpsPerChunk;
    stateVecSize = chunkSize*qureg->numChunks;
    REAL normFactor = 1.0/sqrt((REAL)stateVecSize/2.0);

    // Can't use qureg->stateVec as a private OMP var
    REAL *stateVecReal = qureg->stateVec.real;
    REAL *stateVecImag = qureg->stateVec.imag;

    // initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkSize, stateVecReal, stateVecImag, normFactor, qubitId, outcome) \
    private  (index, bit) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<chunkSize; index++) {
            bit = extractBit(qubitId, index+chunkId*chunkSize);
            if (bit==outcome) {
                stateVecReal[index] = normFactor;
                stateVecImag[index] = 0.0;
            } else {
                stateVecReal[index] = 0.0;
                stateVecImag[index] = 0.0;
            }
        }
    }
}


/**
 * Initialise the state vector of probability amplitudes to an (unphysical) state with
 * each component of each probability amplitude a unique floating point value. For debugging processes
 * @param[in,out] qureg object representing the set of qubits to be initialised
 */
void statevec_initStateDebug (QubitRegister qureg)
{
    long long int chunkSize;
    long long int index;
	long long int indexOffset;

    // dimension of the state vector
    chunkSize = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

	indexOffset = chunkSize * qureg.chunkId;

    // initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkSize, stateVecReal, stateVecImag, indexOffset) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<chunkSize; index++) {
            stateVecReal[index] = ((indexOffset + index)*2.0)/10.0;
            stateVecImag[index] = ((indexOffset + index)*2.0+1.0)/10.0;
        }
    }
}

void statevec_initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env){
    long long int chunkSize, stateVecSize;
    long long int indexInChunk, totalIndex;

    chunkSize = qureg->numAmpsPerChunk;
    stateVecSize = chunkSize*qureg->numChunks;

    REAL *stateVecReal = qureg->stateVec.real;
    REAL *stateVecImag = qureg->stateVec.imag;

    FILE *fp;
    char line[200];

    for (int rank=0; rank<(qureg->numChunks); rank++){
        if (rank==qureg->chunkId){
            fp = fopen(filename, "r");
            QuESTAssert(fp!=NULL, 11, __func__);
            indexInChunk = 0; totalIndex = 0;
            while (fgets(line, sizeof(char)*200, fp) != NULL && totalIndex<stateVecSize){
                if (line[0]!='#'){
                    int chunkId = totalIndex/chunkSize;
                    if (chunkId==qureg->chunkId){
                        # if QuEST_PREC==1
                        sscanf(line, "%f, %f", &(stateVecReal[indexInChunk]), 
                                &(stateVecImag[indexInChunk]));	
						# elif QuEST_PREC==2					
                        sscanf(line, "%lf, %lf", &(stateVecReal[indexInChunk]), 
                                &(stateVecImag[indexInChunk]));
						# elif QuEST_PREC==4
		                sscanf(line, "%Lf, %Lf", &(stateVecReal[indexInChunk]), 
		                        &(stateVecImag[indexInChunk]));
						# endif
                        indexInChunk += 1;
                    }
                    totalIndex += 1;
                }
            }	
            fclose(fp);
        }
        syncQuESTEnv(env);
    }
}

int statevec_compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision){
    REAL diff;
    int chunkSize = mq1.numAmpsPerChunk;
	
    for (int i=0; i<chunkSize; i++){
        diff = absReal(mq1.stateVec.real[i] - mq2.stateVec.real[i]);
        if (diff>precision) return 0;
        diff = absReal(mq1.stateVec.imag[i] - mq2.stateVec.imag[i]);
        if (diff>precision) return 0;
    }
    return 1;
}







void statevec_compactUnitaryLocal (QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;
    REAL alphaImag=alpha.imag, alphaReal=alpha.real;
    REAL betaImag=beta.imag, betaReal=beta.real;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            // store current state vector values in temp variables
            stateRealUp = stateVecReal[indexUp];
            stateImagUp = stateVecImag[indexUp];

            stateRealLo = stateVecReal[indexLo];
            stateImagLo = stateVecImag[indexLo];

            // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
            stateVecReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp 
                - betaReal*stateRealLo - betaImag*stateImagLo;
            stateVecImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp 
                - betaReal*stateImagLo + betaImag*stateRealLo;

            // state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
            stateVecReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp 
                + alphaReal*stateRealLo + alphaImag*stateImagLo;
            stateVecImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp 
                + alphaReal*stateImagLo - alphaImag*stateRealLo;
        } 
    }

} 

void statevec_unitaryLocal(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            // store current state vector values in temp variables
            stateRealUp = stateVecReal[indexUp];
            stateImagUp = stateVecImag[indexUp];

            stateRealLo = stateVecReal[indexLo];
            stateImagLo = stateVecImag[indexLo];


            // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
            stateVecReal[indexUp] = u.r0c0.real*stateRealUp - u.r0c0.imag*stateImagUp 
                + u.r0c1.real*stateRealLo - u.r0c1.imag*stateImagLo;
            stateVecImag[indexUp] = u.r0c0.real*stateImagUp + u.r0c0.imag*stateRealUp 
                + u.r0c1.real*stateImagLo + u.r0c1.imag*stateRealLo;

            // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
            stateVecReal[indexLo] = u.r1c0.real*stateRealUp  - u.r1c0.imag*stateImagUp 
                + u.r1c1.real*stateRealLo  -  u.r1c1.imag*stateImagLo;
            stateVecImag[indexLo] = u.r1c0.real*stateImagUp + u.r1c0.imag*stateRealUp 
                + u.r1c1.real*stateImagLo + u.r1c1.imag*stateRealLo;

        } 
    }
} 

/** Rotate a single qubit in the state vector of probability amplitudes, 
 * given two complex numbers alpha and beta, 
 * and a subset of the state vector with upper and lower block values stored seperately.
 *                                                                       
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_compactUnitaryDistributed (QubitRegister qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;

    REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
    REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
    REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // store current state vector values in temp variables
            stateRealUp = stateVecRealUp[thisTask];
            stateImagUp = stateVecImagUp[thisTask];

            stateRealLo = stateVecRealLo[thisTask];
            stateImagLo = stateVecImagLo[thisTask];

            // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
            stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp + rot2Real*stateRealLo + rot2Imag*stateImagLo;
            stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp + rot2Real*stateImagLo - rot2Imag*stateRealLo;
        }
    }
}

/** Apply a unitary operation to a single qubit
 *  given a subset of the state vector with upper and lower block values 
 * stored seperately.
 *
 *  @remarks Qubits are zero-based and the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] u unitary matrix to apply
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_unitaryDistributed (QubitRegister qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;

    REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
    REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
    REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;


# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real, rot1Imag, rot2Real, rot2Imag) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // store current state vector values in temp variables
            stateRealUp = stateVecRealUp[thisTask];
            stateImagUp = stateVecImagUp[thisTask];

            stateRealLo = stateVecRealLo[thisTask];
            stateImagLo = stateVecImagLo[thisTask];

            stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp 
                + rot2Real*stateRealLo - rot2Imag*stateImagLo;
            stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp 
                + rot2Real*stateImagLo + rot2Imag*stateRealLo;
        }
    }
}

void statevec_controlledCompactUnitaryLocal (QubitRegister qureg, const int controlQubit, const int targetQubit, 
        Complex alpha, Complex beta)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;
    REAL alphaImag=alpha.imag, alphaReal=alpha.real;
    REAL betaImag=beta.imag, betaReal=beta.real;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            controlBit = extractBit (controlQubit, indexUp+chunkId*chunkSize);
            if (controlBit){
                // store current state vector values in temp variables
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];

                // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
                stateVecReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp 
                    - betaReal*stateRealLo - betaImag*stateImagLo;
                stateVecImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp 
                    - betaReal*stateImagLo + betaImag*stateRealLo;

                // state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
                stateVecReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp 
                    + alphaReal*stateRealLo + alphaImag*stateImagLo;
                stateVecImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp 
                    + alphaReal*stateImagLo - alphaImag*stateRealLo;
            }
        } 
    }

} 

void statevec_multiControlledUnitaryLocal(QubitRegister qureg, const int targetQubit, 
        long long int mask, ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u, mask) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            if (mask == (mask & (indexUp+chunkId*chunkSize)) ){
                // store current state vector values in temp variables
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];


                // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
                stateVecReal[indexUp] = u.r0c0.real*stateRealUp - u.r0c0.imag*stateImagUp 
                    + u.r0c1.real*stateRealLo - u.r0c1.imag*stateImagLo;
                stateVecImag[indexUp] = u.r0c0.real*stateImagUp + u.r0c0.imag*stateRealUp 
                    + u.r0c1.real*stateImagLo + u.r0c1.imag*stateRealLo;

                // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
                stateVecReal[indexLo] = u.r1c0.real*stateRealUp  - u.r1c0.imag*stateImagUp 
                    + u.r1c1.real*stateRealLo  -  u.r1c1.imag*stateImagLo;
                stateVecImag[indexLo] = u.r1c0.real*stateImagUp + u.r1c0.imag*stateRealUp 
                    + u.r1c1.real*stateImagLo + u.r1c1.imag*stateRealLo;
            }
        } 
    }

}

void statevec_controlledUnitaryLocal(QubitRegister qureg, const int controlQubit, const int targetQubit, 
        ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            controlBit = extractBit (controlQubit, indexUp+chunkId*chunkSize);
            if (controlBit){
                // store current state vector values in temp variables
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];


                // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
                stateVecReal[indexUp] = u.r0c0.real*stateRealUp - u.r0c0.imag*stateImagUp 
                    + u.r0c1.real*stateRealLo - u.r0c1.imag*stateImagLo;
                stateVecImag[indexUp] = u.r0c0.real*stateImagUp + u.r0c0.imag*stateRealUp 
                    + u.r0c1.real*stateImagLo + u.r0c1.imag*stateRealLo;

                // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
                stateVecReal[indexLo] = u.r1c0.real*stateRealUp  - u.r1c0.imag*stateImagUp 
                    + u.r1c1.real*stateRealLo  -  u.r1c1.imag*stateImagLo;
                stateVecImag[indexLo] = u.r1c0.real*stateImagUp + u.r1c0.imag*stateRealUp 
                    + u.r1c1.real*stateImagLo + u.r1c1.imag*stateRealLo;
            }
        } 
    }

}

/** Rotate a single qubit in the state vector of probability amplitudes, given two complex 
 * numbers alpha and beta and a subset of the state vector with upper and lower block values 
 * stored seperately. Only perform the rotation where the control qubit is one.
 *                                               
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_controlledCompactUnitaryDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
    REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
    REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                // store current state vector values in temp variables
                stateRealUp = stateVecRealUp[thisTask];
                stateImagUp = stateVecImagUp[thisTask];

                stateRealLo = stateVecRealLo[thisTask];
                stateImagLo = stateVecImagLo[thisTask];

                // state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
                stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp + rot2Real*stateRealLo + rot2Imag*stateImagLo;
                stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp + rot2Real*stateImagLo - rot2Imag*stateRealLo;
            }
        }
    }
}

/** Rotate a single qubit in the state vector of probability amplitudes, given two complex 
 *  numbers alpha and beta and a subset of the state vector with upper and lower block values 
 *  stored seperately. Only perform the rotation where the control qubit is one.
 *                                                 
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_controlledUnitaryDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
    REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
    REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                // store current state vector values in temp variables
                stateRealUp = stateVecRealUp[thisTask];
                stateImagUp = stateVecImagUp[thisTask];

                stateRealLo = stateVecRealLo[thisTask];
                stateImagLo = stateVecImagLo[thisTask];

                stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp 
                    + rot2Real*stateRealLo - rot2Imag*stateImagLo;
                stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp 
                    + rot2Real*stateImagLo + rot2Imag*stateRealLo;
            }
        }
    }
}

/** Apply a unitary operation to a single qubit in the state vector of probability amplitudes, given
 *  a subset of the state vector with upper and lower block values 
 stored seperately. Only perform the rotation where all the control qubits are 1.
 *                                                 
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_multiControlledUnitaryDistributed (QubitRegister qureg, 
        const int targetQubit, 
        long long int mask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
    REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
    REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag, mask) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            if (mask == (mask & (thisTask+chunkId*chunkSize)) ){
                // store current state vector values in temp variables
                stateRealUp = stateVecRealUp[thisTask];
                stateImagUp = stateVecImagUp[thisTask];

                stateRealLo = stateVecRealLo[thisTask];
                stateImagLo = stateVecImagLo[thisTask];

                stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp 
                    + rot2Real*stateRealLo - rot2Imag*stateImagLo;
                stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp 
                    + rot2Real*stateImagLo + rot2Imag*stateRealLo;
            }
        }
    }
}

void statevec_sigmaXLocal(QubitRegister qureg, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            stateRealUp = stateVecReal[indexUp];
            stateImagUp = stateVecImag[indexUp];

            stateVecReal[indexUp] = stateVecReal[indexLo];
            stateVecImag[indexUp] = stateVecImag[indexLo];

            stateVecReal[indexLo] = stateRealUp;
            stateVecImag[indexLo] = stateImagUp;
        } 
    }

}

/** Rotate a single qubit by {{0,1},{1,0}.
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *  
 *  @remarks Qubits are zero-based and the                     
 *  the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_sigmaXDistributed (QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut)
{

    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;

    REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut) \
    private  (thisTask)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            stateVecRealOut[thisTask] = stateVecRealIn[thisTask];
            stateVecImagOut[thisTask] = stateVecImagIn[thisTask];
        }
    }
} 

void statevec_controlledNotLocal(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,controlBit) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            controlBit = extractBit(controlQubit, indexUp+chunkId*chunkSize);
            if (controlBit){
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                stateVecReal[indexUp] = stateVecReal[indexLo];
                stateVecImag[indexUp] = stateVecImag[indexLo];

                stateVecReal[indexLo] = stateRealUp;
                stateVecImag[indexLo] = stateImagUp;
            }
        } 
    }
}

/** Rotate a single qubit by {{0,1},{1,0}.
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk. Only perform the rotation
 *  for elements where controlQubit is one.
 *                                          
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_controlledNotDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut)
{

    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut) \
    private  (thisTask,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                stateVecRealOut[thisTask] = stateVecRealIn[thisTask];
                stateVecImagOut[thisTask] = stateVecImagIn[thisTask];
            }
        }
    }
} 

void statevec_sigmaYLocal(QubitRegister qureg, const int targetQubit, const int conjFac)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            stateRealUp = stateVecReal[indexUp];
            stateImagUp = stateVecImag[indexUp];

            stateVecReal[indexUp] = conjFac * stateVecImag[indexLo];
            stateVecImag[indexUp] = conjFac * -stateVecReal[indexLo];
            stateVecReal[indexLo] = conjFac * -stateImagUp;
            stateVecImag[indexLo] = conjFac * stateRealUp;
        } 
    }
}

/** Rotate a single qubit by +-{{0,-i},{i,0}.
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *  
 *  @remarks Qubits are zero-based and the                     
 *  the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_sigmaYDistributed(QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, 
        int updateUpper, const int conjFac)
{

    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;

    REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

    int realSign=1, imagSign=1;
    if (updateUpper) imagSign=-1;
    else realSign = -1;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut,realSign,imagSign) \
    private  (thisTask)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            stateVecRealOut[thisTask] = conjFac * realSign * stateVecImagIn[thisTask];
            stateVecImagOut[thisTask] = conjFac * imagSign * stateVecRealIn[thisTask];
        }
    }
} 




void statevec_controlledSigmaYLocal(QubitRegister qureg, const int controlQubit, const int targetQubit, const int conjFac)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,controlBit) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            controlBit = extractBit(controlQubit, indexUp+chunkId*chunkSize);
            if (controlBit){
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

				// update under +-{{0, -i}, {i, 0}}
			    stateVecReal[indexUp] = conjFac * stateVecImag[indexLo];
			    stateVecImag[indexUp] = conjFac * -stateVecReal[indexLo];
			    stateVecReal[indexLo] = conjFac * -stateImagUp;
			    stateVecImag[indexLo] = conjFac * stateRealUp;
            }
        } 
    }
}


void statevec_controlledSigmaYDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, const int conjFac)
{

    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    int controlBit;

    REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut) \
    private  (thisTask,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                stateVecRealOut[thisTask] = conjFac * stateVecRealIn[thisTask];
                stateVecImagOut[thisTask] = conjFac * stateVecImagIn[thisTask];
            }
        }
    }
} 







void statevec_hadamardLocal(QubitRegister qureg, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

    REAL recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, recRoot2) \
    private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;

            stateRealUp = stateVecReal[indexUp];
            stateImagUp = stateVecImag[indexUp];

            stateRealLo = stateVecReal[indexLo];
            stateImagLo = stateVecImag[indexLo];

            stateVecReal[indexUp] = recRoot2*(stateRealUp + stateRealLo);
            stateVecImag[indexUp] = recRoot2*(stateImagUp + stateImagLo);

            stateVecReal[indexLo] = recRoot2*(stateRealUp - stateRealLo);
            stateVecImag[indexLo] = recRoot2*(stateImagUp - stateImagLo);
        } 
    }
}

/** Rotate a single qubit by {{1,1},{1,-1}}/sqrt2.
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *                                          
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_hadamardDistributed(QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut,
        int updateUpper)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=qureg.numAmpsPerChunk;

    int sign;
    if (updateUpper) sign=1;
    else sign=-1;

    REAL recRoot2 = 1.0/sqrt(2);

    REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            recRoot2, sign) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // store current state vector values in temp variables
            stateRealUp = stateVecRealUp[thisTask];
            stateImagUp = stateVecImagUp[thisTask];

            stateRealLo = stateVecRealLo[thisTask];
            stateImagLo = stateVecImagLo[thisTask];

            stateVecRealOut[thisTask] = recRoot2*(stateRealUp + sign*stateRealLo);
            stateVecImagOut[thisTask] = recRoot2*(stateImagUp + sign*stateImagLo);
        }
    }
}

void statevec_phaseShiftByTerm (QubitRegister qureg, const int targetQubit, Complex term)
{
	QuESTAssert(validateUnitComplex(term), 16, __func__);
	
    long long int index;
    long long int stateVecSize;
    int targetBit;
	
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;
	
	REAL stateRealLo, stateImagLo;
	const REAL cosAngle = term.real;
	const REAL sinAngle = term.imag;

# ifdef _OPENMP
# pragma omp parallel for \
    default  (none)			     \
    shared   (stateVecSize, stateVecReal,stateVecImag ) \
    private  (index,targetBit,stateRealLo,stateImagLo)		       \
    schedule (static)
# endif
    for (index=0; index<stateVecSize; index++) {
		
		// update the coeff of the |1> state of the target qubit
        targetBit = extractBit (targetQubit, index+chunkId*chunkSize);
        if (targetBit) {
			
        	stateRealLo = stateVecReal[index];
        	stateImagLo = stateVecImag[index];
			
			stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
			stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;	
        }
    }
}

void statevec_controlledPhaseShift (QubitRegister qureg, const int idQubit1, const int idQubit2, REAL angle)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;
	
    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;
	
	REAL stateRealLo, stateImagLo;
	const REAL cosAngle = cos(angle);
	const REAL sinAngle = sin(angle);

# ifdef _OPENMP
# pragma omp parallel for \
    default  (none)			     \
    shared   (stateVecSize, stateVecReal,stateVecImag ) \
    private  (index,bit1,bit2,stateRealLo,stateImagLo)		       \
    schedule (static)
# endif
    for (index=0; index<stateVecSize; index++) {
        bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
        bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
        if (bit1 && bit2) {
			
        	stateRealLo = stateVecReal[index];
        	stateImagLo = stateVecImag[index];
			
			stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
			stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;	
        }
    }
}

void statevec_multiControlledPhaseShift(QubitRegister qureg, int *controlQubits, int numControlQubits, REAL angle)
{
    long long int index;
    long long int stateVecSize;

    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    QuESTAssert(numControlQubits > 0 && numControlQubits <= qureg.numQubits, 4, __func__);
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) 
		mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<qureg.numQubits)-1, 2, __func__);

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;
	
	REAL stateRealLo, stateImagLo;
	const REAL cosAngle = cos(angle);
	const REAL sinAngle = sin(angle);

# ifdef _OPENMP
# pragma omp parallel \
    default  (none)			     \
    shared   (stateVecSize, stateVecReal, stateVecImag, mask) \
    private  (index, stateRealLo, stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            if (mask == (mask & (index+chunkId*chunkSize)) ){
				
	        	stateRealLo = stateVecReal[index];
	        	stateImagLo = stateVecImag[index];
			
				stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
				stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;	
            }
        }
    }
}

REAL densmatr_findProbabilityOfZeroLocal(QubitRegister qureg, const int measureQubit) {
	
	// computes first local index containing a diagonal element
	long long int localNumAmps = qureg.numAmpsPerChunk;
	long long int densityDim = (1LL << qureg.numDensityQubits);
	long long int diagSpacing = 1LL + densityDim;
	long long int maxNumDiagsPerChunk = localNumAmps / diagSpacing;
	long long int numPrevDiags = (qureg.chunkId * localNumAmps) / diagSpacing;
	long long int globalIndNextDiag = diagSpacing * numPrevDiags;
	long long int localIndNextDiag = globalIndNextDiag % localNumAmps;
	
	// computes how many diagonals are contained in this chunk
	long long int numDiagsInThisChunk = maxNumDiagsPerChunk;
	if (localIndNextDiag + numDiagsInThisChunk*diagSpacing >= localNumAmps)
		numDiagsInThisChunk -= 1;
	
	long long int visitedDiags;		// number of visited diagonals in this chunk so far
	long long int basisStateInd;	// current diagonal index being considered
	long long int index;			// index in the local chunk
	
	REAL zeroProb = 0;
    REAL *stateVecReal = qureg.stateVec.real;
	
# ifdef _OPENMP
# pragma omp parallel \
    shared    (localIndNextDiag, numPrevDiags, diagSpacing, stateVecReal) \
    private   (visitedDiags, basisStateInd, index) \
    reduction ( +:zeroProb )
# endif	
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
		// sums the diagonal elems of the density matrix where measureQubit=0
		for (visitedDiags = 0; visitedDiags < numDiagsInThisChunk; visitedDiags++) {
			
			basisStateInd = numPrevDiags + visitedDiags;
			index = localIndNextDiag + diagSpacing * visitedDiags;

			if (extractBit(measureQubit, basisStateInd++) == 0)
            	zeroProb += stateVecReal[index]; // assume imag[diagonls] ~ 0

		}
    }
    return zeroProb;
}

/** Measure the total probability of a specified qubit being in the zero state across all amplitudes in this chunk.
 *  Size of regions to skip is less than the size of one chunk.                   
 *  
 *  @param[in] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @return probability of qubit measureQubit being zero
 */
REAL statevec_findProbabilityOfZeroLocal (QubitRegister qureg,
        const int measureQubit)
{
    // ----- sizes
    long long int sizeBlock,                                  // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                  // current block
         index;                                               // current index for first half block
    // ----- measured probability
    REAL   totalProbability;                                  // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                         // size of blocks (pairs of measure and skip entries)

    // initialise returned value
    totalProbability = 0.0;

    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
    private   (thisTask,thisBlock,index) \
    reduction ( +:totalProbability )
# endif	
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            thisBlock = thisTask / sizeHalfBlock;
            index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;

            totalProbability += stateVecReal[index]*stateVecReal[index]
                + stateVecImag[index]*stateVecImag[index];
        }
    }
    return totalProbability;
}

/** Measure the probability of a specified qubit being in the zero state across all amplitudes held in this chunk.
 * Size of regions to skip is a multiple of chunkSize.
 * The results are communicated and aggregated by the caller
 *  
 *  @param[in] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @return probability of qubit measureQubit being zero
 */
REAL statevec_findProbabilityOfZeroDistributed (QubitRegister qureg,
        const int measureQubit)
{
    // ----- measured probability
    REAL   totalProbability;                                  // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    long long int numTasks=qureg.numAmpsPerChunk;

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    // initialise returned value
    totalProbability = 0.0;

    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (numTasks,stateVecReal,stateVecImag) \
    private   (thisTask) \
    reduction ( +:totalProbability )
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            totalProbability += stateVecReal[thisTask]*stateVecReal[thisTask]
                + stateVecImag[thisTask]*stateVecImag[thisTask];
        }
    }

    return totalProbability;
}



void statevec_controlledPhaseFlip (QubitRegister qureg, const int idQubit1, const int idQubit2)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;

    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    QuESTAssert(idQubit1 >= 0 && idQubit1 < qureg.numQubits, 2, __func__);
    QuESTAssert(idQubit2 >= 0 && idQubit2 < qureg.numQubits, 1, __func__);
    QuESTAssert(idQubit1 != idQubit2, 3, __func__);

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel for \
    default  (none)			     \
    shared   (stateVecSize, stateVecReal,stateVecImag ) \
    private  (index,bit1,bit2)		       \
    schedule (static)
# endif
    for (index=0; index<stateVecSize; index++) {
        bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
        bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
        if (bit1 && bit2) {
            stateVecReal [index] = - stateVecReal [index];
            stateVecImag [index] = - stateVecImag [index];
        }
    }
}

void statevec_multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits)
{
    long long int index;
    long long int stateVecSize;

    const long long int chunkSize=qureg.numAmpsPerChunk;
    const long long int chunkId=qureg.chunkId;

    QuESTAssert(numControlQubits > 0 && numControlQubits <= qureg.numQubits, 4, __func__);
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<qureg.numQubits)-1, 2, __func__);

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none)			     \
    shared   (stateVecSize, stateVecReal,stateVecImag, mask ) \
    private  (index)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            if (mask == (mask & (index+chunkId*chunkSize)) ){
                stateVecReal [index] = - stateVecReal [index];
                stateVecImag [index] = - stateVecImag [index];
            }
        }
    }
}

/** Update the state vector to be consistent with measuring measureQubit=0 if outcome=0 and measureQubit=1
 *  if outcome=1.
 *  Performs an irreversible change to the state vector: it updates the vector according
 *  to the event that an outcome have been measured on the qubit indicated by measureQubit (where 
 *  this label starts from 0, of course). It achieves this by setting all inconsistent 
 *  amplitudes to 0 and 
 *  then renormalising based on the total probability of measuring measureQubit=0 or 1 according to the 
 *  value of outcome. 
 *  In the local version, one or more blocks (with measureQubit=0 in the first half of the block and
 *  measureQubit=1 in the second half of the block) fit entirely into one chunk. 
 *  
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @param[in] totalProbability probability of qubit measureQubit being either zero or one
 *  @param[in] outcome to measure the probability of and set the state to -- either zero or one
 */
void statevec_collapseToOutcomeLocal(QubitRegister qureg, int measureQubit, REAL totalProbability, int outcome)
{
    // ----- sizes
    long long int sizeBlock,                                  // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                  // current block
         index;                                               // current index for first half block
    // ----- measured probability
    REAL   renorm;                                            // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    // (good for shared memory parallelism)
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                         // size of blocks (pairs of measure and skip entries)

    renorm=1/sqrt(totalProbability);
    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;


# ifdef _OPENMP
# pragma omp parallel \
    default (none) \
    shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag,renorm,outcome) \
    private   (thisTask,thisBlock,index)
# endif
    {
        if (outcome==0){
            // measure qubit is 0
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                thisBlock = thisTask / sizeHalfBlock;
                index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
                stateVecReal[index]=stateVecReal[index]*renorm;
                stateVecImag[index]=stateVecImag[index]*renorm;

                stateVecReal[index+sizeHalfBlock]=0;
                stateVecImag[index+sizeHalfBlock]=0;
            }
        } else {
            // measure qubit is 1
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                thisBlock = thisTask / sizeHalfBlock;
                index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
                stateVecReal[index]=0;
                stateVecImag[index]=0;

                stateVecReal[index+sizeHalfBlock]=stateVecReal[index+sizeHalfBlock]*renorm;
                stateVecImag[index+sizeHalfBlock]=stateVecImag[index+sizeHalfBlock]*renorm;
            }
        }
    }

}

/** Renormalise parts of the state vector where measureQubit=0 or 1, based on the total probability of that qubit being
 *  in state 0 or 1.
 *  Measure in Zero performs an irreversible change to the state vector: it updates the vector according
 *  to the event that the value 'outcome' has been measured on the qubit indicated by measureQubit (where 
 *  this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
 *  then renormalising based on the total probability of measuring measureQubit=0 if outcome=0 and
 *  measureQubit=1 if outcome=1.
 *  In the distributed version, one block (with measureQubit=0 in the first half of the block and
 *  measureQubit=1 in the second half of the block) is spread over multiple chunks, meaning that each chunks performs
 *  only renormalisation or only setting amplitudes to 0. This function handles the renormalisation.
 *  
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @param[in] totalProbability probability of qubit measureQubit being zero
 */
REAL statevec_collapseToOutcomeDistributedRenorm (QubitRegister qureg, const int measureQubit, const REAL totalProbability)
{
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=qureg.numAmpsPerChunk;

    REAL renorm=1/sqrt(totalProbability);

    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (numTasks,stateVecReal,stateVecImag) \
    private   (thisTask)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            stateVecReal[thisTask] = stateVecReal[thisTask]*renorm;
            stateVecImag[thisTask] = stateVecImag[thisTask]*renorm;
        }
    }
    return totalProbability;
}

/** Set all amplitudes in one chunk to 0. 
 *  Measure in Zero performs an irreversible change to the state vector: it updates the vector according
 *  to the event that a zero have been measured on the qubit indicated by measureQubit (where 
 *  this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
 *  then renormalising based on the total probability of measuring measureQubit=0 or 1.
 *  In the distributed version, one block (with measureQubit=0 in the first half of the block and
 *  measureQubit=1 in the second half of the block) is spread over multiple chunks, meaning that each chunks performs
 *  only renormalisation or only setting amplitudes to 0. This function handles setting amplitudes to 0.
 *  
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 */
void statevec_collapseToOutcomeDistributedSetZero(QubitRegister qureg, const int measureQubit)
{
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=qureg.numAmpsPerChunk;

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    REAL *stateVecReal = qureg.stateVec.real;
    REAL *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (numTasks,stateVecReal,stateVecImag) \
    private   (thisTask)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            stateVecReal[thisTask] = 0;
            stateVecImag[thisTask] = 0;
        }
    }
}


