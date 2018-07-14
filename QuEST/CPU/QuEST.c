// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt 
// for details 

/** @file qubits.c
 * The core of the QuEST Library.
 */

# include "../QuEST.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"
# include "QuEST_internal.h"

# define _BSD_SOURCE
# include <unistd.h>
# include <math.h>  
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include <sys/param.h>
# include <sys/types.h> 
# include <sys/time.h>

# ifdef _OPENMP
# include <omp.h>
# endif

const char* errorCodes[] = {
    "Success",                                              // 0
    "Invalid target qubit. Note qubits are zero indexed.",  // 1
    "Invalid control qubit. Note qubits are zero indexed.", // 2 
    "Control qubit cannot equal target qubit.",             // 3
    "Invalid number of control qubits",                     // 4
    "Invalid unitary matrix.",                              // 5
    "Invalid rotation arguments.",                          // 6
    "Invalid system size. Cannot print output for systems greater than 5 qubits.", // 7
    "Can't collapse to state with zero probability.", // 8
    "Invalid number of qubits.", // 9
    "Invalid measurement outcome -- must be either 0 or 1.", // 10
    "Could not open file" // 11
};

static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber);

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env)
{
    QuESTAssert(numQubits>0, 9, __func__);
    long long int numAmps = 1L << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;

    multiQubit->stateVec.real = malloc(numAmpsPerRank * sizeof(*(multiQubit->stateVec.real)));
    multiQubit->stateVec.imag = malloc(numAmpsPerRank * sizeof(*(multiQubit->stateVec.imag)));
    if (env.numRanks>1){
        multiQubit->pairStateVec.real = malloc(numAmpsPerRank * sizeof(*(multiQubit->pairStateVec.real)));
        multiQubit->pairStateVec.imag = malloc(numAmpsPerRank * sizeof(*(multiQubit->pairStateVec.imag)));
    }

    if ( (!(multiQubit->stateVec.real) || !(multiQubit->stateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!");
        exit (EXIT_FAILURE);
    }

    if ( env.numRanks>1 && (!(multiQubit->pairStateVec.real) || !(multiQubit->pairStateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!");
        exit (EXIT_FAILURE);
    }

    multiQubit->numQubits = numQubits;
    multiQubit->numAmpsPerChunk = numAmpsPerRank;
    multiQubit->chunkId = env.rank;
    multiQubit->numChunks = env.numRanks;

}

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env){
    free(multiQubit.stateVec.real);
    free(multiQubit.stateVec.imag);
    if (env.numRanks>1){
        free(multiQubit.pairStateVec.real);
        free(multiQubit.pairStateVec.imag);
    }
}


void reportState(MultiQubit multiQubit){
    FILE *state;
    char filename[100];
    long long int index;
    sprintf(filename, "state_rank_%d.csv", multiQubit.chunkId);
    state = fopen(filename, "w");
    QuESTAssert(state!=NULL, 11, __func__);
    if (multiQubit.chunkId==0) fprintf(state, "real, imag\n");

    for(index=0; index<multiQubit.numAmpsPerChunk; index++){
        fprintf(state, REAL_STRING_FORMAT "," REAL_STRING_FORMAT "\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
    }
    fclose(state);
}

void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank){
    long long int index;
    int rank;
    if (multiQubit.numQubits<=5){
        for (rank=0; rank<multiQubit.numChunks; rank++){
            if (multiQubit.chunkId==rank){
                if (reportRank) {
                    printf("Reporting state from rank %d [\n", multiQubit.chunkId);
                    printf("real, imag\n");
                } else if (rank==0) {
                    printf("Reporting state [\n");
                    printf("real, imag\n");
                }

                for(index=0; index<multiQubit.numAmpsPerChunk; index++){
                    printf(REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
                }
                if (reportRank || rank==multiQubit.numChunks-1) printf("]\n");
            }
            syncQuESTEnv(env);
        }
    } else printf("Error: reportStateToScreen will not print output for systems of more than 5 qubits.\n");
}

void reportMultiQubitParams(MultiQubit multiQubit){
    long long int numAmps = 1L << multiQubit.numQubits;
    long long int numAmpsPerRank = numAmps/multiQubit.numChunks;
    if (multiQubit.chunkId==0){
        printf("QUBITS:\n");
        printf("Number of qubits is %d.\n", multiQubit.numQubits);
        printf("Number of amps is %lld.\n", numAmps);
        printf("Number of amps per rank is %lld.\n", numAmpsPerRank);
    }
}

int getNumQubits(MultiQubit multiQubit){
    return multiQubit.numQubits;
}

int getNumAmps(MultiQubit multiQubit){
    return multiQubit.numAmpsPerChunk*multiQubit.numChunks;
}


void getEnvironmentString(QuESTEnv env, MultiQubit multiQubit, char str[200]){
    int numThreads=1;
# ifdef _OPENMP
    numThreads=omp_get_max_threads(); 
# endif
    sprintf(str, "%dqubits_CPU_%dranksx%dthreads", multiQubit.numQubits, env.numRanks, numThreads);
}

void initStateZero (MultiQubit multiQubit)
{
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = multiQubit.numAmpsPerChunk;

    // Can't use multiQubit->stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

    // initialise the state to |0000..0000>
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

    if (multiQubit.chunkId==0){
        // zero state |0000..0000> has probability 1
        stateVecReal[0] = 1.0;
        stateVecImag[0] = 0.0;
    }
}

void initStatePlus (MultiQubit multiQubit)
{
    long long int chunkSize, stateVecSize;
    long long int index;

    // dimension of the state vector
    chunkSize = multiQubit.numAmpsPerChunk;
    stateVecSize = chunkSize*multiQubit.numChunks;
    REAL normFactor = 1.0/sqrt((REAL)stateVecSize);

    // Can't use multiQubit->stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

    // initialise the state to |0000..0000>
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

/* Tyson Jones, 16th May 2018 4pm */
void initClassicalState (MultiQubit multiQubit, long long int stateInd)
{
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = multiQubit.numAmpsPerChunk;

    // Can't use multiQubit->stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

    // initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateInd, stateVecSize, stateVecReal, stateVecImag) \
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
    if (multiQubit.chunkId == stateInd/stateVecSize){
        stateVecReal[stateInd % stateVecSize] = 1.0;
        stateVecImag[stateInd % stateVecSize] = 0.0;
    }
}

/**
 * Initialise the state vector of probability amplitudes such that one qubit is set to 'outcome' and all other qubits are in an equal superposition of zero and one.
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 * @param[in] qubitId id of qubit to set to state 'outcome'
 * @param[in] value of qubit 'qubitId'
 */
void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome)
{
    long long int chunkSize, stateVecSize;
    long long int index;
    int bit;
    const long long int chunkId=multiQubit->chunkId;

    // dimension of the state vector
    chunkSize = multiQubit->numAmpsPerChunk;
    stateVecSize = chunkSize*multiQubit->numChunks;
    REAL normFactor = 1.0/sqrt((REAL)stateVecSize/2.0);

    // Can't use multiQubit->stateVec as a private OMP var
    REAL *stateVecReal = multiQubit->stateVec.real;
    REAL *stateVecImag = multiQubit->stateVec.imag;

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
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 */
void initStateDebug (MultiQubit multiQubit)
{
    long long int chunkSize;
    long long int index;
	long long int indexOffset;

    // dimension of the state vector
    chunkSize = multiQubit.numAmpsPerChunk;

    // Can't use multiQubit->stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

	indexOffset = chunkSize * multiQubit.chunkId;

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

void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env){
    long long int chunkSize, stateVecSize;
    long long int indexInChunk, totalIndex;

    chunkSize = multiQubit->numAmpsPerChunk;
    stateVecSize = chunkSize*multiQubit->numChunks;

    REAL *stateVecReal = multiQubit->stateVec.real;
    REAL *stateVecImag = multiQubit->stateVec.imag;

    FILE *fp;
    char line[200];

    for (int rank=0; rank<(multiQubit->numChunks); rank++){
        if (rank==multiQubit->chunkId){
            fp = fopen(filename, "r");
            QuESTAssert(fp!=NULL, 11, __func__);
            indexInChunk = 0; totalIndex = 0;
            while (fgets(line, sizeof(char)*200, fp) != NULL && totalIndex<stateVecSize){
                if (line[0]!='#'){
                    int chunkId = totalIndex/chunkSize;
                    if (chunkId==multiQubit->chunkId){
                        //! fix -- format needs to work for single precision values
                        sscanf(line, "%lf, %lf", &(stateVecReal[indexInChunk]), 
                                &(stateVecImag[indexInChunk]));
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

int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision){
    REAL diff;
    int chunkSize = mq1.numAmpsPerChunk;
    for (int i=0; i<chunkSize; i++){
        diff = fabs(mq1.stateVec.real[i] - mq2.stateVec.real[i]);
        if (diff>precision) return 0;
        diff = fabs(mq1.stateVec.imag[i] - mq2.stateVec.imag[i]);
        if (diff>precision) return 0;
    }
    return 1;
}

int validateMatrixIsUnitary(ComplexMatrix2 u){

    if ( fabs(u.r0c0.real*u.r0c0.real 
                + u.r0c0.imag*u.r0c0.imag
                + u.r1c0.real*u.r1c0.real
                + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    if ( fabs(u.r0c1.real*u.r0c1.real 
                + u.r0c1.imag*u.r0c1.imag
                + u.r1c1.real*u.r1c1.real
                + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;

    if ( fabs(u.r0c0.real*u.r0c1.real 
                + u.r0c0.imag*u.r0c1.imag
                + u.r1c0.real*u.r1c1.real
                + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;

    if ( fabs(u.r0c1.real*u.r0c0.imag
                - u.r0c0.real*u.r0c1.imag
                + u.r1c1.real*u.r1c0.imag
                - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;

    return 1;
}

int validateAlphaBeta(Complex alpha, Complex beta){
    if ( fabs(alpha.real*alpha.real 
                + alpha.imag*alpha.imag
                + beta.real*beta.real 
                + beta.imag*beta.imag - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateUnitVector(REAL ux, REAL uy, REAL uz){
    if ( fabs(sqrt(ux*ux + uy*uy + uz*uz) - 1) > REAL_EPS ) return 0;
    else return 1;
}

void rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector axis){

    double mag = sqrt(pow(axis.x,2) + pow(axis.y,2) + pow(axis.z,2));
    Vector unitAxis = {axis.x/mag, axis.y/mag, axis.z/mag};

    Complex alpha, beta;
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;	
    beta.real = sin(angle/2.0)*unitAxis.y;
    beta.imag = -sin(angle/2.0)*unitAxis.x;
    compactUnitary(multiQubit, rotQubit, alpha, beta);
}

void rotateX(MultiQubit multiQubit, const int rotQubit, REAL angle){

    Vector unitAxis = {1, 0, 0};
    rotateAroundAxis(multiQubit, rotQubit, angle, unitAxis);
}

void rotateY(MultiQubit multiQubit, const int rotQubit, REAL angle){

    Vector unitAxis = {0, 1, 0};
    rotateAroundAxis(multiQubit, rotQubit, angle, unitAxis);
}

void rotateZ(MultiQubit multiQubit, const int rotQubit, REAL angle){

    Vector unitAxis = {0, 0, 1};
    rotateAroundAxis(multiQubit, rotQubit, angle, unitAxis);
}

void controlledRotateAroundAxis(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle, Vector axis){

    double mag = sqrt(pow(axis.x,2) + pow(axis.y,2) + pow(axis.z,2));
    Vector unitAxis = {axis.x/mag, axis.y/mag, axis.z/mag};

    Complex alpha, beta;
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;	
    beta.real = sin(angle/2.0)*unitAxis.y;
    beta.imag = -sin(angle/2.0)*unitAxis.x;
    controlledCompactUnitary(multiQubit, controlQubit, targetQubit, alpha, beta);
}

void controlledRotateX(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {1, 0, 0};
    controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, unitAxis);
}

void controlledRotateY(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {0, 1, 0};
    controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, unitAxis);
}

void controlledRotateZ(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {0, 0, 1};
    controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, unitAxis);
}

void compactUnitaryLocal (MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;
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

void unitaryLocal(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void compactUnitaryDistributed (MultiQubit multiQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] u unitary matrix to apply
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void unitaryDistributed (MultiQubit multiQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;

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

void controlledCompactUnitaryLocal (MultiQubit multiQubit, const int controlQubit, const int targetQubit, 
        Complex alpha, Complex beta)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;
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

void multiControlledUnitaryLocal(MultiQubit multiQubit, const int targetQubit, 
        long long int mask, ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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

void controlledUnitaryLocal(MultiQubit multiQubit, const int controlQubit, const int targetQubit, 
        ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void controlledCompactUnitaryDistributed (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void controlledUnitaryDistributed (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void multiControlledUnitaryDistributed (MultiQubit multiQubit, 
        const int targetQubit, 
        long long int mask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

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

void sigmaXLocal(MultiQubit multiQubit, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void sigmaXDistributed (MultiQubit multiQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut)
{

    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;

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

void controlledNotLocal(MultiQubit multiQubit, const int controlQubit, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 


    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void controlledNotDistributed (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut)
{

    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;
    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

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

void sigmaYLocal(MultiQubit multiQubit, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateImagUp;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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

            stateVecReal[indexUp] = stateVecImag[indexLo];
            stateVecImag[indexUp] = -stateVecReal[indexLo];

            stateVecReal[indexLo] = -stateImagUp;
            stateVecImag[indexLo] = stateRealUp;
        } 
    }
}

/** Rotate a single qubit by {{0,-i},{i,0}.
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *  
 *  @remarks Qubits are zero-based and the                     
 *  the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void sigmaYDistributed(MultiQubit multiQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, 
        int updateUpper)
{

    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;

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
            stateVecRealOut[thisTask] = realSign*stateVecImagIn[thisTask];
            stateVecImagOut[thisTask] = imagSign*stateVecRealIn[thisTask];
        }
    }
} 

void hadamardLocal(MultiQubit multiQubit, const int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void hadamardDistributed(MultiQubit multiQubit, const int targetQubit,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut,
        int updateUpper)
{

    REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    const long long int numTasks=multiQubit.numAmpsPerChunk;

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

void phaseGateLocal(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    REAL stateRealLo,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

    REAL recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock,stateVecReal,stateVecImag,recRoot2,type) \
    private  (thisTask,thisBlock,indexUp,indexLo,stateRealLo,stateImagLo) 
# endif
    {
        if (type==SIGMA_Z){
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                //! fix -- can i rewrite this to not use mod?
                thisBlock   = thisTask / sizeHalfBlock;
                indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
                indexLo     = indexUp + sizeHalfBlock;

                stateVecReal[indexLo] = -stateVecReal[indexLo];
                stateVecImag[indexLo] = -stateVecImag[indexLo];
            } 
        } 

        else if (type==S_GATE){
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                //! fix -- can i rewrite this to not use mod?
                thisBlock   = thisTask / sizeHalfBlock;
                indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
                indexLo     = indexUp + sizeHalfBlock;
                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];

                stateVecReal[indexLo] = -stateImagLo;
                stateVecImag[indexLo] = stateRealLo;
            } 
        } else if (type==T_GATE){
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                //! fix -- can i rewrite this to not use mod?
                thisBlock   = thisTask / sizeHalfBlock;
                indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
                indexLo     = indexUp + sizeHalfBlock;
                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];

                stateVecReal[indexLo] = recRoot2 * (stateRealLo - stateImagLo);
                stateVecImag[indexLo] = recRoot2 * (stateRealLo + stateImagLo);
            } 
        } else printf("Type %d is an invalid phase gate\n", type);
    }
}

void phaseGateDistributed(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type)
{
    REAL stateRealLo,stateImagLo;
    long long int thisTask;         
    const long long int numTasks=multiQubit.numAmpsPerChunk;

    // Can't use multiQubit.stateVec as a private OMP var
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

    REAL recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecReal,stateVecImag, recRoot2, type) \
    private  (thisTask,stateRealLo,stateImagLo) 
# endif
    {
        if (type==SIGMA_Z){
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                stateVecReal[thisTask] = -stateVecReal[thisTask];
                stateVecImag[thisTask] = -stateVecImag[thisTask];
            } 
        } else if (type==S_GATE){
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                stateRealLo = stateVecReal[thisTask];
                stateImagLo = stateVecImag[thisTask];

                stateVecReal[thisTask] = -stateImagLo;
                stateVecImag[thisTask] = stateRealLo;
            } 
        } else if (type==T_GATE){
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
            for (thisTask=0; thisTask<numTasks; thisTask++) {
                stateRealLo = stateVecReal[thisTask];
                stateImagLo = stateVecImag[thisTask];

                stateVecReal[thisTask] = recRoot2 * (stateRealLo - stateImagLo);
                stateVecImag[thisTask] = recRoot2 * (stateRealLo + stateImagLo);
            } 
        } else printf("Type %d is an invalid phase gate\n", type);
    }
}

void sigmaZ(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, SIGMA_Z);
}

void sGate(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, S_GATE);
} 

void tGate(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, T_GATE);
}

/** Measure the total probability of a specified qubit being in the zero state across all amplitudes in this chunk.
 *  Size of regions to skip is less than the size of one chunk.                   
 *  
 *  @param[in] multiQubit object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @return probability of qubit measureQubit being zero
 */
REAL findProbabilityOfZeroLocal (MultiQubit multiQubit,
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
    long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                         // size of blocks (pairs of measure and skip entries)

    // initialise returned value
    totalProbability = 0.0;

    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  Size of regions to skip is a multiple of chunkSize.
 *  
 *  @param[in] multiQubit object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @return probability of qubit measureQubit being zero
 */
REAL findProbabilityOfZeroDistributed (MultiQubit multiQubit,
        const int measureQubit)
{
    // ----- measured probability
    REAL   totalProbability;                                  // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    long long int numTasks=multiQubit.numAmpsPerChunk;

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    // initialise returned value
    totalProbability = 0.0;

    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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

void controlledPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;

    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

    QuESTAssert(idQubit1 >= 0 && idQubit1 < multiQubit.numQubits, 2, __func__);
    QuESTAssert(idQubit2 >= 0 && idQubit2 < multiQubit.numQubits, 1, __func__);
    QuESTAssert(idQubit1 != idQubit2, 3, __func__);

    // dimension of the state vector
    stateVecSize = multiQubit.numAmpsPerChunk;
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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

void multiControlledPhaseGate(MultiQubit multiQubit, int *controlQubits, int numControlQubits)
{
    long long int index;
    long long int stateVecSize;

    const long long int chunkSize=multiQubit.numAmpsPerChunk;
    const long long int chunkId=multiQubit.chunkId;

    QuESTAssert(numControlQubits > 0 && numControlQubits <= multiQubit.numQubits, 4, __func__);
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<multiQubit.numQubits)-1, 2, __func__);

    stateVecSize = multiQubit.numAmpsPerChunk;
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @param[in] totalProbability probability of qubit measureQubit being either zero or one
 *  @param[in] outcome to measure the probability of and set the state to -- either zero or one
 */
void collapseToOutcomeLocal(MultiQubit multiQubit, int measureQubit, REAL totalProbability, int outcome)
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
    long long int numTasks=multiQubit.numAmpsPerChunk>>1;

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                         // size of blocks (pairs of measure and skip entries)

    renorm=1/sqrt(totalProbability);
    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;


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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @param[in] totalProbability probability of qubit measureQubit being zero
 */
REAL collapseToOutcomeDistributedRenorm (MultiQubit multiQubit, const int measureQubit, const REAL totalProbability)
{
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=multiQubit.numAmpsPerChunk;

    REAL renorm=1/sqrt(totalProbability);

    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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
 *  @param[in,out] multiQubit object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 */
void collapseToOutcomeDistributedSetZero(MultiQubit multiQubit, const int measureQubit)
{
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=multiQubit.numAmpsPerChunk;

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    REAL *stateVecReal = multiQubit.stateVec.real;
    REAL *stateVecImag = multiQubit.stateVec.imag;

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

REAL getProbEl(MultiQubit multiQubit, long long int index){
    REAL real;
    REAL imag;
    real = getRealAmpEl(multiQubit, index);
    imag = getImagAmpEl(multiQubit, index);
    return real*real + imag*imag;
}

void seedQuESTDefault(){
    // init MT random number generator with three keys -- time, pid and a hash of hostname 
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process

    struct timeval  tv;
    gettimeofday(&tv, NULL);

    double time_in_mill = 
        (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ; // convert tv_sec & tv_usec to millisecond

    unsigned long int pid = getpid();
    unsigned long int msecs = (unsigned long int) time_in_mill;
    char hostName[MAXHOSTNAMELEN+1];
    gethostname(hostName, sizeof(hostName));
    unsigned long int hostNameInt = hashString(hostName);

    unsigned long int key[3];
    key[0] = msecs; key[1] = pid; key[2] = hostNameInt;
    init_by_array(key, 3); 
}

/** 
 * numSeeds <= 64
 */
void seedQuEST(unsigned long int *seedArray, int numSeeds){
    // init MT random number generator with user defined list of seeds
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process
    init_by_array(seedArray, numSeeds); 
}

unsigned long int hashString(char *str){
    unsigned long int hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;    
}


