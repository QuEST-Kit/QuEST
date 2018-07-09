// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
  An implementation of the API in qubits.h for a local (non-MPI) environment.
 */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"

# include "mt19937ar.h" // MT random number generation

# define REDUCE_SHARED_SIZE 512
# define DEBUG 0

static __device__ int extractBit (int locationOfBitFromRight, long long int theEncodedNumber)
{
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

#ifdef __cplusplus
extern "C" {
#endif

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env)
{
    QuESTAssert(numQubits>0, 9, __func__);
    // Allocate CPU memory
    long long int numAmps = 1L << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;

    multiQubit->stateVec.real = (REAL*) malloc(numAmpsPerRank * sizeof(multiQubit->stateVec.real));
    multiQubit->stateVec.imag = (REAL*) malloc(numAmpsPerRank * sizeof(multiQubit->stateVec.imag));
    if (env.numRanks>1){
        multiQubit->pairStateVec.real = (REAL*) malloc(numAmpsPerRank * sizeof(multiQubit->pairStateVec.real));
        multiQubit->pairStateVec.imag = (REAL*) malloc(numAmpsPerRank * sizeof(multiQubit->pairStateVec.imag));
    }

    if ( (!(multiQubit->stateVec.real) || !(multiQubit->stateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!\n");
        exit (EXIT_FAILURE);
    }

    if ( env.numRanks>1 && (!(multiQubit->pairStateVec.real) || !(multiQubit->pairStateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!\n");
        exit (EXIT_FAILURE);
    }

    multiQubit->numQubits = numQubits;
    multiQubit->numAmpsDividedByNumChunks = numAmpsPerRank;
    multiQubit->chunkId = env.rank;
    multiQubit->numChunks = env.numRanks;

    // Allocate GPU memory
    cudaMalloc(&(multiQubit->deviceStateVec.real), multiQubit->numAmpsDividedByNumChunks*sizeof(*(multiQubit->deviceStateVec.real)));
    cudaMalloc(&(multiQubit->deviceStateVec.imag), multiQubit->numAmpsDividedByNumChunks*sizeof(*(multiQubit->deviceStateVec.imag)));
    cudaMalloc(&(multiQubit->firstLevelReduction), ceil(multiQubit->numAmpsDividedByNumChunks/(REAL)REDUCE_SHARED_SIZE)*sizeof(REAL));
    cudaMalloc(&(multiQubit->secondLevelReduction), ceil(multiQubit->numAmpsDividedByNumChunks/(REAL)(REDUCE_SHARED_SIZE*REDUCE_SHARED_SIZE))*
            sizeof(REAL));

    if (!(multiQubit->deviceStateVec.real) || !(multiQubit->deviceStateVec.imag)){
        printf("Could not allocate memory on GPU!\n");
        exit (EXIT_FAILURE);
    }

}

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env)
{
    // Free CPU memory
    free(multiQubit.stateVec.real);
    free(multiQubit.stateVec.imag);
    if (env.numRanks>1){
        free(multiQubit.pairStateVec.real);
        free(multiQubit.pairStateVec.imag);
    }

    // Free GPU memory
    cudaFree(multiQubit.deviceStateVec.real);
    cudaFree(multiQubit.deviceStateVec.imag);
}

int GPUExists(void){
    int deviceCount, device;
    int gpuDeviceCount = 0;
    struct cudaDeviceProp properties;
    cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);
    if (cudaResultCode != cudaSuccess) deviceCount = 0;
    /* machines with no GPUs can still report one emulation device */
    for (device = 0; device < deviceCount; ++device) {
        cudaGetDeviceProperties(&properties, device);
        if (properties.major != 9999) { /* 9999 means emulation only */
            ++gpuDeviceCount;
        }
    }
    if (gpuDeviceCount) return 1;
    else return 0;
}

void initQuESTEnv(QuESTEnv *env){
    // init MPI environment
    if (!GPUExists()){
        printf("Trying to run GPU code with no GPU available\n");
        exit(EXIT_FAILURE);
    }
    env->rank=0;
    env->numRanks=1;
}

void syncQuESTEnv(QuESTEnv env){
    cudaDeviceSynchronize();
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void closeQuESTEnv(QuESTEnv env){
    // MPI finalize goes here in MPI version. Call this function anyway for consistency
}

void reportQuESTEnv(QuESTEnv env){
    printf("EXECUTION ENVIRONMENT:\n");
    printf("Running locally on one node with GPU\n");
    printf("Number of ranks is %d\n", env.numRanks);
# ifdef _OPENMP
    printf("OpenMP enabled\n");
    printf("Number of threads available is %d\n", omp_get_max_threads());
# else
    printf("OpenMP disabled\n");
# endif
}

void getEnvironmentString(QuESTEnv env, MultiQubit multiQubit, char str[200]){
    sprintf(str, "%dqubits_GPU_noMpi_noOMP", multiQubit.numQubits);	
}

void copyStateToGPU(MultiQubit multiQubit)
{
    if (DEBUG) printf("Copying data to GPU\n");
    cudaMemcpy(multiQubit.deviceStateVec.real, multiQubit.stateVec.real, 
            multiQubit.numAmpsDividedByNumChunks*sizeof(*(multiQubit.deviceStateVec.real)), cudaMemcpyHostToDevice);
    cudaMemcpy(multiQubit.deviceStateVec.imag, multiQubit.stateVec.imag, 
            multiQubit.numAmpsDividedByNumChunks*sizeof(*(multiQubit.deviceStateVec.imag)), cudaMemcpyHostToDevice);
    if (DEBUG) printf("Finished copying data to GPU\n");
}

void copyStateFromGPU(MultiQubit multiQubit)
{
    cudaDeviceSynchronize();
    if (DEBUG) printf("Copying data from GPU\n");
    cudaMemcpy(multiQubit.stateVec.real, multiQubit.deviceStateVec.real, 
            multiQubit.numAmpsDividedByNumChunks*sizeof(*(multiQubit.deviceStateVec.real)), cudaMemcpyDeviceToHost);
    cudaMemcpy(multiQubit.stateVec.imag, multiQubit.deviceStateVec.imag, 
            multiQubit.numAmpsDividedByNumChunks*sizeof(*(multiQubit.deviceStateVec.imag)), cudaMemcpyDeviceToHost);
    if (DEBUG) printf("Finished copying data from GPU\n");
}

/** Print the current state vector of probability amplitudes for a set of qubits to standard out. 
  For debugging purposes. Each rank should print output serially. Only print output for systems <= 5 qubits
 */
void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank){
    long long int index;
    int rank;
    copyStateFromGPU(multiQubit); 
    if (multiQubit.numQubits<=5){
        for (rank=0; rank<multiQubit.numChunks; rank++){
            if (multiQubit.chunkId==rank){
                if (reportRank) {
                    printf("Reporting state from rank %d [\n", multiQubit.chunkId);
                    //printf("\trank, index, real, imag\n");
                    printf("real, imag\n");
                } else if (rank==0) {
                    printf("Reporting state [\n");
                    printf("real, imag\n");
                }

                for(index=0; index<multiQubit.numAmpsDividedByNumChunks; index++){
                    printf(REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
                }
                if (reportRank || rank==multiQubit.numChunks-1) printf("]\n");
            }
            syncQuESTEnv(env);
        }
    }
}

void __global__ initStateZeroKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
    long long int index;

    // initialise the state to |0000..0000>
    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;
    stateVecReal[index] = 0.0;
    stateVecImag[index] = 0.0;

    if (index==0){
        // zero state |0000..0000> has probability 1
        stateVecReal[0] = 1.0;
        stateVecImag[0] = 0.0;
    }
}

void initStateZero(MultiQubit *multiQubit)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit->numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    initStateZeroKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit->numAmpsDividedByNumChunks, multiQubit->deviceStateVec.real, 
            multiQubit->deviceStateVec.imag);
}

void __global__ initStatePlusKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    REAL normFactor = 1.0/sqrt((REAL)stateVecSize);
    stateVecReal[index] = normFactor;
    stateVecImag[index] = 0.0;
}

void initStatePlus(MultiQubit *multiQubit)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit->numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    initStatePlusKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit->numAmpsDividedByNumChunks, multiQubit->deviceStateVec.real, 
            multiQubit->deviceStateVec.imag);
}

/* Tyson Jones, 4th July 2018 8pm */
void __global__ initClassicalStateKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag, long long int stateInd){
    long long int index;

    // initialise the state to |stateInd>
    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;
    stateVecReal[index] = 0.0;
    stateVecImag[index] = 0.0;

    if (index==stateInd){
        // classical state has probability 1
        stateVecReal[stateInd] = 1.0;
        stateVecImag[stateInd] = 0.0;
    }
}
void initClassicalState(MultiQubit *multiQubit, long long int stateInd)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit->numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    initClassicalStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit->numAmpsDividedByNumChunks, multiQubit->deviceStateVec.real, 
            multiQubit->deviceStateVec.imag, stateInd);
}

void __global__ initStateDebugKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    stateVecReal[index] = (index*2.0)/10.0;
    stateVecImag[index] = (index*2.0+1.0)/10.0;
}

void initStateDebug(MultiQubit *multiQubit)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit->numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    initStateDebugKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit->numAmpsDividedByNumChunks, multiQubit->deviceStateVec.real, 
            multiQubit->deviceStateVec.imag);
}

void __global__ initStateOfSingleQubitKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag, int qubitId, int outcome){
    long long int index;
    int bit;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    REAL normFactor = 1.0/sqrt((REAL)stateVecSize/2);
    bit = extractBit(qubitId, index);
    if (bit==outcome) {
        stateVecReal[index] = normFactor;
        stateVecImag[index] = 0.0;
    } else {
        stateVecReal[index] = 0.0;
        stateVecImag[index] = 0.0;
    }
}

void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit->numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    initStateOfSingleQubitKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit->numAmpsDividedByNumChunks, multiQubit->deviceStateVec.real, multiQubit->deviceStateVec.imag, qubitId, outcome);
}

void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env){
    long long int chunkSize, stateVecSize;
    long long int indexInChunk, totalIndex;

    chunkSize = multiQubit->numAmpsDividedByNumChunks;
    stateVecSize = chunkSize*multiQubit->numChunks;

    REAL *stateVecReal = multiQubit->stateVec.real;
    REAL *stateVecImag = multiQubit->stateVec.imag;

    FILE *fp;
    char line[200];

    fp = fopen(filename, "r");
    indexInChunk = 0; totalIndex = 0;
    while (fgets(line, sizeof(char)*200, fp) != NULL && totalIndex<stateVecSize){
        if (line[0]!='#'){
            int chunkId = totalIndex/chunkSize;
            if (chunkId==multiQubit->chunkId){
                //! fix -- hacky
                if (QuEST_PREC==1){
                    sscanf(line, "%f, %f", &(stateVecReal[indexInChunk]),
                            &(stateVecImag[indexInChunk]));
                } else {
                    sscanf(line, "%lf, %lf", &(stateVecReal[indexInChunk]),
                            &(stateVecImag[indexInChunk]));
                }
                indexInChunk += 1;
            }
            totalIndex += 1;
        }
    }
    fclose(fp);

    copyStateToGPU(*multiQubit);
}

int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision){
    REAL diff;
    int chunkSize = mq1.numAmpsDividedByNumChunks;

    copyStateFromGPU(mq1);
    copyStateFromGPU(mq2);

    for (int i=0; i<chunkSize; i++){
        diff = mq1.stateVec.real[i] - mq2.stateVec.real[i];
        if (diff<0) diff *= -1;
        if (diff>precision) return 0;
        diff = mq1.stateVec.imag[i] - mq2.stateVec.imag[i];
        if (diff<0) diff *= -1;
        if (diff>precision) return 0;
    }
    return 1;
}


REAL calcTotalProbability(MultiQubit multiQubit){
    /* IJB - implemented using Kahan summation for greater accuracy at a slight floating
       point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    /* Don't change the bracketing in this routine! */
    REAL pTotal=0;
    REAL y, t, c;
    long long int index;
    long long int numAmpsPerRank = multiQubit.numAmpsDividedByNumChunks;

    copyStateFromGPU(multiQubit);

    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){
        /* Perform pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index]; by Kahan */
        // pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];

        y = multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index] - c;
        t = pTotal + y;
        c = ( t - pTotal ) - y;
        pTotal = t;

        /* Perform pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index]; by Kahan */
        //pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];


        y = multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index] - c;
        t = pTotal + y;
        c = ( t - pTotal ) - y;
        pTotal = t;


    }
    return pTotal;
}


__global__ void compactUnitaryKernel (MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    sizeHalfBlock = 1LL << rotQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;
    REAL alphaImag=alpha.imag, alphaReal=alpha.real;
    REAL betaImag=beta.imag, betaReal=beta.real;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

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

void compactUnitary(MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    compactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, targetQubit, alpha, beta);
}

__global__ void controlledCompactUnitaryKernel (MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;
    int controlBit;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;
    REAL alphaImag=alpha.imag, alphaReal=alpha.real;
    REAL betaImag=beta.imag, betaReal=beta.real;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    controlBit = extractBit(controlQubit, indexUp);
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

void controlledCompactUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    controlledCompactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, controlQubit, targetQubit, alpha, beta);
}

__global__ void unitaryKernel(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

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

void unitary(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    unitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, targetQubit, u);
}

__global__ void controlledUnitaryKernel(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    int controlBit;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    // store current state vector values in temp variables
    stateRealUp = stateVecReal[indexUp];
    stateImagUp = stateVecImag[indexUp];

    stateRealLo = stateVecReal[indexLo];
    stateImagLo = stateVecImag[indexLo];

    controlBit = extractBit(controlQubit, indexUp);
    if (controlBit){
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

void controlledUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    controlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, controlQubit, targetQubit, u);
}

__global__ void multiControlledUnitaryKernel(MultiQubit multiQubit, long long int mask, const int targetQubit, ComplexMatrix2 u){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;


    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    if (mask == (mask & indexUp) ){
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

void multiControlledUnitary(MultiQubit multiQubit, int *controlQubits, int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(numControlQubits > 0 && numControlQubits <= multiQubit.numQubits, 4, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<multiQubit.numQubits)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<targetQubit)) != (1LL<<targetQubit), 3, __func__);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);

    multiControlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, mask, targetQubit, u);
}

__global__ void sigmaXKernel(MultiQubit multiQubit, const int targetQubit){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,                             // storage for previous state values
           stateImagUp;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    // store current state vector values in temp variables
    stateRealUp = stateVecReal[indexUp];
    stateImagUp = stateVecImag[indexUp];

    stateVecReal[indexUp] = stateVecReal[indexLo];
    stateVecImag[indexUp] = stateVecImag[indexLo];

    stateVecReal[indexLo] = stateRealUp;
    stateVecImag[indexLo] = stateImagUp;
}

void sigmaX(MultiQubit multiQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    sigmaXKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, targetQubit);
}


__global__ void sigmaYKernel(MultiQubit multiQubit, const int targetQubit){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,                             // storage for previous state values
           stateImagUp;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    // store current state vector values in temp variables
    stateRealUp = stateVecReal[indexUp];
    stateImagUp = stateVecImag[indexUp];

    stateVecReal[indexUp] = stateVecImag[indexLo];
    stateVecImag[indexUp] = -stateVecReal[indexLo];

    stateVecReal[indexLo] = -stateImagUp;
    stateVecImag[indexLo] = stateRealUp;
}

void sigmaY(MultiQubit multiQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    sigmaYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, targetQubit);
}

__global__ void phaseGateKernel(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealLo,                             // storage for previous state values
           stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    REAL recRoot2 = 1.0/sqrt(2.0);

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    if (type==SIGMA_Z){
        stateVecReal[indexLo] = -stateVecReal[indexLo];
        stateVecImag[indexLo] = -stateVecImag[indexLo];
    } else if (type==S_GATE){
        stateRealLo = stateVecReal[indexLo];
        stateImagLo = stateVecImag[indexLo];

        stateVecReal[indexLo] = -stateImagLo;
        stateVecImag[indexLo] = stateRealLo;
    } else if (type==T_GATE){
        stateRealLo = stateVecReal[indexLo];
        stateImagLo = stateVecImag[indexLo];

        stateVecReal[indexLo] = recRoot2 * (stateRealLo - stateImagLo);
        stateVecImag[indexLo] = recRoot2 * (stateRealLo + stateImagLo);
    }

}

void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    phaseGateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, targetQubit, type);
}

__global__ void controlledPhaseGateKernel(MultiQubit multiQubit, const int idQubit1, const int idQubit2)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;

    stateVecSize = multiQubit.numAmpsDividedByNumChunks;
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    bit1 = extractBit (idQubit1, index);
    bit2 = extractBit (idQubit2, index);
    if (bit1 && bit2) {
        stateVecReal [index] = - stateVecReal [index];
        stateVecImag [index] = - stateVecImag [index];
    }
}

void controlledPhaseGate(MultiQubit multiQubit, const int idQubit1, const int idQubit2)
{
    QuESTAssert(idQubit1 >= 0 && idQubit1 < multiQubit.numQubits, 2, __func__);
    QuESTAssert(idQubit2 >= 0 && idQubit2 < multiQubit.numQubits, 1, __func__);
    QuESTAssert(idQubit1 != idQubit2, 3, __func__);

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    controlledPhaseGateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, idQubit1, idQubit2);
}

__global__ void multiControlledPhaseGateKernel(MultiQubit multiQubit, long long int mask)
{
    long long int index;
    long long int stateVecSize;

    stateVecSize = multiQubit.numAmpsDividedByNumChunks;
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    if (mask == (mask & index) ){
        stateVecReal [index] = - stateVecReal [index];
        stateVecImag [index] = - stateVecImag [index];
    }
}

void multiControlledPhaseGate(MultiQubit multiQubit, int *controlQubits, int numControlQubits)
{
    QuESTAssert(numControlQubits > 0 && numControlQubits <= multiQubit.numQubits, 4, __func__);

    int threadsPerCUDABlock, CUDABlocks;
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<multiQubit.numQubits)-1, 2, __func__);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    multiControlledPhaseGateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, mask);
}


__global__ void hadamardKernel (MultiQubit multiQubit, const int targetQubit){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    REAL   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    const long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    REAL recRoot2 = 1.0/sqrt(2.0);

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    // store current state vector values in temp variables
    stateRealUp = stateVecReal[indexUp];
    stateImagUp = stateVecImag[indexUp];

    stateRealLo = stateVecReal[indexLo];
    stateImagLo = stateVecImag[indexLo];

    stateVecReal[indexUp] = recRoot2*(stateRealUp + stateRealLo);
    stateVecImag[indexUp] = recRoot2*(stateImagUp + stateImagLo);

    stateVecReal[indexLo] = recRoot2*(stateRealUp - stateRealLo);
    stateVecImag[indexLo] = recRoot2*(stateImagUp - stateImagLo);
}

void hadamard(MultiQubit multiQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    hadamardKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, targetQubit);
}

__global__ void controlledNotKernel(MultiQubit multiQubit, const int controlQubit, const int targetQubit)
{
    long long int index;
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    long long int stateVecSize;
    int controlBit;

    // ----- temp variables
    REAL   stateRealUp,                             // storage for previous state values
           stateImagUp;                             // (used in updates)
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block
    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    stateVecSize = multiQubit.numAmpsDividedByNumChunks;
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=(stateVecSize>>1)) return;
    thisBlock   = index / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + index%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    controlBit = extractBit(controlQubit, indexUp);
    if (controlBit){
        stateRealUp = stateVecReal[indexUp];
        stateImagUp = stateVecImag[indexUp];

        stateVecReal[indexUp] = stateVecReal[indexLo];
        stateVecImag[indexUp] = stateVecImag[indexLo];

        stateVecReal[indexLo] = stateRealUp;
        stateVecImag[indexLo] = stateImagUp;
    }
}

void controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks)/threadsPerCUDABlock);
    controlledNotKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, controlQubit, targetQubit);
}

__device__ __host__ unsigned int log2Int( unsigned int x )
{
    unsigned int ans = 0 ;
    while( x>>=1 ) ans++;
    return ans ;
}

__device__ void reduceBlock(REAL *arrayIn, REAL *reducedArray, int length){
    int i, l, r;
    int threadMax, maxDepth;
    threadMax = length/2;
    maxDepth = log2Int(length/2);

    for (i=0; i<maxDepth+1; i++){
        if (threadIdx.x<threadMax){
            l = threadIdx.x;
            r = l + threadMax;
            arrayIn[l] = arrayIn[r] + arrayIn[l];
        }
        threadMax = threadMax >> 1;
        __syncthreads(); // optimise -- use warp shuffle instead
    }

    if (threadIdx.x==0) reducedArray[blockIdx.x] = arrayIn[0];
}

__global__ void copySharedReduceBlock(REAL*arrayIn, REAL *reducedArray, int length){
    extern __shared__ REAL tempReductionArray[];
    int blockOffset = blockIdx.x*length;
    tempReductionArray[threadIdx.x*2] = arrayIn[blockOffset + threadIdx.x*2];
    tempReductionArray[threadIdx.x*2+1] = arrayIn[blockOffset + threadIdx.x*2+1];
    __syncthreads();
    reduceBlock(tempReductionArray, reducedArray, length);
}

__global__ void findProbabilityOfZeroKernel(MultiQubit multiQubit,
        const int measureQubit, REAL *reducedArray)
{
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         index;                                               // current index for first half block
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;
    // (good for shared memory parallelism)

    extern __shared__ REAL tempReductionArray[];

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    //
    // --- task-based shared-memory parallel implementation
    //

    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock = thisTask / sizeHalfBlock;
    index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    REAL realVal, imagVal;
    realVal = stateVecReal[index];
    imagVal = stateVecImag[index]; 	
    tempReductionArray[threadIdx.x] = realVal*realVal + imagVal*imagVal;
    __syncthreads();

    if (threadIdx.x<blockDim.x/2){
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
    }
}

int getNumReductionLevels(long long int numValuesToReduce, int numReducedPerLevel){
    int levels=0;
    while (numValuesToReduce){
        numValuesToReduce = numValuesToReduce/numReducedPerLevel;
        levels++;
    }
    return levels;
}

void swapDouble(REAL **a, REAL **b){
    REAL *temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

REAL findProbabilityOfZero(MultiQubit multiQubit,
        const int measureQubit)
{
    long long int numValuesToReduce = multiQubit.numAmpsDividedByNumChunks>>1;
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    REAL stateProb=0;
    int firstTime=1;
    int maxReducedPerLevel = REDUCE_SHARED_SIZE;

    while(numValuesToReduce>1){	
        if (numValuesToReduce<maxReducedPerLevel){
            // Need less than one CUDA block to reduce values
            valuesPerCUDABlock = numValuesToReduce;
            numCUDABlocks = 1;
        } else {
            // Use full CUDA blocks, with block size constrained by shared mem usage
            valuesPerCUDABlock = maxReducedPerLevel;
            numCUDABlocks = ceil((REAL)numValuesToReduce/valuesPerCUDABlock);
        }
        sharedMemSize = valuesPerCUDABlock*sizeof(REAL);

        if (firstTime){
            findProbabilityOfZeroKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                    multiQubit, measureQubit, multiQubit.firstLevelReduction);
            firstTime=0;
        } else {
            cudaDeviceSynchronize();	
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    multiQubit.firstLevelReduction, 
                    multiQubit.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();	
            swapDouble(&(multiQubit.firstLevelReduction), &(multiQubit.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    cudaMemcpy(&stateProb, multiQubit.firstLevelReduction, sizeof(REAL), cudaMemcpyDeviceToHost);
    return stateProb;
}

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 1, __func__);
    REAL stateProb=0;
    stateProb = findProbabilityOfZero(multiQubit, measureQubit);
    if (outcome==1) stateProb = 1.0 - stateProb;
    return stateProb;
}

__global__ void collapseToOutcomeKernel(MultiQubit multiQubit, int measureQubit, REAL totalProbability, int outcome)
{
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         index;                                               // current index for first half block
    // ----- measured probability
    REAL   renorm;                                    // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    // (good for shared memory parallelism)
    long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    // ---------------------------------------------------------------- //
    //            tests                                                 //
    // ---------------------------------------------------------------- //

    //! fix -- this should report an error
    if (!(measureQubit >= 0 && measureQubit < multiQubit.numQubits)) return;
    if (!(totalProbability != 0)) return;
    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    //
    // --- task-based shared-memory parallel implementation
    //
    renorm=1/sqrt(totalProbability);
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    thisBlock = thisTask / sizeHalfBlock;
    index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;

    if (outcome==0){
        stateVecReal[index]=stateVecReal[index]*renorm;
        stateVecImag[index]=stateVecImag[index]*renorm;

        stateVecReal[index+sizeHalfBlock]=0;
        stateVecImag[index+sizeHalfBlock]=0;
    } else if (outcome==1){
        stateVecReal[index]=0;
        stateVecImag[index]=0;

        stateVecReal[index+sizeHalfBlock]=stateVecReal[index+sizeHalfBlock]*renorm;
        stateVecImag[index+sizeHalfBlock]=stateVecImag[index+sizeHalfBlock]*renorm;
    }

}

REAL collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{        
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert((outcome==0 || outcome==1), 10, __func__);
    REAL stateProb;
    stateProb = findProbabilityOfOutcome(multiQubit, measureQubit, outcome);
    QuESTAssert(fabs(stateProb)>REAL_EPS, 8, __func__);

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    if (stateProb!=0) collapseToOutcomeKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, measureQubit, stateProb, outcome);
    return stateProb;
}

__global__ void measureInZeroKernel(MultiQubit multiQubit, int measureQubit, REAL totalProbability)
{
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         index;                                               // current index for first half block
    // ----- measured probability
    REAL   renorm;                                    // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    // (good for shared memory parallelism)
    long long int numTasks=multiQubit.numAmpsDividedByNumChunks>>1;

    // ---------------------------------------------------------------- //
    //            tests                                                 //
    // ---------------------------------------------------------------- //
    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    //
    // --- task-based shared-memory parallel implementation
    //
    renorm=1/sqrt(totalProbability);
    REAL *stateVecReal = multiQubit.deviceStateVec.real;
    REAL *stateVecImag = multiQubit.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    thisBlock = thisTask / sizeHalfBlock;
    index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    stateVecReal[index]=stateVecReal[index]*renorm;
    stateVecImag[index]=stateVecImag[index]*renorm;

    stateVecReal[index+sizeHalfBlock]=0;
    stateVecImag[index+sizeHalfBlock]=0;
}

REAL measureInZero(MultiQubit multiQubit, const int measureQubit)
{        
    REAL stateProb;
    stateProb = findProbabilityOfZero(multiQubit, measureQubit);

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    if (stateProb!=0) measureInZeroKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, measureQubit, stateProb);
    return stateProb;
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

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(multiQubit.numAmpsDividedByNumChunks>>1)/threadsPerCUDABlock);
    collapseToOutcomeKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, measureQubit, stateProbInternal, outcome);

    *stateProb = stateProbInternal;
    return outcome;
}

void exitWithError(int errorCode, const char* func){
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, errorCodes[errorCode]);
    printf("!!!\n");
    printf("exiting..\n");
    exit(errorCode);
}

void QuESTAssert(int isValid, int errorCode, const char* func){
    if (!isValid) exitWithError(errorCode, func);
}

#ifdef __cplusplus
}
#endif
