// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * An implementation of the pure backend in ../QuEST_ops_pure.h for a GPU environment.
 */

# include "../QuEST.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# define REDUCE_SHARED_SIZE 512
# define DEBUG 0


static __device__ int extractBit (int locationOfBitFromRight, long long int theEncodedNumber)
{
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

#ifdef __cplusplus
extern "C" {
#endif

/** works for both statevectors and density matrices */
void statevec_cloneQubitRegister(QubitRegister targetQureg, QubitRegister copyQureg) {
    
    // copy copyQureg's GPU statevec to targetQureg's GPU statevec
    cudaMemcpy(
        copyQureg.deviceStateVec.real, 
        targetQureg.deviceStateVec.real, 
        targetQureg.numAmpsPerChunk*sizeof(*(targetQureg.deviceStateVec.real)), 
        cudaMemcpyDeviceToDevice);
    cudaMemcpy(
        copyQureg.deviceStateVec.imag, 
        targetQureg.deviceStateVec.imag, 
        targetQureg.numAmpsPerChunk*sizeof(*(targetQureg.deviceStateVec.imag)), 
        cudaMemcpyDeviceToDevice);
}

__global__ void densmatr_initPureStateKernel(
    long long int numPureAmps,
    REAL *targetVecReal, REAL *targetVecImag, 
    REAL *copyVecReal, REAL *copyVecImag) 
{
    // this is a particular index of the pure copyQureg
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=numPureAmps) return;
    
    REAL realRow = copyVecReal[index];
    REAL imagRow = copyVecImag[index];
    for (long long int col=0; col < numPureAmps; col++) {
        REAL realCol =   copyVecReal[col];
        REAL imagCol = - copyVecImag[col]; // minus for conjugation
        targetVecReal[col*numPureAmps + index] = realRow*realCol - imagRow*imagCol;
        targetVecImag[col*numPureAmps + index] = realRow*imagCol + imagRow*realCol;
    }
}

void densmatr_initPureState(QubitRegister targetQureg, QubitRegister copyQureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(copyQureg.numAmpsPerChunk)/threadsPerCUDABlock);
    densmatr_initPureStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        copyQureg.numAmpsPerChunk,
        targetQureg.deviceStateVec.real, targetQureg.deviceStateVec.imag,
        copyQureg.deviceStateVec.real,   copyQureg.deviceStateVec.imag);
}

__global__ void densmatr_initStatePlusKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    REAL probFactor = 1.0/((REAL)stateVecSize);
    stateVecReal[index] = probFactor;
    stateVecImag[index] = 0.0;
}

void densmatr_initStatePlus(QubitRegister qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    densmatr_initStatePlusKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void densmatr_initClassicalStateKernel(
    long long int densityNumElems, 
    REAL *densityReal, REAL *densityImag, 
    long long int densityInd)
{
    // initialise the state to all zeros
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= densityNumElems) return;
    
    densityReal[index] = 0.0;
    densityImag[index] = 0.0;
    
    if (index==densityInd){
        // classical state has probability 1
        densityReal[densityInd] = 1.0;
        densityImag[densityInd] = 0.0;
    }
}

void densmatr_initClassicalState(QubitRegister qureg, long long int stateInd)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    
    // index of the desired state in the flat density matrix
    long long int densityDim = 1LL << qureg.numQubitsRepresented;
    long long int densityInd = (densityDim + 1)*stateInd;
    
    // identical to pure version
    densmatr_initClassicalStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag, densityInd);
}

void statevec_createQubitRegister(QubitRegister *qureg, int numQubits, QuESTEnv env)
{   
    // allocate CPU memory
    long long int numAmps = 1L << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;
    qureg->stateVec.real = (REAL*) malloc(numAmpsPerRank * sizeof(qureg->stateVec.real));
    qureg->stateVec.imag = (REAL*) malloc(numAmpsPerRank * sizeof(qureg->stateVec.imag));
    if (env.numRanks>1){
        qureg->pairStateVec.real = (REAL*) malloc(numAmpsPerRank * sizeof(qureg->pairStateVec.real));
        qureg->pairStateVec.imag = (REAL*) malloc(numAmpsPerRank * sizeof(qureg->pairStateVec.imag));
    }

    // check cpu memory allocation was successful
    if ( (!(qureg->stateVec.real) || !(qureg->stateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!\n");
        exit (EXIT_FAILURE);
    }
    if ( env.numRanks>1 && (!(qureg->pairStateVec.real) || !(qureg->pairStateVec.imag))
            && numAmpsPerRank ) {
        printf("Could not allocate memory!\n");
        exit (EXIT_FAILURE);
    }

    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsPerChunk = numAmpsPerRank;
    qureg->numAmpsTotal = numAmps;
    qureg->chunkId = env.rank;
    qureg->numChunks = env.numRanks;
    qureg->isDensityMatrix = 0;

    // allocate GPU memory
    cudaMalloc(&(qureg->deviceStateVec.real), qureg->numAmpsPerChunk*sizeof(*(qureg->deviceStateVec.real)));
    cudaMalloc(&(qureg->deviceStateVec.imag), qureg->numAmpsPerChunk*sizeof(*(qureg->deviceStateVec.imag)));
    cudaMalloc(&(qureg->firstLevelReduction), ceil(qureg->numAmpsPerChunk/(REAL)REDUCE_SHARED_SIZE)*sizeof(REAL));
    cudaMalloc(&(qureg->secondLevelReduction), ceil(qureg->numAmpsPerChunk/(REAL)(REDUCE_SHARED_SIZE*REDUCE_SHARED_SIZE))*
            sizeof(REAL));

    // check gpu memory allocation was successful
    if (!(qureg->deviceStateVec.real) || !(qureg->deviceStateVec.imag)){
        printf("Could not allocate memory on GPU!\n");
        exit (EXIT_FAILURE);
    }

}

void statevec_destroyQubitRegister(QubitRegister qureg, QuESTEnv env)
{
    // Free CPU memory
    free(qureg.stateVec.real);
    free(qureg.stateVec.imag);
    if (env.numRanks>1){
        free(qureg.pairStateVec.real);
        free(qureg.pairStateVec.imag);
    }

    // Free GPU memory
    cudaFree(qureg.deviceStateVec.real);
    cudaFree(qureg.deviceStateVec.imag);
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
    
    seedQuESTDefault();
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

void getEnvironmentString(QuESTEnv env, QubitRegister qureg, char str[200]){
    sprintf(str, "%dqubits_GPU_noMpi_noOMP", qureg.numQubitsInStateVec);    
}

void copyStateToGPU(QubitRegister qureg)
{
    if (DEBUG) printf("Copying data to GPU\n");
    cudaMemcpy(qureg.deviceStateVec.real, qureg.stateVec.real, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyHostToDevice);
    cudaMemcpy(qureg.deviceStateVec.real, qureg.stateVec.real, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyHostToDevice);
    cudaMemcpy(qureg.deviceStateVec.imag, qureg.stateVec.imag, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyHostToDevice);
    cudaMemcpy(qureg.deviceStateVec.imag, qureg.stateVec.imag, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyHostToDevice);
    if (DEBUG) printf("Finished copying data to GPU\n");
}

void copyStateFromGPU(QubitRegister qureg)
{
    cudaDeviceSynchronize();
    if (DEBUG) printf("Copying data from GPU\n");
    cudaMemcpy(qureg.stateVec.real, qureg.deviceStateVec.real, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyDeviceToHost);
    cudaMemcpy(qureg.stateVec.imag, qureg.deviceStateVec.imag, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyDeviceToHost);
    if (DEBUG) printf("Finished copying data from GPU\n");
}

/** Print the current state vector of probability amplitudes for a set of qubits to standard out. 
  For debugging purposes. Each rank should print output serially. Only print output for systems <= 5 qubits
 */
void statevec_reportStateToScreen(QubitRegister qureg, QuESTEnv env, int reportRank){
    long long int index;
    int rank;
    copyStateFromGPU(qureg); 
    if (qureg.numQubitsInStateVec<=5){
        for (rank=0; rank<qureg.numChunks; rank++){
            if (qureg.chunkId==rank){
                if (reportRank) {
                    printf("Reporting state from rank %d [\n", qureg.chunkId);
                    //printf("\trank, index, real, imag\n");
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
    }
}

REAL statevec_getRealAmpEl(QubitRegister qureg, long long int index){
    REAL el=0;
    cudaMemcpy(&el, &(qureg.deviceStateVec.real[index]), 
            sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyDeviceToHost);
    return el;
}

REAL statevec_getImagAmpEl(QubitRegister qureg, long long int index){
    REAL el=0;
    cudaMemcpy(&el, &(qureg.deviceStateVec.imag[index]), 
            sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyDeviceToHost);
    return el;
}

__global__ void statevec_initStateZeroKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
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

void statevec_initStateZero(QubitRegister qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initStateZeroKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_initStatePlusKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    REAL normFactor = 1.0/sqrt((REAL)stateVecSize);
    stateVecReal[index] = normFactor;
    stateVecImag[index] = 0.0;
}

void statevec_initStatePlus(QubitRegister qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initStatePlusKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_initClassicalStateKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag, long long int stateInd){
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
void statevec_initClassicalState(QubitRegister qureg, long long int stateInd)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initClassicalStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag, stateInd);
}

__global__ void statevec_initStateDebugKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    stateVecReal[index] = (index*2.0)/10.0;
    stateVecImag[index] = (index*2.0+1.0)/10.0;
}

void statevec_initStateDebug(QubitRegister qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initStateDebugKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk,
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_initStateOfSingleQubitKernel(long long int stateVecSize, REAL *stateVecReal, REAL *stateVecImag, int qubitId, int outcome){
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

void statevec_initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg->numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initStateOfSingleQubitKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg->numAmpsPerChunk, qureg->deviceStateVec.real, qureg->deviceStateVec.imag, qubitId, outcome);
}

// returns 1 if successful, else 0
int statevec_initStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env){
    long long int chunkSize, stateVecSize;
    long long int indexInChunk, totalIndex;

    chunkSize = qureg->numAmpsPerChunk;
    stateVecSize = chunkSize*qureg->numChunks;

    REAL *stateVecReal = qureg->stateVec.real;
    REAL *stateVecImag = qureg->stateVec.imag;

    FILE *fp;
    char line[200];

    fp = fopen(filename, "r");
    if (fp == NULL)
        return 0;
    
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
                    sscanf(line, "%lf, %lf", &(stateVecReal[indexInChunk]),
                            &(stateVecImag[indexInChunk]));
                # endif
                indexInChunk += 1;
            }
            totalIndex += 1;
        }
    }
    fclose(fp);
    copyStateToGPU(*qureg);
    
    // indicate success
    return 1;
}

int statevec_compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision){
    REAL diff;
    int chunkSize = mq1.numAmpsPerChunk;

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

__global__ void statevec_compactUnitaryKernel (QubitRegister qureg, const int rotQubit, Complex alpha, Complex beta){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    sizeHalfBlock = 1LL << rotQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;
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

void statevec_compactUnitary(QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_compactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, alpha, beta);
}

__global__ void statevec_controlledCompactUnitaryKernel (QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;
    int controlBit;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;
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

void statevec_controlledCompactUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledCompactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, alpha, beta);
}

__global__ void statevec_unitaryKernel(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_unitary(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_unitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, u);
}

__global__ void statevec_controlledUnitaryKernel(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    int controlBit;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_controlledUnitary(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, u);
}

__global__ void statevec_multiControlledUnitaryKernel(QubitRegister qureg, long long int mask, const int targetQubit, ComplexMatrix2 u){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;


    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_multiControlledUnitary(QubitRegister qureg, int *controlQubits, int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_multiControlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask, targetQubit, u);
}

__global__ void statevec_sigmaXKernel(QubitRegister qureg, const int targetQubit){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_sigmaX(QubitRegister qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_sigmaXKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit);
}

__global__ void statevec_sigmaYKernel(QubitRegister qureg, const int targetQubit, const int conjFac){

    long long int sizeHalfBlock = 1LL << targetQubit;
    long long int sizeBlock     = 2LL * sizeHalfBlock;
    long long int numTasks      = qureg.numAmpsPerChunk >> 1;
    long long int thisTask      = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    
    long long int thisBlock     = thisTask / sizeHalfBlock;
    long long int indexUp       = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    long long int indexLo       = indexUp + sizeHalfBlock;
    REAL  stateRealUp, stateImagUp;

    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;
    stateRealUp = stateVecReal[indexUp];
    stateImagUp = stateVecImag[indexUp];

    // update under +-{{0, -i}, {i, 0}}
    stateVecReal[indexUp] = conjFac * stateVecImag[indexLo];
    stateVecImag[indexUp] = conjFac * -stateVecReal[indexLo];
    stateVecReal[indexLo] = conjFac * -stateImagUp;
    stateVecImag[indexLo] = conjFac * stateRealUp;
}

void statevec_sigmaY(QubitRegister qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_sigmaYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, 1);
}

void statevec_sigmaYConj(QubitRegister qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_sigmaYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, -1);
}

__global__ void statevec_controlledSigmaYKernel(QubitRegister qureg, const int controlQubit, const int targetQubit, const int conjFac)
{
    long long int index;
    long long int sizeBlock, sizeHalfBlock;
    long long int stateVecSize;
    int controlBit;

    REAL   stateRealUp, stateImagUp; 
    long long int thisBlock, indexUp, indexLo;                                     
    sizeHalfBlock = 1LL << targetQubit;
    sizeBlock     = 2LL * sizeHalfBlock;

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=(stateVecSize>>1)) return;
    thisBlock   = index / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + index%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    controlBit = extractBit(controlQubit, indexUp);
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

void statevec_controlledSigmaY(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    int conjFactor = 1;
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledSigmaYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, conjFactor);
}

void statevec_controlledSigmaYConj(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    int conjFactor = 1;
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledSigmaYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, conjFactor);
}

__global__ void statevec_phaseShiftByTermKernel(QubitRegister qureg, const int targetQubit, REAL cosAngle, REAL sinAngle) {

    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, indexUp,indexLo;

    REAL stateRealLo, stateImagLo;             
    long long int thisTask; 
    const long long int numTasks = qureg.numAmpsPerChunk >> 1;

    sizeHalfBlock = 1LL << targetQubit;
    sizeBlock     = 2LL * sizeHalfBlock;

    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    thisBlock   = thisTask / sizeHalfBlock;
    indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    indexLo     = indexUp + sizeHalfBlock;

    stateRealLo = stateVecReal[indexLo];
    stateImagLo = stateVecImag[indexLo];

    stateVecReal[indexLo] = cosAngle*stateRealLo - sinAngle*stateImagLo;
    stateVecImag[indexLo] = sinAngle*stateRealLo + cosAngle*stateImagLo;
}

void statevec_phaseShiftByTerm(QubitRegister qureg, const int targetQubit, Complex term)
{   
    REAL cosAngle = term.real;
    REAL sinAngle = term.imag;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_phaseShiftByTermKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, cosAngle, sinAngle);
}

__global__ void statevec_controlledPhaseShiftKernel(QubitRegister qureg, const int idQubit1, const int idQubit2, REAL cosAngle, REAL sinAngle)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;
    REAL stateRealLo, stateImagLo;

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    bit1 = extractBit (idQubit1, index);
    bit2 = extractBit (idQubit2, index);
    if (bit1 && bit2) {
        stateRealLo = stateVecReal[index];
        stateImagLo = stateVecImag[index];
        
        stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
        stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;
    }
}

void statevec_controlledPhaseShift(QubitRegister qureg, const int idQubit1, const int idQubit2, REAL angle)
{
    REAL cosAngle = cos(angle);
    REAL sinAngle = sin(angle);
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledPhaseShiftKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, idQubit1, idQubit2, cosAngle, sinAngle);
}

__global__ void statevec_multiControlledPhaseShiftKernel(QubitRegister qureg, long long int mask, REAL cosAngle, REAL sinAngle) {
    REAL stateRealLo, stateImagLo;
    long long int index;
    long long int stateVecSize;

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;
    
    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    if (mask == (mask & index) ){
        stateRealLo = stateVecReal[index];
        stateImagLo = stateVecImag[index];
        stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
        stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;
    }
}

void statevec_multiControlledPhaseShift(QubitRegister qureg, int *controlQubits, int numControlQubits, REAL angle)
{   
    REAL cosAngle = cos(angle);
    REAL sinAngle = sin(angle);

    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) 
        mask = mask | (1LL<<controlQubits[i]);
        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledPhaseShiftKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask, cosAngle, sinAngle);
}

REAL densmatr_calcTotalProbability(QubitRegister qureg) {
    
    // computes the trace using Kahan summation
    REAL pTotal=0;
    REAL y, t, c;
    c = 0;
    
    long long int numCols = 1LL << qureg.numQubitsRepresented;
    long long diagIndex;
    
    copyStateFromGPU(qureg);
    
    for (int col=0; col< numCols; col++) {
        diagIndex = col*(numCols + 1);
        y = qureg.stateVec.real[diagIndex] - c;
        t = pTotal + y;
        c = ( t - pTotal ) - y; // brackets are important
        pTotal = t;
    }
    
    return pTotal;
}

REAL statevec_calcTotalProbability(QubitRegister qureg){
    /* IJB - implemented using Kahan summation for greater accuracy at a slight floating
       point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    /* Don't change the bracketing in this routine! */
    REAL pTotal=0;
    REAL y, t, c;
    long long int index;
    long long int numAmpsPerRank = qureg.numAmpsPerChunk;

    copyStateFromGPU(qureg);

    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){
        /* Perform pTotal+=qureg.stateVec.real[index]*qureg.stateVec.real[index]; by Kahan */
        // pTotal+=qureg.stateVec.real[index]*qureg.stateVec.real[index];
        y = qureg.stateVec.real[index]*qureg.stateVec.real[index] - c;
        t = pTotal + y;
        c = ( t - pTotal ) - y;
        pTotal = t;

        /* Perform pTotal+=qureg.stateVec.imag[index]*qureg.stateVec.imag[index]; by Kahan */
        //pTotal+=qureg.stateVec.imag[index]*qureg.stateVec.imag[index];
        y = qureg.stateVec.imag[index]*qureg.stateVec.imag[index] - c;
        t = pTotal + y;
        c = ( t - pTotal ) - y;
        pTotal = t;


    }
    return pTotal;
}

__global__ void statevec_controlledPhaseFlipKernel(QubitRegister qureg, const int idQubit1, const int idQubit2)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    bit1 = extractBit (idQubit1, index);
    bit2 = extractBit (idQubit2, index);
    if (bit1 && bit2) {
        stateVecReal [index] = - stateVecReal [index];
        stateVecImag [index] = - stateVecImag [index];
    }
}

void statevec_controlledPhaseFlip(QubitRegister qureg, const int idQubit1, const int idQubit2)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledPhaseFlipKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, idQubit1, idQubit2);
}

__global__ void statevec_multiControlledPhaseFlipKernel(QubitRegister qureg, long long int mask)
{
    long long int index;
    long long int stateVecSize;

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    if (mask == (mask & index) ){
        stateVecReal [index] = - stateVecReal [index];
        stateVecImag [index] = - stateVecImag [index];
    }
}

void statevec_multiControlledPhaseFlip(QubitRegister qureg, int *controlQubits, int numControlQubits)
{
    int threadsPerCUDABlock, CUDABlocks;
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledPhaseFlipKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask);
}


__global__ void statevec_hadamardKernel (QubitRegister qureg, const int targetQubit){
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
    const long long int numTasks=qureg.numAmpsPerChunk>>1;

    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    // ---------------------------------------------------------------- //
    //            rotate                                                //
    // ---------------------------------------------------------------- //

    //! fix -- no necessary for GPU version
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_hadamard(QubitRegister qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_hadamardKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit);
}

__global__ void statevec_controlledNotKernel(QubitRegister qureg, const int controlQubit, const int targetQubit)
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

    stateVecSize = qureg.numAmpsPerChunk;
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_controlledNot(QubitRegister qureg, const int controlQubit, const int targetQubit)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledNotKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit);
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

__global__ void densmatr_findProbabilityOfZeroKernel(
    QubitRegister qureg, const int measureQubit, REAL *reducedArray
) {
    // run by each thread
    // use of block here refers to contiguous amplitudes where measureQubit = 0, 
    // (then =1) and NOT the CUDA block, which is the partitioning of CUDA threads
    
    long long int densityDim    = 1LL << qureg.numQubitsRepresented;
    long long int numTasks      = densityDim >> 1;
    long long int sizeHalfBlock = 1LL << (measureQubit);
    long long int sizeBlock     = 2LL * sizeHalfBlock;
    
    long long int thisBlock;    // which block this thread is processing
    long long int thisTask;     // which part of the block this thread is processing
    long long int basisIndex;   // index of this thread's computational basis state
    long long int densityIndex; // " " index of |basis><basis| in the flat density matrix
    
    // array of each thread's collected probability, to be summed
    extern __shared__ REAL tempReductionArray[];
    
    // figure out which density matrix prob that this thread is assigned
    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    thisBlock = thisTask / sizeHalfBlock;
    basisIndex = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    densityIndex = (densityDim + 1) * basisIndex;
    
    // record the probability in the CUDA-BLOCK-wide array
    REAL prob = qureg.deviceStateVec.real[densityIndex];   // im[densityIndex] assumed ~ 0
    tempReductionArray[threadIdx.x] = prob;
    
    // sum the probs collected by this CUDA-BLOCK's threads into a per-CUDA-BLOCK array
    __syncthreads();
    if (threadIdx.x<blockDim.x/2){
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
    }
}

__global__ void statevec_findProbabilityOfZeroKernel(
        QubitRegister qureg, const int measureQubit, REAL *reducedArray
) {
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         index;                                               // current index for first half block
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    long long int numTasks=qureg.numAmpsPerChunk>>1;
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

    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

__global__ void densmatr_calcPurityKernel(REAL* vecReal, REAL* vecImag, long long int numAmpsToSum, REAL *reducedArray) {
    
    // figure out which density matrix term this thread is assigned
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= numAmpsToSum) return;
    
    REAL term = vecReal[index]*vecReal[index] + vecImag[index]*vecImag[index];
    
    // array of each thread's collected probability, to be summed
    extern __shared__ REAL tempReductionArray[];
    tempReductionArray[threadIdx.x] = term;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

/** Computes the trace of the density matrix squared */
REAL densmatr_calcPurity(QubitRegister qureg) {
    
    // we're summing the square of every term in the density matrix
    long long int numValuesToReduce = qureg.numAmpsPerChunk;
    
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    int maxReducedPerLevel = REDUCE_SHARED_SIZE;
    int firstTime = 1;
    
    while (numValuesToReduce > 1) {
        
        // need less than one CUDA-BLOCK to reduce
        if (numValuesToReduce < maxReducedPerLevel) {
            valuesPerCUDABlock = numValuesToReduce;
            numCUDABlocks = 1;
        }
        // otherwise use only full CUDA-BLOCKS
        else {
            valuesPerCUDABlock = maxReducedPerLevel; // constrained by shared memory
            numCUDABlocks = ceil((REAL)numValuesToReduce/valuesPerCUDABlock);
        }
        // dictates size of reduction array
        sharedMemSize = valuesPerCUDABlock*sizeof(REAL);
        
        // spawn threads to sum the probs in each block
        if (firstTime) {
             densmatr_calcPurityKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                 qureg.deviceStateVec.real, qureg.deviceStateVec.imag, 
                 numValuesToReduce, qureg.firstLevelReduction);
            firstTime = 0;
            
        // sum the block probs
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    REAL traceDensSquared;
    cudaMemcpy(&traceDensSquared, qureg.firstLevelReduction, sizeof(REAL), cudaMemcpyDeviceToHost);
    return traceDensSquared;
}

REAL densmatr_findProbabilityOfZero(QubitRegister qureg, const int measureQubit)
{
    long long int densityDim = 1LL << qureg.numQubitsRepresented;
    long long int numValuesToReduce = densityDim >> 1;  // half of the diagonal has measureQubit=0
    
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    int maxReducedPerLevel = REDUCE_SHARED_SIZE;
    int firstTime = 1;
    
    while (numValuesToReduce > 1) {
        
        // need less than one CUDA-BLOCK to reduce
        if (numValuesToReduce < maxReducedPerLevel) {
            valuesPerCUDABlock = numValuesToReduce;
            numCUDABlocks = 1;
        }
        // otherwise use only full CUDA-BLOCKS
        else {
            valuesPerCUDABlock = maxReducedPerLevel; // constrained by shared memory
            numCUDABlocks = ceil((REAL)numValuesToReduce/valuesPerCUDABlock);
        }
        
        sharedMemSize = valuesPerCUDABlock*sizeof(REAL);
        
        // spawn threads to sum the probs in each block
        if (firstTime) {
            densmatr_findProbabilityOfZeroKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                qureg, measureQubit, qureg.firstLevelReduction);
            firstTime = 0;
            
        // sum the block probs
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    REAL zeroProb;
    cudaMemcpy(&zeroProb, qureg.firstLevelReduction, sizeof(REAL), cudaMemcpyDeviceToHost);
    return zeroProb;
}

REAL statevec_findProbabilityOfZero(QubitRegister qureg, const int measureQubit)
{
    long long int numValuesToReduce = qureg.numAmpsPerChunk>>1;
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
            statevec_findProbabilityOfZeroKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                    qureg, measureQubit, qureg.firstLevelReduction);
            firstTime=0;
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    cudaMemcpy(&stateProb, qureg.firstLevelReduction, sizeof(REAL), cudaMemcpyDeviceToHost);
    return stateProb;
}

REAL statevec_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    REAL outcomeProb = statevec_findProbabilityOfZero(qureg, measureQubit);
    if (outcome==1)
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

REAL densmatr_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome)
{
    REAL outcomeProb = densmatr_findProbabilityOfZero(qureg, measureQubit);
    if (outcome==1) 
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

__global__ void statevec_collapseToKnownProbOutcomeKernel(QubitRegister qureg, int measureQubit, int outcome, REAL totalProbability)
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;

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
    REAL *stateVecReal = qureg.deviceStateVec.real;
    REAL *stateVecImag = qureg.deviceStateVec.imag;

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

/*
 * outcomeProb must accurately be the probability of that qubit outcome in the state-vector, or
 * else the state-vector will lose normalisation
 */
void statevec_collapseToKnownProbOutcome(QubitRegister qureg, const int measureQubit, int outcome, REAL outcomeProb)
{        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((REAL)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_collapseToKnownProbOutcomeKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, measureQubit, outcome, outcomeProb);
}

/** Maps thread ID to a |..0..><..0..| state and then locates |0><1|, |1><0| and |1><1| */
__global__ void densmatr_collapseToKnownProbOutcomeKernel(
    REAL outcomeProb, REAL* vecReal, REAL *vecImag, long long int numBasesToVisit,
    long long int part1, long long int part2, long long int part3, 
    long long int rowBit, long long int colBit, long long int desired, long long int undesired) 
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= numBasesToVisit) return;
    
    long long int base = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2);
    
    // renormalise desired outcome
    vecReal[base + desired] /= outcomeProb;
    vecImag[base + desired] /= outcomeProb;
    
    // kill undesired outcome
    vecReal[base + undesired] = 0;
    vecImag[base + undesired] = 0;
    
    // kill |..0..><..1..| states
    vecReal[base + colBit] = 0;
    vecImag[base + colBit] = 0;
    vecReal[base + rowBit] = 0;
    vecImag[base + rowBit] = 0;
}

/** This involves finding |...i...><...j...| states and killing those where i!=j */
void densmatr_collapseToKnownProbOutcome(QubitRegister qureg, const int measureQubit, int outcome, REAL outcomeProb) {
    
	int rowQubit = measureQubit + qureg.numQubitsRepresented;
    
    int colBit = 1LL << measureQubit;
    int rowBit = 1LL << rowQubit;

    long long int numBasesToVisit = qureg.numAmpsPerChunk/4;
	long long int part1 = colBit -1;	
	long long int part2 = (rowBit >> 1) - colBit;
	long long int part3 = numBasesToVisit - (rowBit >> 1);
    
    long long int desired, undesired;
    if (outcome == 0) {
        desired = 0;
        undesired = colBit | rowBit;
    } else {
        desired = colBit | rowBit;
        undesired = 0;
    }
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numBasesToVisit / (REAL) threadsPerCUDABlock);
    densmatr_collapseToKnownProbOutcomeKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        outcomeProb, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numBasesToVisit,
        part1, part2, part3, rowBit, colBit, desired, undesired);
}

__global__ void densmatr_combineDensityMatricesKernel(REAL combineProb, QubitRegister combineQureg, REAL otherProb, QubitRegister otherQureg, long long int numAmpsToVisit) {
    
    long long int ampInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (ampInd >= numAmpsToVisit) return;
    
    combineQureg.deviceStateVec.real[ampInd] *= combineProb;
    combineQureg.deviceStateVec.imag[ampInd] *= combineProb;
  
    combineQureg.deviceStateVec.real[ampInd] += otherProb*otherQureg.deviceStateVec.real[ampInd];
    combineQureg.deviceStateVec.imag[ampInd] += otherProb*otherQureg.deviceStateVec.imag[ampInd];
}

void densmatr_combineDensityMatrices(REAL combineProb, QubitRegister combineQureg, REAL otherProb, QubitRegister otherQureg) {
    
    long long int numAmpsToVisit = combineQureg.numAmpsPerChunk;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (REAL) threadsPerCUDABlock);
    densmatr_combineDensityMatricesKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        combineProb, combineQureg, otherProb, otherQureg, numAmpsToVisit
    );
}

/** Called once for every 4 amplitudes in density matrix 
 * Works by establishing the |..0..><..0..| state (for its given index) then 
 * visiting |..1..><..0..| and |..0..><..1..|. Labels |part1 X pa><rt2 NOT(X) part3|
 * From the brain of Simon Benjamin
 */
__global__ void densmatr_oneQubitDephaseKernel(
    REAL fac, REAL* vecReal, REAL *vecImag, long long int numAmpsToVisit,
    long long int part1, long long int part2, long long int part3, 
    long long int colBit, long long int rowBit)
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= numAmpsToVisit) return;
    
    long long int ampInd = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2);
    vecReal[ampInd + colBit] *= fac;
    vecImag[ampInd + colBit] *= fac;
    vecReal[ampInd + rowBit] *= fac;
    vecImag[ampInd + rowBit] *= fac;
}

void densmatr_oneQubitDephase(QubitRegister qureg, const int targetQubit, REAL dephase) {
    
    long long int numAmpsToVisit = qureg.numAmpsPerChunk/4;
    
    int rowQubit = targetQubit + qureg.numQubitsRepresented;
    long long int colBit = 1LL << targetQubit;
    long long int rowBit = 1LL << rowQubit;
    
    long long int part1 = colBit - 1;
    long long int part2 = (rowBit >> 1) - colBit;
    long long int part3 = numAmpsToVisit - (rowBit >> 1);
    REAL dephFac = 1 - dephase;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (REAL) threadsPerCUDABlock);
    densmatr_oneQubitDephaseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        dephFac, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, colBit, rowBit);
}

/** Called 12 times for every 16 amplitudes in density matrix 
 * Each sums from the |..0..0..><..0..0..| index to visit either
 * |..0..0..><..0..1..|,  |..0..0..><..1..0..|,  |..0..0..><..1..1..|,  |..0..1..><..0..0..|
 * etc and so on to |..1..1..><..1..0|. Labels |part1 0 part2 0 par><t3 0 part4 0 part5|.
 * From the brain of Simon Benjamin
 */
__global__ void densmatr_twoQubitDephaseKernel(
    REAL fac, REAL* vecReal, REAL *vecImag, long long int numBackgroundStates, long long int numAmpsToVisit,
    long long int part1, long long int part2, long long int part3, long long int part4, long long int part5,
    long long int colBit1, long long int rowBit1, long long int colBit2, long long int rowBit2) 
{
    long long int outerInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (outerInd >= numAmpsToVisit) return;
    
    // sets meta in 1...14 excluding 5, 10, creating bit string DCBA for |..D..C..><..B..A|
    int meta = 1 + (outerInd/numBackgroundStates);
    if (meta > 4) meta++;
    if (meta > 9) meta++;
    
    long long int shift = rowBit2*((meta>>3)%2) + rowBit1*((meta>>2)%2) + colBit2*((meta>>1)%2) + colBit1*(meta%2);
    long long int scanInd = outerInd % numBackgroundStates;
    long long int stateInd = (
        shift + 
        (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2) + ((scanInd&part4)<<3) + ((scanInd&part5)<<4));
    
    vecReal[stateInd] *= fac;
    vecImag[stateInd] *= fac;
}

// @TODO is separating these 12 amplitudes really faster than letting every 16th base modify 12 elems?
void densmatr_twoQubitDephase(QubitRegister qureg, int qubit1, int qubit2, REAL dephase) {
    if (dephase == 0)
        return;
    
    // ensure qubit2 is further left than qubit1
    if (qubit1 > qubit2)  {
        int tmp = qubit1;
        qubit1 = qubit2;
        qubit2 = tmp;
    }
    
    int rowQubit1 = qubit1 + qureg.numQubitsRepresented;
    int rowQubit2 = qubit2 + qureg.numQubitsRepresented;
    
    long long int colBit1 = 1LL << qubit1;
    long long int rowBit1 = 1LL << rowQubit1;
    long long int colBit2 = 1LL << qubit2;
    long long int rowBit2 = 1LL << rowQubit2;
    
    long long int part1 = colBit1 - 1;
    long long int part2 = (colBit2 >> 1) - colBit1;
    long long int part3 = (rowBit1 >> 2) - (colBit2 >> 1);
    long long int part4 = (rowBit2 >> 3) - (rowBit1 >> 2);
    long long int part5 = (qureg.numAmpsPerChunk/16) - (rowBit2 >> 3);
    REAL dephFac = 1 - dephase;
    
    // refers to states |a 0 b 0 c><d 0 e 0 f| (target qubits are fixed)
    long long int numBackgroundStates = qureg.numAmpsPerChunk/16;
    
    // 12 of these states experience dephasing
    long long int numAmpsToVisit = 12 * numBackgroundStates;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (REAL) threadsPerCUDABlock);
    densmatr_twoQubitDephaseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        dephFac, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numBackgroundStates, numAmpsToVisit,
        part1, part2, part3, part4, part5, colBit1, rowBit1, colBit2, rowBit2);
}

/** Works like oneQubitDephase but modifies every other element, and elements are averaged in pairs */
__global__ void densmatr_oneQubitDepolariseKernel(
    REAL depolLevel, REAL* vecReal, REAL *vecImag, long long int numAmpsToVisit,
    long long int part1, long long int part2, long long int part3, 
    long long int bothBits)
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= numAmpsToVisit) return;
    
    long long int baseInd = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2);
    long long int targetInd = baseInd + bothBits;
    
    REAL realAvDepol = depolLevel * 0.5 * (vecReal[baseInd] + vecReal[targetInd]);
    REAL imagAvDepol = depolLevel * 0.5 * (vecImag[baseInd] + vecImag[targetInd]);
    
    vecReal[baseInd]   *= 1 - depolLevel;
    vecImag[baseInd]   *= 1 - depolLevel;
    vecReal[targetInd] *= 1 - depolLevel;
    vecImag[targetInd] *= 1 - depolLevel;
    
    vecReal[baseInd]   += realAvDepol;
    vecImag[baseInd]   += imagAvDepol;
    vecReal[targetInd] += realAvDepol;
    vecImag[targetInd] += imagAvDepol;
}

void densmatr_oneQubitDepolarise(QubitRegister qureg, const int targetQubit, REAL depolLevel) {
    
    if (depolLevel == 0)
        return;
    
    densmatr_oneQubitDephase(qureg, targetQubit, depolLevel);
    
    long long int numAmpsToVisit = qureg.numAmpsPerChunk/4;
    int rowQubit = targetQubit + qureg.numQubitsRepresented;
    
    long long int colBit = 1LL << targetQubit;
    long long int rowBit = 1LL << rowQubit;
    long long int bothBits = colBit | rowBit;
    
    long long int part1 = colBit - 1;
    long long int part2 = (rowBit >> 1) - colBit;
    long long int part3 = numAmpsToVisit - (rowBit >> 1);
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (REAL) threadsPerCUDABlock);
    densmatr_oneQubitDepolariseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        depolLevel, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, bothBits);
}

/** Called once for every 16 amplitudes */
__global__ void densmatr_twoQubitDepolariseKernel(
    REAL depolLevel, REAL* vecReal, REAL *vecImag, long long int numAmpsToVisit,
    long long int part1, long long int part2, long long int part3, 
    long long int part4, long long int part5,
    long long int rowCol1, long long int rowCol2)
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= numAmpsToVisit) return;
    
    // index of |..0..0..><..0..0|
    long long int ind00 = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2) + ((scanInd&part4)<<3) + ((scanInd&part5)<<4);
    long long int ind01 = ind00 + rowCol1;
    long long int ind10 = ind00 + rowCol2;
    long long int ind11 = ind00 + rowCol1 + rowCol2;
    
    REAL realAvDepol = depolLevel * 0.25 * (
        vecReal[ind00] + vecReal[ind01] + vecReal[ind10] + vecReal[ind11]);
    REAL imagAvDepol = depolLevel * 0.25 * (
        vecImag[ind00] + vecImag[ind01] + vecImag[ind10] + vecImag[ind11]);
    
    REAL retain = 1 - depolLevel;
    vecReal[ind00] *= retain; vecImag[ind00] *= retain;
    vecReal[ind01] *= retain; vecImag[ind01] *= retain;
    vecReal[ind10] *= retain; vecImag[ind10] *= retain;
    vecReal[ind11] *= retain; vecImag[ind11] *= retain;

    vecReal[ind00] += realAvDepol; vecImag[ind00] += imagAvDepol;
    vecReal[ind01] += realAvDepol; vecImag[ind01] += imagAvDepol;
    vecReal[ind10] += realAvDepol; vecImag[ind10] += imagAvDepol;
    vecReal[ind11] += realAvDepol; vecImag[ind11] += imagAvDepol;
}

void densmatr_twoQubitDepolarise(QubitRegister qureg, int qubit1, int qubit2, REAL depolLevel) {
    
    if (depolLevel == 0)
        return;
    
    // ensure qubit2 is further left than qubit1
    if (qubit1 > qubit2)  {
        int tmp = qubit1;
        qubit1 = qubit2;
        qubit2 = tmp;
    }
    
    densmatr_twoQubitDephase(qureg, qubit1, qubit2, depolLevel);
    
    int rowQubit1 = qubit1 + qureg.numQubitsRepresented;
    int rowQubit2 = qubit2 + qureg.numQubitsRepresented;
    
    long long int colBit1 = 1LL << qubit1;
    long long int rowBit1 = 1LL << rowQubit1;
    long long int colBit2 = 1LL << qubit2;
    long long int rowBit2 = 1LL << rowQubit2;
    
    long long int rowCol1 = colBit1 | rowBit1;
    long long int rowCol2 = colBit2 | rowBit2;
    
    long long int numAmpsToVisit = qureg.numAmpsPerChunk/16;
    long long int part1 = colBit1 - 1;
    long long int part2 = (colBit2 >> 1) - colBit1;
    long long int part3 = (rowBit1 >> 2) - (colBit2 >> 1);
    long long int part4 = (rowBit2 >> 3) - (rowBit1 >> 2);
    long long int part5 = numAmpsToVisit - (rowBit2 >> 3);
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (REAL) threadsPerCUDABlock);
    densmatr_twoQubitDepolariseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        depolLevel, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, part4, part5, rowCol1, rowCol2);
}

/*

SHOULD REALLY REDUCE A 2*LENGTH ARRAY, not reduce twice
__global__ void statevec_getInnerProductKernel(QubitRegister bra, QubitRegister ket, REAL *reducedArray) {
    
    // work amplitude index
    long long int stateVecSize = bra.numAmpsPerChunk;
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;
    
    REAL braRe = bra.deviceStateVec.real[index];
    REAL braIm = bra.deviceStateVec.imag[index];
    REAL ketRe = ket.deviceStateVec.real[index];
    REAL ketIm = ket.deviceStateVec.imag[index];
    
    // conj(bra_i) * ket_i
    REAL innerProdReal = braRe*ketRe - braIm*ketIm;
    REAL innerProdImag = braRe*ketIm + braIm*ketRe;
    
    // array of each thread's collected innerProdReal and innerProdImag, alternating
    // real in even inds, imag in off inds: this array has length 2*numThreads
    extern __shared__ REAL innerProdReductionArray[];
    innerProdReductionArray[2*thread.x    ] = innerProdReal;
    innerProdReductionArray[2*thread.x + 1] = innerProdImag;
    __syncthreads(); 
}
*/


#ifdef __cplusplus
}
#endif
