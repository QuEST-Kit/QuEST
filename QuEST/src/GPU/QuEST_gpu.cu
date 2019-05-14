// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * An implementation of the backend in ../QuEST_internal.h for a GPU environment.
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_internal.h"    // purely to resolve getQuESTDefaultSeedKey
# include "mt19937ar.h"

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
    





void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps) {
    
    cudaDeviceSynchronize();
    cudaMemcpy(
        qureg.deviceStateVec.real + startInd, 
        reals,
        numAmps * sizeof(*(qureg.deviceStateVec.real)), 
        cudaMemcpyHostToDevice);
    cudaMemcpy(
        qureg.deviceStateVec.imag + startInd,
        imags,
        numAmps * sizeof(*(qureg.deviceStateVec.real)), 
        cudaMemcpyHostToDevice);
}


/** works for both statevectors and density matrices */
void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg) {
    
    // copy copyQureg's GPU statevec to targetQureg's GPU statevec
    cudaDeviceSynchronize();
    cudaMemcpy(
        targetQureg.deviceStateVec.real, 
        copyQureg.deviceStateVec.real, 
        targetQureg.numAmpsPerChunk*sizeof(*(targetQureg.deviceStateVec.real)), 
        cudaMemcpyDeviceToDevice);
    cudaMemcpy(
        targetQureg.deviceStateVec.imag, 
        copyQureg.deviceStateVec.imag, 
        targetQureg.numAmpsPerChunk*sizeof(*(targetQureg.deviceStateVec.imag)), 
        cudaMemcpyDeviceToDevice);
}

__global__ void densmatr_initPureStateKernel(
    long long int numPureAmps,
    qreal *targetVecReal, qreal *targetVecImag, 
    qreal *copyVecReal, qreal *copyVecImag) 
{
    // this is a particular index of the pure copyQureg
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=numPureAmps) return;
    
    qreal realRow = copyVecReal[index];
    qreal imagRow = copyVecImag[index];
    for (long long int col=0; col < numPureAmps; col++) {
        qreal realCol =   copyVecReal[col];
        qreal imagCol = - copyVecImag[col]; // minus for conjugation
        targetVecReal[col*numPureAmps + index] = realRow*realCol - imagRow*imagCol;
        targetVecImag[col*numPureAmps + index] = realRow*imagCol + imagRow*realCol;
    }
}

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(copyQureg.numAmpsPerChunk)/threadsPerCUDABlock);
    densmatr_initPureStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        copyQureg.numAmpsPerChunk,
        targetQureg.deviceStateVec.real, targetQureg.deviceStateVec.imag,
        copyQureg.deviceStateVec.real,   copyQureg.deviceStateVec.imag);
}

__global__ void densmatr_initPlusStateKernel(long long int stateVecSize, qreal probFactor, qreal *stateVecReal, qreal *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    stateVecReal[index] = probFactor;
    stateVecImag[index] = 0.0;
}

void densmatr_initPlusState(Qureg qureg)
{
    qreal probFactor = 1.0/((qreal) (1LL << qureg.numQubitsRepresented));
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    densmatr_initPlusStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        probFactor,
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void densmatr_initClassicalStateKernel(
    long long int densityNumElems, 
    qreal *densityReal, qreal *densityImag, 
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

void densmatr_initClassicalState(Qureg qureg, long long int stateInd)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    
    // index of the desired state in the flat density matrix
    long long int densityDim = 1LL << qureg.numQubitsRepresented;
    long long int densityInd = (densityDim + 1)*stateInd;
    
    // identical to pure version
    densmatr_initClassicalStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag, densityInd);
}

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env)
{   
    // allocate CPU memory
    long long int numAmps = 1L << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;
    qureg->stateVec.real = (qreal*) malloc(numAmpsPerRank * sizeof(qureg->stateVec.real));
    qureg->stateVec.imag = (qreal*) malloc(numAmpsPerRank * sizeof(qureg->stateVec.imag));
    if (env.numRanks>1){
        qureg->pairStateVec.real = (qreal*) malloc(numAmpsPerRank * sizeof(qureg->pairStateVec.real));
        qureg->pairStateVec.imag = (qreal*) malloc(numAmpsPerRank * sizeof(qureg->pairStateVec.imag));
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
    cudaMalloc(&(qureg->firstLevelReduction), ceil(qureg->numAmpsPerChunk/(qreal)REDUCE_SHARED_SIZE)*sizeof(qreal));
    cudaMalloc(&(qureg->secondLevelReduction), ceil(qureg->numAmpsPerChunk/(qreal)(REDUCE_SHARED_SIZE*REDUCE_SHARED_SIZE))*
            sizeof(qreal));

    // check gpu memory allocation was successful
    if (!(qureg->deviceStateVec.real) || !(qureg->deviceStateVec.imag)){
        printf("Could not allocate memory on GPU!\n");
        exit (EXIT_FAILURE);
    }

}

void statevec_destroyQureg(Qureg qureg, QuESTEnv env)
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

QuESTEnv createQuESTEnv(void) {
    // init MPI environment
    if (!GPUExists()){
        printf("Trying to run GPU code with no GPU available\n");
        exit(EXIT_FAILURE);
    }
    
    QuESTEnv env;
    env.rank=0;
    env.numRanks=1;
    
    seedQuESTDefault();
    
    return env;
}

void syncQuESTEnv(QuESTEnv env){
    cudaDeviceSynchronize();
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void destroyQuESTEnv(QuESTEnv env){
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

void getEnvironmentString(QuESTEnv env, Qureg qureg, char str[200]){
    sprintf(str, "%dqubits_GPU_noMpi_noOMP", qureg.numQubitsInStateVec);    
}

void copyStateToGPU(Qureg qureg)
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

void copyStateFromGPU(Qureg qureg)
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
void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank){
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

qreal statevec_getRealAmp(Qureg qureg, long long int index){
    qreal el=0;
    cudaMemcpy(&el, &(qureg.deviceStateVec.real[index]), 
            sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyDeviceToHost);
    return el;
}

qreal statevec_getImagAmp(Qureg qureg, long long int index){
    qreal el=0;
    cudaMemcpy(&el, &(qureg.deviceStateVec.imag[index]), 
            sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyDeviceToHost);
    return el;
}

__global__ void statevec_initZeroStateKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag){
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

void statevec_initZeroState(Qureg qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initZeroStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_initPlusStateKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    qreal normFactor = 1.0/sqrt((qreal)stateVecSize);
    stateVecReal[index] = normFactor;
    stateVecImag[index] = 0.0;
}

void statevec_initPlusState(Qureg qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initPlusStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_initClassicalStateKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag, long long int stateInd){
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
void statevec_initClassicalState(Qureg qureg, long long int stateInd)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initClassicalStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag, stateInd);
}

__global__ void statevec_initStateDebugKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    stateVecReal[index] = (index*2.0)/10.0;
    stateVecImag[index] = (index*2.0+1.0)/10.0;
}

void statevec_initStateDebug(Qureg qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initStateDebugKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk,
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_initStateOfSingleQubitKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag, int qubitId, int outcome){
    long long int index;
    int bit;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    qreal normFactor = 1.0/sqrt((qreal)stateVecSize/2);
    bit = extractBit(qubitId, index);
    if (bit==outcome) {
        stateVecReal[index] = normFactor;
        stateVecImag[index] = 0.0;
    } else {
        stateVecReal[index] = 0.0;
        stateVecImag[index] = 0.0;
    }
}

void statevec_initStateOfSingleQubit(Qureg *qureg, int qubitId, int outcome)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg->numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initStateOfSingleQubitKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg->numAmpsPerChunk, qureg->deviceStateVec.real, qureg->deviceStateVec.imag, qubitId, outcome);
}

// returns 1 if successful, else 0
int statevec_initStateFromSingleFile(Qureg *qureg, char filename[200], QuESTEnv env){
    long long int chunkSize, stateVecSize;
    long long int indexInChunk, totalIndex;

    chunkSize = qureg->numAmpsPerChunk;
    stateVecSize = chunkSize*qureg->numChunks;

    qreal *stateVecReal = qureg->stateVec.real;
    qreal *stateVecImag = qureg->stateVec.imag;

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

int statevec_compareStates(Qureg mq1, Qureg mq2, qreal precision){
    qreal diff;
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

__global__ void statevec_compactUnitaryKernel (Qureg qureg, const int rotQubit, Complex alpha, Complex beta){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,stateRealLo,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;
    qreal alphaImag=alpha.imag, alphaReal=alpha.real;
    qreal betaImag=beta.imag, betaReal=beta.real;

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

void statevec_compactUnitary(Qureg qureg, const int targetQubit, Complex alpha, Complex beta) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_compactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, alpha, beta);
}

__global__ void statevec_controlledCompactUnitaryKernel (Qureg qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,stateRealLo,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;
    qreal alphaImag=alpha.imag, alphaReal=alpha.real;
    qreal betaImag=beta.imag, betaReal=beta.real;

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

void statevec_controlledCompactUnitary(Qureg qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledCompactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, alpha, beta);
}

__global__ void statevec_unitaryKernel(Qureg qureg, const int targetQubit, ComplexMatrix2 u){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,stateRealLo,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_unitary(Qureg qureg, const int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_unitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, u);
}

__global__ void statevec_controlledUnitaryKernel(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,stateRealLo,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_controlledUnitary(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, u);
}

__global__ void statevec_multiControlledUnitaryKernel(Qureg qureg, long long int mask, const int targetQubit, ComplexMatrix2 u){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,stateRealLo,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_multiControlledUnitary(Qureg qureg, int *controlQubits, int numControlQubits, const int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_multiControlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask, targetQubit, u);
}

__global__ void statevec_pauliXKernel(Qureg qureg, const int targetQubit){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_pauliX(Qureg qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_pauliXKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit);
}

__global__ void statevec_pauliYKernel(Qureg qureg, const int targetQubit, const int conjFac){

    long long int sizeHalfBlock = 1LL << targetQubit;
    long long int sizeBlock     = 2LL * sizeHalfBlock;
    long long int numTasks      = qureg.numAmpsPerChunk >> 1;
    long long int thisTask      = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    
    long long int thisBlock     = thisTask / sizeHalfBlock;
    long long int indexUp       = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    long long int indexLo       = indexUp + sizeHalfBlock;
    qreal  stateRealUp, stateImagUp;

    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;
    stateRealUp = stateVecReal[indexUp];
    stateImagUp = stateVecImag[indexUp];

    // update under +-{{0, -i}, {i, 0}}
    stateVecReal[indexUp] = conjFac * stateVecImag[indexLo];
    stateVecImag[indexUp] = conjFac * -stateVecReal[indexLo];
    stateVecReal[indexLo] = conjFac * -stateImagUp;
    stateVecImag[indexLo] = conjFac * stateRealUp;
}

void statevec_pauliY(Qureg qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_pauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, 1);
}

void statevec_pauliYConj(Qureg qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_pauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, -1);
}

__global__ void statevec_controlledPauliYKernel(Qureg qureg, const int controlQubit, const int targetQubit, const int conjFac)
{
    long long int index;
    long long int sizeBlock, sizeHalfBlock;
    long long int stateVecSize;
    int controlBit;

    qreal   stateRealUp, stateImagUp; 
    long long int thisBlock, indexUp, indexLo;                                     
    sizeHalfBlock = 1LL << targetQubit;
    sizeBlock     = 2LL * sizeHalfBlock;

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_controlledPauliY(Qureg qureg, const int controlQubit, const int targetQubit)
{
    int conjFactor = 1;
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledPauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, conjFactor);
}

void statevec_controlledPauliYConj(Qureg qureg, const int controlQubit, const int targetQubit)
{
    int conjFactor = -1;
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledPauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, conjFactor);
}

__global__ void statevec_phaseShiftByTermKernel(Qureg qureg, const int targetQubit, qreal cosAngle, qreal sinAngle) {

    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, indexUp,indexLo;

    qreal stateRealLo, stateImagLo;             
    long long int thisTask; 
    const long long int numTasks = qureg.numAmpsPerChunk >> 1;

    sizeHalfBlock = 1LL << targetQubit;
    sizeBlock     = 2LL * sizeHalfBlock;

    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_phaseShiftByTerm(Qureg qureg, const int targetQubit, Complex term)
{   
    qreal cosAngle = term.real;
    qreal sinAngle = term.imag;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_phaseShiftByTermKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, cosAngle, sinAngle);
}

__global__ void statevec_controlledPhaseShiftKernel(Qureg qureg, const int idQubit1, const int idQubit2, qreal cosAngle, qreal sinAngle)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;
    qreal stateRealLo, stateImagLo;

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_controlledPhaseShift(Qureg qureg, const int idQubit1, const int idQubit2, qreal angle)
{
    qreal cosAngle = cos(angle);
    qreal sinAngle = sin(angle);
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledPhaseShiftKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, idQubit1, idQubit2, cosAngle, sinAngle);
}

__global__ void statevec_multiControlledPhaseShiftKernel(Qureg qureg, long long int mask, qreal cosAngle, qreal sinAngle) {
    qreal stateRealLo, stateImagLo;
    long long int index;
    long long int stateVecSize;

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;
    
    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    if (mask == (mask & index) ){
        stateRealLo = stateVecReal[index];
        stateImagLo = stateVecImag[index];
        stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
        stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;
    }
}

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle)
{   
    qreal cosAngle = cos(angle);
    qreal sinAngle = sin(angle);

    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) 
        mask = mask | (1LL<<controlQubits[i]);
        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledPhaseShiftKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask, cosAngle, sinAngle);
}

qreal densmatr_calcTotalProb(Qureg qureg) {
    
    // computes the trace using Kahan summation
    qreal pTotal=0;
    qreal y, t, c;
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

qreal statevec_calcTotalProb(Qureg qureg){
    /* IJB - implemented using Kahan summation for greater accuracy at a slight floating
       point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    /* Don't change the bracketing in this routine! */
    qreal pTotal=0;
    qreal y, t, c;
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

__global__ void statevec_controlledPhaseFlipKernel(Qureg qureg, const int idQubit1, const int idQubit2)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    bit1 = extractBit (idQubit1, index);
    bit2 = extractBit (idQubit2, index);
    if (bit1 && bit2) {
        stateVecReal [index] = - stateVecReal [index];
        stateVecImag [index] = - stateVecImag [index];
    }
}

void statevec_controlledPhaseFlip(Qureg qureg, const int idQubit1, const int idQubit2)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledPhaseFlipKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, idQubit1, idQubit2);
}

__global__ void statevec_multiControlledPhaseFlipKernel(Qureg qureg, long long int mask)
{
    long long int index;
    long long int stateVecSize;

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    if (mask == (mask & index) ){
        stateVecReal [index] = - stateVecReal [index];
        stateVecImag [index] = - stateVecImag [index];
    }
}

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits)
{
    int threadsPerCUDABlock, CUDABlocks;
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledPhaseFlipKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask);
}


__global__ void statevec_hadamardKernel (Qureg qureg, const int targetQubit){
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block

    // ----- temp variables
    qreal   stateRealUp,stateRealLo,                             // storage for previous state values
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

    qreal recRoot2 = 1.0/sqrt(2.0);

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

void statevec_hadamard(Qureg qureg, const int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_hadamardKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit);
}

__global__ void statevec_controlledNotKernel(Qureg qureg, const int controlQubit, const int targetQubit)
{
    long long int index;
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    long long int stateVecSize;
    int controlBit;

    // ----- temp variables
    qreal   stateRealUp,                             // storage for previous state values
           stateImagUp;                             // (used in updates)
    long long int thisBlock,                                           // current block
         indexUp,indexLo;                                     // current index and corresponding index in lower half block
    sizeHalfBlock = 1LL << targetQubit;                               // size of blocks halved
    sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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

void statevec_controlledNot(Qureg qureg, const int controlQubit, const int targetQubit)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledNotKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit);
}

__device__ __host__ unsigned int log2Int( unsigned int x )
{
    unsigned int ans = 0 ;
    while( x>>=1 ) ans++;
    return ans ;
}

__device__ void reduceBlock(qreal *arrayIn, qreal *reducedArray, int length){
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

__global__ void copySharedReduceBlock(qreal*arrayIn, qreal *reducedArray, int length){
    extern __shared__ qreal tempReductionArray[];
    int blockOffset = blockIdx.x*length;
    tempReductionArray[threadIdx.x*2] = arrayIn[blockOffset + threadIdx.x*2];
    tempReductionArray[threadIdx.x*2+1] = arrayIn[blockOffset + threadIdx.x*2+1];
    __syncthreads();
    reduceBlock(tempReductionArray, reducedArray, length);
}

__global__ void densmatr_findProbabilityOfZeroKernel(
    Qureg qureg, const int measureQubit, qreal *reducedArray
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
    extern __shared__ qreal tempReductionArray[];
    
    // figure out which density matrix prob that this thread is assigned
    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    thisBlock = thisTask / sizeHalfBlock;
    basisIndex = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    densityIndex = (densityDim + 1) * basisIndex;
    
    // record the probability in the CUDA-BLOCK-wide array
    qreal prob = qureg.deviceStateVec.real[densityIndex];   // im[densityIndex] assumed ~ 0
    tempReductionArray[threadIdx.x] = prob;
    
    // sum the probs collected by this CUDA-BLOCK's threads into a per-CUDA-BLOCK array
    __syncthreads();
    if (threadIdx.x<blockDim.x/2){
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
    }
}

__global__ void statevec_findProbabilityOfZeroKernel(
        Qureg qureg, const int measureQubit, qreal *reducedArray
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

    extern __shared__ qreal tempReductionArray[];

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

    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

    thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;

    thisBlock = thisTask / sizeHalfBlock;
    index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
    qreal realVal, imagVal;
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

void swapDouble(qreal **a, qreal **b){
    qreal *temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

qreal densmatr_findProbabilityOfZero(Qureg qureg, const int measureQubit)
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
            numCUDABlocks = ceil((qreal)numValuesToReduce/valuesPerCUDABlock);
        }
        
        sharedMemSize = valuesPerCUDABlock*sizeof(qreal);
        
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
    
    qreal zeroProb;
    cudaMemcpy(&zeroProb, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    return zeroProb;
}

qreal statevec_findProbabilityOfZero(Qureg qureg, const int measureQubit)
{
    long long int numValuesToReduce = qureg.numAmpsPerChunk>>1;
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    qreal stateProb=0;
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
            numCUDABlocks = ceil((qreal)numValuesToReduce/valuesPerCUDABlock);
        }
        sharedMemSize = valuesPerCUDABlock*sizeof(qreal);

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
    cudaMemcpy(&stateProb, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    return stateProb;
}

qreal statevec_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome)
{
    qreal outcomeProb = statevec_findProbabilityOfZero(qureg, measureQubit);
    if (outcome==1)
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome)
{
    qreal outcomeProb = densmatr_findProbabilityOfZero(qureg, measureQubit);
    if (outcome==1) 
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}


/** computes either a real or imag term in the inner product */
__global__ void statevec_calcInnerProductKernel(
    int getRealComp,
    qreal* vecReal1, qreal* vecImag1, qreal* vecReal2, qreal* vecImag2, 
    long long int numTermsToSum, qreal* reducedArray) 
{
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= numTermsToSum) return;
    
    // choose whether to calculate the real or imaginary term of the inner product
    qreal innerProdTerm;
    if (getRealComp)
        innerProdTerm = vecReal1[index]*vecReal2[index] + vecImag1[index]*vecImag2[index];
    else
        innerProdTerm = vecReal1[index]*vecImag2[index] - vecImag1[index]*vecReal2[index];
    
    // array of each thread's collected probability, to be summed
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = innerProdTerm;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

/** Terrible code which unnecessarily individually computes and sums the real and imaginary components of the
 * inner product, so as to not have to worry about keeping the sums separated during reduction.
 * Truly disgusting, probably doubles runtime, please fix.
 * @TODO could even do the kernel twice, storing real in bra.reduc and imag in ket.reduc?
 */
Complex statevec_calcInnerProduct(Qureg bra, Qureg ket) {
    
    qreal innerProdReal, innerProdImag;
    
    int getRealComp;
    long long int numValuesToReduce;
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    int maxReducedPerLevel;
    int firstTime;
    
    // compute real component of inner product
    getRealComp = 1;
    numValuesToReduce = bra.numAmpsPerChunk;
    maxReducedPerLevel = REDUCE_SHARED_SIZE;
    firstTime = 1;
    while (numValuesToReduce > 1) {
        if (numValuesToReduce < maxReducedPerLevel) {
            valuesPerCUDABlock = numValuesToReduce;
            numCUDABlocks = 1;
        }
        else {
            valuesPerCUDABlock = maxReducedPerLevel;
            numCUDABlocks = ceil((qreal)numValuesToReduce/valuesPerCUDABlock);
        }
        sharedMemSize = valuesPerCUDABlock*sizeof(qreal);
        if (firstTime) {
             statevec_calcInnerProductKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                 getRealComp,
                 bra.deviceStateVec.real, bra.deviceStateVec.imag, 
                 ket.deviceStateVec.real, ket.deviceStateVec.imag, 
                 numValuesToReduce, 
                 bra.firstLevelReduction);
            firstTime = 0;
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    bra.firstLevelReduction, 
                    bra.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(bra.firstLevelReduction), &(bra.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    cudaMemcpy(&innerProdReal, bra.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    // compute imag component of inner product
    getRealComp = 0;
    numValuesToReduce = bra.numAmpsPerChunk;
    maxReducedPerLevel = REDUCE_SHARED_SIZE;
    firstTime = 1;
    while (numValuesToReduce > 1) {
        if (numValuesToReduce < maxReducedPerLevel) {
            valuesPerCUDABlock = numValuesToReduce;
            numCUDABlocks = 1;
        }
        else {
            valuesPerCUDABlock = maxReducedPerLevel;
            numCUDABlocks = ceil((qreal)numValuesToReduce/valuesPerCUDABlock);
        }
        sharedMemSize = valuesPerCUDABlock*sizeof(qreal);
        if (firstTime) {
             statevec_calcInnerProductKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                 getRealComp,
                 bra.deviceStateVec.real, bra.deviceStateVec.imag, 
                 ket.deviceStateVec.real, ket.deviceStateVec.imag, 
                 numValuesToReduce, 
                 bra.firstLevelReduction);
            firstTime = 0;
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    bra.firstLevelReduction, 
                    bra.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(bra.firstLevelReduction), &(bra.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    cudaMemcpy(&innerProdImag, bra.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    // return complex
    Complex innerProd;
    innerProd.real = innerProdReal;
    innerProd.imag = innerProdImag;
    return innerProd;
}

/** computes one term of (vec^*T) dens * vec */
__global__ void densmatr_calcFidelityKernel(Qureg dens, Qureg vec, long long int dim, qreal* reducedArray) {

    // figure out which density matrix row to consider
    long long int col;
    long long int row = blockIdx.x*blockDim.x + threadIdx.x;
    if (row >= dim) return;
    
    qreal* densReal = dens.deviceStateVec.real;
    qreal* densImag = dens.deviceStateVec.imag;
    qreal* vecReal  = vec.deviceStateVec.real;
    qreal* vecImag  = vec.deviceStateVec.imag;
    
    // compute the row-th element of the product dens*vec
    qreal prodReal = 0;
    qreal prodImag = 0;
    for (col=0LL; col < dim; col++) {
        qreal densElemReal = densReal[dim*col + row];
        qreal densElemImag = densImag[dim*col + row];
        
        prodReal += densElemReal*vecReal[col] - densElemImag*vecImag[col];
        prodImag += densElemReal*vecImag[col] + densElemImag*vecReal[col];
    }
    
    // multiply with row-th elem of (vec^*)
    qreal termReal = prodImag*vecImag[row] + prodReal*vecReal[row];
    
    // imag of every term should be zero, because each is a valid fidelity calc of an eigenstate
    //qreal termImag = prodImag*vecReal[row] - prodReal*vecImag[row];
    
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = termReal;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

// @TODO implement
qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState) {
    
    // we're summing the square of every term in the density matrix
    long long int densityDim = 1LL << qureg.numQubitsRepresented;
    long long int numValuesToReduce = densityDim;
    
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
            numCUDABlocks = ceil((qreal)numValuesToReduce/valuesPerCUDABlock);
        }
        // dictates size of reduction array
        sharedMemSize = valuesPerCUDABlock*sizeof(qreal);
        
        // spawn threads to sum the probs in each block
        // store the reduction in the pureState array
        if (firstTime) {
             densmatr_calcFidelityKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                 qureg, pureState, densityDim, pureState.firstLevelReduction);
            firstTime = 0;
            
        // sum the block probs
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    pureState.firstLevelReduction, 
                    pureState.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(pureState.firstLevelReduction), &(pureState.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal fidelity;
    cudaMemcpy(&fidelity, pureState.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    return fidelity;
}


__global__ void densmatr_calcPurityKernel(qreal* vecReal, qreal* vecImag, long long int numAmpsToSum, qreal *reducedArray) {
    
    // figure out which density matrix term this thread is assigned
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= numAmpsToSum) return;
    
    qreal term = vecReal[index]*vecReal[index] + vecImag[index]*vecImag[index];
    
    // array of each thread's collected probability, to be summed
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = term;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

/** Computes the trace of the density matrix squared */
qreal densmatr_calcPurity(Qureg qureg) {
    
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
            numCUDABlocks = ceil((qreal)numValuesToReduce/valuesPerCUDABlock);
        }
        // dictates size of reduction array
        sharedMemSize = valuesPerCUDABlock*sizeof(qreal);
        
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
    
    qreal traceDensSquared;
    cudaMemcpy(&traceDensSquared, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    return traceDensSquared;
}

__global__ void statevec_collapseToKnownProbOutcomeKernel(Qureg qureg, int measureQubit, int outcome, qreal totalProbability)
{
    // ----- sizes
    long long int sizeBlock,                                           // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                           // current block
         index;                                               // current index for first half block
    // ----- measured probability
    qreal   renorm;                                    // probability (returned) value
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
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;

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
void statevec_collapseToKnownProbOutcome(Qureg qureg, const int measureQubit, int outcome, qreal outcomeProb)
{        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_collapseToKnownProbOutcomeKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, measureQubit, outcome, outcomeProb);
}

/** Maps thread ID to a |..0..><..0..| state and then locates |0><1|, |1><0| and |1><1| */
__global__ void densmatr_collapseToKnownProbOutcomeKernel(
    qreal outcomeProb, qreal* vecReal, qreal *vecImag, long long int numBasesToVisit,
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
void densmatr_collapseToKnownProbOutcome(Qureg qureg, const int measureQubit, int outcome, qreal outcomeProb) {
    
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
    CUDABlocks = ceil(numBasesToVisit / (qreal) threadsPerCUDABlock);
    densmatr_collapseToKnownProbOutcomeKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        outcomeProb, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numBasesToVisit,
        part1, part2, part3, rowBit, colBit, desired, undesired);
}

__global__ void densmatr_addDensityMatrixKernel(Qureg combineQureg, qreal otherProb, Qureg otherQureg, long long int numAmpsToVisit) {
    
    long long int ampInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (ampInd >= numAmpsToVisit) return;
    
    combineQureg.deviceStateVec.real[ampInd] *= 1-otherProb;
    combineQureg.deviceStateVec.imag[ampInd] *= 1-otherProb;
  
    combineQureg.deviceStateVec.real[ampInd] += otherProb*otherQureg.deviceStateVec.real[ampInd];
    combineQureg.deviceStateVec.imag[ampInd] += otherProb*otherQureg.deviceStateVec.imag[ampInd];
}

void densmatr_addDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg) {
    
    long long int numAmpsToVisit = combineQureg.numAmpsPerChunk;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_addDensityMatrixKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        combineQureg, otherProb, otherQureg, numAmpsToVisit
    );
}

/** Called once for every 4 amplitudes in density matrix 
 * Works by establishing the |..0..><..0..| state (for its given index) then 
 * visiting |..1..><..0..| and |..0..><..1..|. Labels |part1 X pa><rt2 NOT(X) part3|
 * From the brain of Simon Benjamin
 */
__global__ void densmatr_oneQubitDephaseKernel(
    qreal fac, qreal* vecReal, qreal *vecImag, long long int numAmpsToVisit,
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


void densmatr_oneQubitDegradeOffDiagonal(Qureg qureg, const int targetQubit, qreal dephFac) {
    
    long long int numAmpsToVisit = qureg.numAmpsPerChunk/4;
    
    int rowQubit = targetQubit + qureg.numQubitsRepresented;
    long long int colBit = 1LL << targetQubit;
    long long int rowBit = 1LL << rowQubit;
    
    long long int part1 = colBit - 1;
    long long int part2 = (rowBit >> 1) - colBit;
    long long int part3 = numAmpsToVisit - (rowBit >> 1);
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_oneQubitDephaseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        dephFac, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, colBit, rowBit);
}

void densmatr_oneQubitDephase(Qureg qureg, const int targetQubit, qreal dephase) {
    
    if (dephase == 0)
        return;
    
    qreal dephFac = 1 - dephase;
    densmatr_oneQubitDegradeOffDiagonal(qureg, targetQubit, dephFac);
}

/** Called 12 times for every 16 amplitudes in density matrix 
 * Each sums from the |..0..0..><..0..0..| index to visit either
 * |..0..0..><..0..1..|,  |..0..0..><..1..0..|,  |..0..0..><..1..1..|,  |..0..1..><..0..0..|
 * etc and so on to |..1..1..><..1..0|. Labels |part1 0 part2 0 par><t3 0 part4 0 part5|.
 * From the brain of Simon Benjamin
 */
__global__ void densmatr_twoQubitDephaseKernel(
    qreal fac, qreal* vecReal, qreal *vecImag, long long int numBackgroundStates, long long int numAmpsToVisit,
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
void densmatr_twoQubitDephase(Qureg qureg, int qubit1, int qubit2, qreal dephase) {
    
    if (dephase == 0)
        return;
    
    // assumes qubit2 > qubit1
    
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
    qreal dephFac = 1 - dephase;
    
    // refers to states |a 0 b 0 c><d 0 e 0 f| (target qubits are fixed)
    long long int numBackgroundStates = qureg.numAmpsPerChunk/16;
    
    // 12 of these states experience dephasing
    long long int numAmpsToVisit = 12 * numBackgroundStates;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_twoQubitDephaseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        dephFac, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numBackgroundStates, numAmpsToVisit,
        part1, part2, part3, part4, part5, colBit1, rowBit1, colBit2, rowBit2);
}

/** Works like oneQubitDephase but modifies every other element, and elements are averaged in pairs */
__global__ void densmatr_oneQubitDepolariseKernel(
    qreal depolLevel, qreal* vecReal, qreal *vecImag, long long int numAmpsToVisit,
    long long int part1, long long int part2, long long int part3, 
    long long int bothBits)
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= numAmpsToVisit) return;
    
    long long int baseInd = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2);
    long long int targetInd = baseInd + bothBits;
    
    qreal realAvDepol = depolLevel * 0.5 * (vecReal[baseInd] + vecReal[targetInd]);
    qreal imagAvDepol = depolLevel * 0.5 * (vecImag[baseInd] + vecImag[targetInd]);
    
    vecReal[baseInd]   *= 1 - depolLevel;
    vecImag[baseInd]   *= 1 - depolLevel;
    vecReal[targetInd] *= 1 - depolLevel;
    vecImag[targetInd] *= 1 - depolLevel;
    
    vecReal[baseInd]   += realAvDepol;
    vecImag[baseInd]   += imagAvDepol;
    vecReal[targetInd] += realAvDepol;
    vecImag[targetInd] += imagAvDepol;
}

/** Works like oneQubitDephase but modifies every other element, and elements are averaged in pairs */
__global__ void densmatr_oneQubitDampingKernel(
    qreal damping, qreal* vecReal, qreal *vecImag, long long int numAmpsToVisit,
    long long int part1, long long int part2, long long int part3, 
    long long int bothBits)
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= numAmpsToVisit) return;
    
    long long int baseInd = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2);
    long long int targetInd = baseInd + bothBits;
    
    qreal realAvDepol = damping  * ( vecReal[targetInd]);
    qreal imagAvDepol = damping  * ( vecImag[targetInd]);
    
    vecReal[targetInd] *= 1 - damping;
    vecImag[targetInd] *= 1 - damping;
    
    vecReal[baseInd]   += realAvDepol;
    vecImag[baseInd]   += imagAvDepol;
}

void densmatr_oneQubitDepolarise(Qureg qureg, const int targetQubit, qreal depolLevel) {
    
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
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_oneQubitDepolariseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        depolLevel, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, bothBits);
}

void densmatr_oneQubitDamping(Qureg qureg, const int targetQubit, qreal damping) {
    
    if (damping == 0)
        return;
    
    qreal dephase = sqrt(1-damping);
    densmatr_oneQubitDegradeOffDiagonal(qureg, targetQubit, dephase);
    
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
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_oneQubitDampingKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        damping, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, bothBits);
}

/** Called once for every 16 amplitudes */
__global__ void densmatr_twoQubitDepolariseKernel(
    qreal depolLevel, qreal* vecReal, qreal *vecImag, long long int numAmpsToVisit,
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
    
    qreal realAvDepol = depolLevel * 0.25 * (
        vecReal[ind00] + vecReal[ind01] + vecReal[ind10] + vecReal[ind11]);
    qreal imagAvDepol = depolLevel * 0.25 * (
        vecImag[ind00] + vecImag[ind01] + vecImag[ind10] + vecImag[ind11]);
    
    qreal retain = 1 - depolLevel;
    vecReal[ind00] *= retain; vecImag[ind00] *= retain;
    vecReal[ind01] *= retain; vecImag[ind01] *= retain;
    vecReal[ind10] *= retain; vecImag[ind10] *= retain;
    vecReal[ind11] *= retain; vecImag[ind11] *= retain;

    vecReal[ind00] += realAvDepol; vecImag[ind00] += imagAvDepol;
    vecReal[ind01] += realAvDepol; vecImag[ind01] += imagAvDepol;
    vecReal[ind10] += realAvDepol; vecImag[ind10] += imagAvDepol;
    vecReal[ind11] += realAvDepol; vecImag[ind11] += imagAvDepol;
}

void densmatr_twoQubitDepolarise(Qureg qureg, int qubit1, int qubit2, qreal depolLevel) {
    
    if (depolLevel == 0)
        return;
    
    // assumes qubit2 > qubit1
    
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
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_twoQubitDepolariseKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        depolLevel, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, part4, part5, rowCol1, rowCol2);
}

void seedQuESTDefault(){
    // init MT random number generator with three keys -- time and pid
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process

    unsigned long int key[2];
    getQuESTDefaultSeedKey(key); 
    init_by_array(key, 2); 
}  


#ifdef __cplusplus
}
#endif
