// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * A custom-kernel implementation of the backend in ../QuEST_internal.h for a GPU environment.
 *
 * @author Ania Brown 
 * @author Tyson Jones
 */

# include "QuEST.h"
# include "QuEST_gpu_common.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "QuEST_internal.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

#ifdef USE_HIP
// Translate CUDA calls into HIP calls 
#include "cuda_to_hip.h"
#endif

# define REDUCE_SHARED_SIZE 512
# define DEBUG 0



/*
 * struct types for concisely passing unitaries to kernels
 */
 
 // hide these from doxygen
 /// \cond HIDDEN_SYMBOLS  
 
 typedef struct ArgMatrix2 {
     Complex r0c0, r0c1;
     Complex r1c0, r1c1;
 } ArgMatrix2;
 
 typedef struct ArgMatrix4
 {
     Complex r0c0, r0c1, r0c2, r0c3;
     Complex r1c0, r1c1, r1c2, r1c3;
     Complex r2c0, r2c1, r2c2, r2c3;
     Complex r3c0, r3c1, r3c2, r3c3;
 } ArgMatrix4;
 
ArgMatrix2 argifyMatrix2(ComplexMatrix2 m) {    
    ArgMatrix2 a;
    a.r0c0.real=m.real[0][0]; a.r0c0.imag=m.imag[0][0];
    a.r0c1.real=m.real[0][1]; a.r0c1.imag=m.imag[0][1];
    a.r1c0.real=m.real[1][0]; a.r1c0.imag=m.imag[1][0];
    a.r1c1.real=m.real[1][1]; a.r1c1.imag=m.imag[1][1];
    return a;
 }

ArgMatrix4 argifyMatrix4(ComplexMatrix4 m) {     
    ArgMatrix4 a;
    a.r0c0.real=m.real[0][0]; a.r0c0.imag=m.imag[0][0];
    a.r0c1.real=m.real[0][1]; a.r0c1.imag=m.imag[0][1];
    a.r0c2.real=m.real[0][2]; a.r0c2.imag=m.imag[0][2];
    a.r0c3.real=m.real[0][3]; a.r0c3.imag=m.imag[0][3];
    a.r1c0.real=m.real[1][0]; a.r1c0.imag=m.imag[1][0];
    a.r1c1.real=m.real[1][1]; a.r1c1.imag=m.imag[1][1];
    a.r1c2.real=m.real[1][2]; a.r1c2.imag=m.imag[1][2];
    a.r1c3.real=m.real[1][3]; a.r1c3.imag=m.imag[1][3];
    a.r2c0.real=m.real[2][0]; a.r2c0.imag=m.imag[2][0];
    a.r2c1.real=m.real[2][1]; a.r2c1.imag=m.imag[2][1];
    a.r2c2.real=m.real[2][2]; a.r2c2.imag=m.imag[2][2];
    a.r2c3.real=m.real[2][3]; a.r2c3.imag=m.imag[2][3];
    a.r3c0.real=m.real[3][0]; a.r3c0.imag=m.imag[3][0];
    a.r3c1.real=m.real[3][1]; a.r3c1.imag=m.imag[3][1];
    a.r3c2.real=m.real[3][2]; a.r3c2.imag=m.imag[3][2];
    a.r3c3.real=m.real[3][3]; a.r3c3.imag=m.imag[3][3];
    return a;
 }
 
 /// \endcond



/*
 * in-kernel bit twiddling functions
 */

__forceinline__ __device__ int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber) {
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

__forceinline__ __device__ int getBitMaskParity(long long int mask) {
    int parity = 0;
    while (mask) {
        parity = !parity;
        mask = mask & (mask-1);
    }
    return parity;
}

__forceinline__ __device__ long long int flipBit(const long long int number, const int bitInd) {
    return (number ^ (1LL << bitInd));
}

__forceinline__ __device__ long long int insertZeroBit(const long long int number, const int index) {
    long long int left, right;
    left = (number >> index) << index;
    right = number - left;
    return (left << 1) ^ right;
}

__forceinline__ __device__ long long int insertTwoZeroBits(const long long int number, const int bit1, const int bit2) {
    int small = (bit1 < bit2)? bit1 : bit2;
    int big = (bit1 < bit2)? bit2 : bit1;
    return insertZeroBit(insertZeroBit(number, small), big);
}

__forceinline__ __device__ long long int insertZeroBits(long long int number, int* inds, const int numInds) {
    /* inserted bit inds must strictly increase, so that their final indices are correct.
     * in-lieu of sorting (avoided since no C++ variable-size arrays, and since we're already 
     * memory bottle-necked so overhead eats this slowdown), we find the next-smallest index each 
     * at each insert. recall every element of inds (a positive or zero number) is unique.
     * This function won't appear in the CPU code, which can use C99 variable-size arrays and 
     * ought to make a sorted array before threading
     */
     int curMin = inds[0];
     int prevMin = -1;
     for (int n=0; n < numInds; n++) {
         
         // find next min
         for (int t=0; t < numInds; t++)
            if (inds[t]>prevMin && inds[t]<curMin)
                curMin = inds[t];
        
        number = insertZeroBit(number, curMin);
        
        // set curMin to an arbitrary non-visited elem
        prevMin = curMin;
        for (int t=0; t < numInds; t++)
            if (inds[t] > curMin) {
                curMin = inds[t];
                break;
            }
     }
     return number;
}



/*
 * state vector and density matrix operations 
 */

#ifdef __cplusplus
extern "C" {
#endif



QuESTEnv createQuESTEnv(void) {
    validateGPUExists(GPUExists(), __func__);
    
    QuESTEnv env;
    env.rank=0;
    env.numRanks=1;
    
    env.seeds = NULL;
    env.numSeeds = 0;
    seedQuESTDefault(&env);

    return env;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);
}

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
        numAmps * sizeof(*(qureg.deviceStateVec.imag)), 
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

    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsPerChunk = numAmpsPerRank;
    qureg->numAmpsTotal = numAmps;
    qureg->chunkId = env.rank;
    qureg->numChunks = env.numRanks;
    qureg->isDensityMatrix = 0;

    // check cpu memory allocation was successful
    validateQuregAllocation(qureg, env, __func__);

    // allocate GPU memory
    cudaMalloc(&(qureg->deviceStateVec.real), qureg->numAmpsPerChunk*sizeof(*(qureg->deviceStateVec.real)));
    cudaMalloc(&(qureg->deviceStateVec.imag), qureg->numAmpsPerChunk*sizeof(*(qureg->deviceStateVec.imag)));
    cudaMalloc(&(qureg->firstLevelReduction), ceil(qureg->numAmpsPerChunk/(qreal)REDUCE_SHARED_SIZE)*sizeof(qreal));
    cudaMalloc(&(qureg->secondLevelReduction), ceil(qureg->numAmpsPerChunk/(qreal)(REDUCE_SHARED_SIZE*REDUCE_SHARED_SIZE))*
            sizeof(qreal));

    // check gpu memory allocation was successful
    validateQuregGPUAllocation(qureg, env, __func__);
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
    cudaFree(qureg.firstLevelReduction);
    cudaFree(qureg.secondLevelReduction);
}

void copyStateToGPU(Qureg qureg)
{
    if (DEBUG) printf("Copying data to GPU\n");
    cudaMemcpy(qureg.deviceStateVec.real, qureg.stateVec.real, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyHostToDevice);
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

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    if (DEBUG) printf("Copying data to GPU\n");
    cudaMemcpy(&(qureg.deviceStateVec.real[startInd]), &(qureg.stateVec.real[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyHostToDevice);
    cudaMemcpy(&(qureg.deviceStateVec.imag[startInd]), &(qureg.stateVec.imag[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyHostToDevice);
    if (DEBUG) printf("Finished copying data to GPU\n");
}

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    cudaDeviceSynchronize();
    if (DEBUG) printf("Copying data from GPU\n");
    cudaMemcpy(&(qureg.stateVec.real[startInd]), &(qureg.deviceStateVec.real[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.real)), cudaMemcpyDeviceToHost);
    cudaMemcpy(&(qureg.stateVec.imag[startInd]), &(qureg.deviceStateVec.imag[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.imag)), cudaMemcpyDeviceToHost);
    if (DEBUG) printf("Finished copying data from GPU\n");
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

__global__ void statevec_initBlankStateKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag){
    long long int index;

    // initialise the statevector to be all-zeros
    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;
    stateVecReal[index] = 0.0;
    stateVecImag[index] = 0.0;
}

void statevec_initBlankState(Qureg qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initBlankStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk, 
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
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

__global__ void statevec_initDebugStateKernel(long long int stateVecSize, qreal *stateVecReal, qreal *stateVecImag){
    long long int index;

    index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;

    stateVecReal[index] = (index*2.0)/10.0;
    stateVecImag[index] = (index*2.0+1.0)/10.0;
}

void statevec_initDebugState(Qureg qureg)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_initDebugStateKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg.numAmpsPerChunk,
        qureg.deviceStateVec.real, 
        qureg.deviceStateVec.imag);
}

__global__ void statevec_compactUnitaryKernel (Qureg qureg, int rotQubit, Complex alpha, Complex beta){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;

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

void statevec_compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_compactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, alpha, beta);
}

__global__ void statevec_controlledCompactUnitaryKernel (Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;
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

void statevec_controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledCompactUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, alpha, beta);
}

__global__ void statevec_unitaryKernel(Qureg qureg, int targetQubit, ArgMatrix2 u){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;

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

void statevec_unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_unitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, argifyMatrix2(u));
}

__global__ void statevec_multiControlledMultiQubitUnitaryKernel(
    Qureg qureg, long long int ctrlMask, int* targs, int numTargs, 
    qreal* uRe, qreal* uIm, long long int* ampInds, qreal* reAmps, qreal* imAmps, long long int numTargAmps)
{
    
    // decide the amplitudes this thread will modify
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;                        
    long long int numTasks = qureg.numAmpsPerChunk >> numTargs; // kernel called on every 1 in 2^numTargs amplitudes
    if (thisTask>=numTasks) return;
    
    // find this task's start index (where all targs are 0)
    long long int ind00 = insertZeroBits(thisTask, targs, numTargs);
    
    // this task only modifies amplitudes if control qubits are 1 for this state
    if (ctrlMask && (ctrlMask&ind00) != ctrlMask)
        return;
        
    qreal *reVec = qureg.deviceStateVec.real;
    qreal *imVec = qureg.deviceStateVec.imag;
    
    /*
    each thread needs:
        long long int ampInds[numAmps];
        qreal reAmps[numAmps];
        qreal imAmps[numAmps];
    but instead has access to shared arrays, with below stride and offset
    */
    size_t stride = gridDim.x*blockDim.x;
    size_t offset = blockIdx.x*blockDim.x + threadIdx.x;
    
    // determine the indices and record values of target amps
    long long int ind;
    for (int i=0; i < numTargAmps; i++) {
        
        // get global index of current target qubit assignment
        ind = ind00;
        for (int t=0; t < numTargs; t++)
            if (extractBit(t, i))
                ind = flipBit(ind, targs[t]);
        
        ampInds[i*stride+offset] = ind;
        reAmps [i*stride+offset] = reVec[ind];
        imAmps [i*stride+offset] = imVec[ind];
    }
    
    // update the amplitudes
    for (int r=0; r < numTargAmps; r++) {
        ind = ampInds[r*stride+offset];
        reVec[ind] = 0;
        imVec[ind] = 0;
        for (int c=0; c < numTargAmps; c++) {
            qreal uReElem = uRe[c + r*numTargAmps];
            qreal uImElem = uIm[c + r*numTargAmps];
            reVec[ind] += reAmps[c*stride+offset]*uReElem - imAmps[c*stride+offset]*uImElem;
            imVec[ind] += reAmps[c*stride+offset]*uImElem + imAmps[c*stride+offset]*uReElem;
        }
    }
}

void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u)
{
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>numTargs)/threadsPerCUDABlock);
    
    // allocate device space for global {targs} (length: numTargs) and populate
    int *d_targs;
    size_t targMemSize = numTargs * sizeof *d_targs;
    cudaMalloc(&d_targs, targMemSize);
    cudaMemcpy(d_targs, targs, targMemSize, cudaMemcpyHostToDevice);
    
    // flatten out the u.real and u.imag lists
    int uNumRows = (1 << u.numQubits);
    qreal* uReFlat = (qreal*) malloc(uNumRows*uNumRows * sizeof *uReFlat);
    qreal* uImFlat = (qreal*) malloc(uNumRows*uNumRows * sizeof *uImFlat);
    long long int i = 0;
    for (int r=0; r < uNumRows; r++)
        for (int c=0; c < uNumRows; c++) {
            uReFlat[i] = u.real[r][c];
            uImFlat[i] = u.imag[r][c];
            i++;
        }
    
    // allocate device space for global u.real and u.imag (flatten by concatenating rows) and populate
    qreal* d_uRe;
    qreal* d_uIm;
    size_t uMemSize = uNumRows*uNumRows * sizeof *d_uRe; // size of each of d_uRe and d_uIm
    cudaMalloc(&d_uRe, uMemSize);
    cudaMalloc(&d_uIm, uMemSize);
    cudaMemcpy(d_uRe, uReFlat, uMemSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d_uIm, uImFlat, uMemSize, cudaMemcpyHostToDevice);
    
    // allocate device Wspace for thread-local {ampInds}, {reAmps}, {imAmps} (length: 1<<numTargs)
    long long int *d_ampInds;
    qreal *d_reAmps;
    qreal *d_imAmps;
    size_t gridSize = (size_t) threadsPerCUDABlock * CUDABlocks;
    int numTargAmps = uNumRows;
    cudaMalloc(&d_ampInds, numTargAmps*gridSize * sizeof *d_ampInds);
    cudaMalloc(&d_reAmps,  numTargAmps*gridSize * sizeof *d_reAmps);
    cudaMalloc(&d_imAmps,  numTargAmps*gridSize * sizeof *d_imAmps);
    
    // call kernel
    statevec_multiControlledMultiQubitUnitaryKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, ctrlMask, d_targs, numTargs, d_uRe, d_uIm, d_ampInds, d_reAmps, d_imAmps, numTargAmps);
        
    // free kernel memory
    free(uReFlat);
    free(uImFlat);
    cudaFree(d_targs);
    cudaFree(d_uRe);
    cudaFree(d_uIm);
    cudaFree(d_ampInds);
    cudaFree(d_reAmps);
    cudaFree(d_imAmps);
}

__global__ void statevec_multiControlledTwoQubitUnitaryKernel(Qureg qureg, long long int ctrlMask, int q1, int q2, ArgMatrix4 u){
    
    // decide the 4 amplitudes this thread will modify
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;                        
    long long int numTasks = qureg.numAmpsPerChunk >> 2; // kernel called on every 1 in 4 amplitudes
    if (thisTask>=numTasks) return;
    
    qreal *reVec = qureg.deviceStateVec.real;
    qreal *imVec = qureg.deviceStateVec.imag;
    
    // find indices of amplitudes to modify (treat q1 as the least significant bit)
    long long int ind00, ind01, ind10, ind11;
    ind00 = insertTwoZeroBits(thisTask, q1, q2);
    
    // modify only if control qubits are 1 for this state
    if (ctrlMask && (ctrlMask&ind00) != ctrlMask)
        return;
    
    ind01 = flipBit(ind00, q1);
    ind10 = flipBit(ind00, q2);
    ind11 = flipBit(ind01, q2);
    
    // extract statevec amplitudes 
    qreal re00, re01, re10, re11;
    qreal im00, im01, im10, im11;
    re00 = reVec[ind00]; im00 = imVec[ind00];
    re01 = reVec[ind01]; im01 = imVec[ind01];
    re10 = reVec[ind10]; im10 = imVec[ind10];
    re11 = reVec[ind11]; im11 = imVec[ind11];
    
    // apply u * {amp00, amp01, amp10, amp11}
    reVec[ind00] = 
        u.r0c0.real*re00 - u.r0c0.imag*im00 +
        u.r0c1.real*re01 - u.r0c1.imag*im01 +
        u.r0c2.real*re10 - u.r0c2.imag*im10 +
        u.r0c3.real*re11 - u.r0c3.imag*im11;
    imVec[ind00] =
        u.r0c0.imag*re00 + u.r0c0.real*im00 +
        u.r0c1.imag*re01 + u.r0c1.real*im01 +
        u.r0c2.imag*re10 + u.r0c2.real*im10 +
        u.r0c3.imag*re11 + u.r0c3.real*im11;
        
    reVec[ind01] = 
        u.r1c0.real*re00 - u.r1c0.imag*im00 +
        u.r1c1.real*re01 - u.r1c1.imag*im01 +
        u.r1c2.real*re10 - u.r1c2.imag*im10 +
        u.r1c3.real*re11 - u.r1c3.imag*im11;
    imVec[ind01] =
        u.r1c0.imag*re00 + u.r1c0.real*im00 +
        u.r1c1.imag*re01 + u.r1c1.real*im01 +
        u.r1c2.imag*re10 + u.r1c2.real*im10 +
        u.r1c3.imag*re11 + u.r1c3.real*im11;
        
    reVec[ind10] = 
        u.r2c0.real*re00 - u.r2c0.imag*im00 +
        u.r2c1.real*re01 - u.r2c1.imag*im01 +
        u.r2c2.real*re10 - u.r2c2.imag*im10 +
        u.r2c3.real*re11 - u.r2c3.imag*im11;
    imVec[ind10] =
        u.r2c0.imag*re00 + u.r2c0.real*im00 +
        u.r2c1.imag*re01 + u.r2c1.real*im01 +
        u.r2c2.imag*re10 + u.r2c2.real*im10 +
        u.r2c3.imag*re11 + u.r2c3.real*im11;    
        
    reVec[ind11] = 
        u.r3c0.real*re00 - u.r3c0.imag*im00 +
        u.r3c1.real*re01 - u.r3c1.imag*im01 +
        u.r3c2.real*re10 - u.r3c2.imag*im10 +
        u.r3c3.real*re11 - u.r3c3.imag*im11;
    imVec[ind11] =
        u.r3c0.imag*re00 + u.r3c0.real*im00 +
        u.r3c1.imag*re01 + u.r3c1.real*im01 +
        u.r3c2.imag*re10 + u.r3c2.real*im10 +
        u.r3c3.imag*re11 + u.r3c3.real*im11;    
}

void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int q1, int q2, ComplexMatrix4 u)
{
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>2)/threadsPerCUDABlock); // one kernel eval for every 4 amplitudes
    statevec_multiControlledTwoQubitUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, ctrlMask, q1, q2, argifyMatrix4(u));
}

__global__ void statevec_controlledUnitaryKernel(Qureg qureg, int controlQubit, int targetQubit, ArgMatrix2 u){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;

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

void statevec_controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_controlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, argifyMatrix2(u));
}

__global__ void statevec_multiControlledUnitaryKernel(
    Qureg qureg, 
    long long int ctrlQubitsMask, long long int ctrlFlipMask, 
    int targetQubit, ArgMatrix2 u
){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;


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

    if (ctrlQubitsMask == (ctrlQubitsMask & (indexUp ^ ctrlFlipMask))) {
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

void statevec_multiControlledUnitary(
    Qureg qureg, 
    long long int ctrlQubitsMask, long long int ctrlFlipMask, 
    int targetQubit, ComplexMatrix2 u
){
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_multiControlledUnitaryKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        qureg, ctrlQubitsMask, ctrlFlipMask, targetQubit, argifyMatrix2(u));
}

__global__ void statevec_pauliXKernel(Qureg qureg, int targetQubit){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;

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

void statevec_pauliX(Qureg qureg, int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_pauliXKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit);
}

__global__ void statevec_pauliYKernel(Qureg qureg, int targetQubit, int conjFac){

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

void statevec_pauliY(Qureg qureg, int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_pauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, 1);
}

void statevec_pauliYConj(Qureg qureg, int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_pauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, -1);
}

__global__ void statevec_controlledPauliYKernel(Qureg qureg, int controlQubit, int targetQubit, int conjFac)
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

void statevec_controlledPauliY(Qureg qureg, int controlQubit, int targetQubit)
{
    int conjFactor = 1;
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledPauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, conjFactor);
}

void statevec_controlledPauliYConj(Qureg qureg, int controlQubit, int targetQubit)
{
    int conjFactor = -1;
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledPauliYKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit, conjFactor);
}

__global__ void statevec_phaseShiftByTermKernel(Qureg qureg, int targetQubit, qreal cosAngle, qreal sinAngle) {

    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, indexUp,indexLo;

    qreal stateRealLo, stateImagLo;             
    long long int thisTask; 
    long long int numTasks = qureg.numAmpsPerChunk >> 1;

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

void statevec_phaseShiftByTerm(Qureg qureg, int targetQubit, Complex term)
{   
    qreal cosAngle = term.real;
    qreal sinAngle = term.imag;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_phaseShiftByTermKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit, cosAngle, sinAngle);
}

__global__ void statevec_controlledPhaseShiftKernel(Qureg qureg, int idQubit1, int idQubit2, qreal cosAngle, qreal sinAngle)
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

void statevec_controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle)
{
    qreal cosAngle = cos(angle);
    qreal sinAngle = sin(angle);
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
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

    long long int mask = getQubitBitMask(controlQubits, numControlQubits);
        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledPhaseShiftKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask, cosAngle, sinAngle);
}

__global__ void statevec_multiRotateZKernel(Qureg qureg, long long int mask, qreal cosAngle, qreal sinAngle) {
    
    long long int stateVecSize = qureg.numAmpsPerChunk;
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;
    
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;
    
    int fac = getBitMaskParity(mask & index)? -1 : 1;
    qreal stateReal = stateVecReal[index];
    qreal stateImag = stateVecImag[index];
    
    stateVecReal[index] = cosAngle*stateReal + fac * sinAngle*stateImag;
    stateVecImag[index] = - fac * sinAngle*stateReal + cosAngle*stateImag;  
}

void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle)
{   
    qreal cosAngle = cos(angle/2.0);
    qreal sinAngle = sin(angle/2.0);
        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiRotateZKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask, cosAngle, sinAngle);
}

__global__ void statevec_multiControlledMultiRotateZKernel(Qureg qureg, long long int ctrlMask, long long int targMask, qreal cosAngle, qreal sinAngle) {
    
    long long int stateVecSize = qureg.numAmpsPerChunk;
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=stateVecSize) return;
    
    // amplitudes corresponding to control qubits not all-in-one are unmodified
    if (ctrlMask && ((ctrlMask & index) != ctrlMask))
        return;
    
    qreal *stateVecReal = qureg.deviceStateVec.real;
    qreal *stateVecImag = qureg.deviceStateVec.imag;
    
    // avoid warp divergence, setting fac = +- 1
    int fac = 1-2*getBitMaskParity(targMask & index);
    qreal stateReal = stateVecReal[index];
    qreal stateImag = stateVecImag[index];
    
    stateVecReal[index] = cosAngle*stateReal + fac * sinAngle*stateImag;
    stateVecImag[index] = - fac * sinAngle*stateReal + cosAngle*stateImag;  
}

void statevec_multiControlledMultiRotateZ(Qureg qureg, long long int ctrlMask, long long int targMask, qreal angle)
{   
    qreal cosAngle = cos(angle/2.0);
    qreal sinAngle = sin(angle/2.0);
        
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledMultiRotateZKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, ctrlMask, targMask, cosAngle, sinAngle);
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

__global__ void statevec_controlledPhaseFlipKernel(Qureg qureg, int idQubit1, int idQubit2)
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

void statevec_controlledPhaseFlip(Qureg qureg, int idQubit1, int idQubit2)
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
    long long int mask = getQubitBitMask(controlQubits, numControlQubits);
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_multiControlledPhaseFlipKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, mask);
}

__global__ void statevec_swapQubitAmpsKernel(Qureg qureg, int qb1, int qb2) {

    qreal *reVec = qureg.deviceStateVec.real;
    qreal *imVec = qureg.deviceStateVec.imag;
    
    long long int numTasks = qureg.numAmpsPerChunk >> 2; // each iteration updates 2 amps and skips 2 amps
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask>=numTasks) return;
    
    long long int ind00, ind01, ind10;
    qreal re01, re10, im01, im10;
  
    // determine ind00 of |..0..0..>, |..0..1..> and |..1..0..>
    ind00 = insertTwoZeroBits(thisTask, qb1, qb2);
    ind01 = flipBit(ind00, qb1);
    ind10 = flipBit(ind00, qb2);

    // extract statevec amplitudes 
    re01 = reVec[ind01]; im01 = imVec[ind01];
    re10 = reVec[ind10]; im10 = imVec[ind10];

    // swap 01 and 10 amps
    reVec[ind01] = re10; reVec[ind10] = re01;
    imVec[ind01] = im10; imVec[ind10] = im01;
}

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>2)/threadsPerCUDABlock);
    statevec_swapQubitAmpsKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, qb1, qb2);
}

__global__ void statevec_hadamardKernel (Qureg qureg, int targetQubit){
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
    long long int numTasks=qureg.numAmpsPerChunk>>1;

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

void statevec_hadamard(Qureg qureg, int targetQubit) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk>>1)/threadsPerCUDABlock);
    statevec_hadamardKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, targetQubit);
}

__global__ void statevec_controlledNotKernel(Qureg qureg, int controlQubit, int targetQubit)
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

void statevec_controlledNot(Qureg qureg, int controlQubit, int targetQubit)
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_controlledNotKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, controlQubit, targetQubit);
}

__global__ void statevec_multiControlledMultiQubitNotKernel(Qureg qureg, int ctrlMask, int targMask) {
    
    qreal* stateRe = qureg.deviceStateVec.real;
    qreal* stateIm = qureg.deviceStateVec.imag;
    
    // althouugh each thread swaps/updates two amplitudes, we still invoke one thread per amp
    long long int ampInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (ampInd >= qureg.numAmpsPerChunk)
        return;

    // modify amplitudes only if control qubits are 1 for this state
    if (ctrlMask && ((ctrlMask & ampInd) != ctrlMask))
        return;
    
    long long int mateInd = ampInd ^ targMask;
    
    // if the mate is lower index, another thread is handling it
    if (mateInd < ampInd)
        return;
        
    /* it may seem wasteful to spawn more threads than are needed, and abort 
     * half of them due to the amp pairing above (and potentially abort
     * an exponential number due to ctrlMask). however, since we are moving 
     * global memory directly in a potentially non-contiguous fashoin, this 
     * method is likely to be memory bandwidth bottlenecked anyway 
     */
    
    qreal mateRe = stateRe[mateInd];
    qreal mateIm = stateIm[mateInd];
    
    // swap amp with mate
    stateRe[mateInd] = stateRe[ampInd];
    stateIm[mateInd] = stateIm[ampInd];
    stateRe[ampInd] = mateRe;
    stateIm[ampInd] = mateIm;
}

void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask) {
    
    int numThreadsPerBlock = 128;
    int numBlocks = ceil(qureg.numAmpsPerChunk / (qreal) numThreadsPerBlock);
    statevec_multiControlledMultiQubitNotKernel<<<numBlocks, numThreadsPerBlock>>>(qureg, ctrlMask, targMask);
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
    Qureg qureg, int measureQubit, qreal *reducedArray
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
        Qureg qureg, int measureQubit, qreal *reducedArray
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

qreal densmatr_findProbabilityOfZero(Qureg qureg, int measureQubit)
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

qreal statevec_findProbabilityOfZero(Qureg qureg, int measureQubit)
{
    qreal stateProb=0;
    
    // 1-qubit edge-case breaks below loop logic
    if (qureg.numQubitsInStateVec == 1) {
        qreal amp;
        cudaMemcpy(&amp, qureg.deviceStateVec.real, sizeof(qreal), cudaMemcpyDeviceToHost);
        stateProb += amp*amp;
        cudaMemcpy(&amp, qureg.deviceStateVec.imag, sizeof(qreal), cudaMemcpyDeviceToHost);
        stateProb += amp*amp;
        return stateProb;
    }
    
    long long int numValuesToReduce = qureg.numAmpsPerChunk>>1;
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
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

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    qreal outcomeProb = statevec_findProbabilityOfZero(qureg, measureQubit);
    if (outcome==1)
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    qreal outcomeProb = densmatr_findProbabilityOfZero(qureg, measureQubit);
    if (outcome==1) 
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

/** computes Tr(conjTrans(a) b) = sum of (a_ij^* b_ij), which is a real number */
__global__ void densmatr_calcInnerProductKernel(
    Qureg a, Qureg b, long long int numTermsToSum, qreal* reducedArray
) {    
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= numTermsToSum) return;
    
    // Re{ conj(a) b } = Re{ (aRe - i aIm)(bRe + i bIm) } = aRe bRe + aIm bIm
    qreal prod = (
          a.deviceStateVec.real[index]*b.deviceStateVec.real[index] 
        + a.deviceStateVec.imag[index]*b.deviceStateVec.imag[index]);
    
    // array of each thread's collected sum term, to be summed
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = prod;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

qreal densmatr_calcInnerProduct(Qureg a, Qureg b) {
    
    // we're summing the square of every term in the density matrix
    long long int numValuesToReduce = a.numAmpsTotal;
    
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
        
        // spawn threads to sum the terms in each block
        // arbitrarily store the reduction in the b qureg's array
        if (firstTime) {
             densmatr_calcInnerProductKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                 a, b, a.numAmpsTotal, b.firstLevelReduction);
            firstTime = 0;
        }    
        // sum the block terms
        else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    b.firstLevelReduction, 
                    b.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(b.firstLevelReduction), &(b.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal innerprod;
    cudaMemcpy(&innerprod, b.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    return innerprod;
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
    
    // array of each thread's collected sum term, to be summed
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
 * @todo could even do the kernel twice, storing real in bra.reduc and imag in ket.reduc?
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

__global__ void densmatr_calcHilbertSchmidtDistanceSquaredKernel(
    qreal* aRe, qreal* aIm, qreal* bRe, qreal* bIm, 
    long long int numAmpsToSum, qreal *reducedArray
) {
    // figure out which density matrix term this thread is assigned
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= numAmpsToSum) return;
    
    // compute this thread's sum term
    qreal difRe = aRe[index] - bRe[index];
    qreal difIm = aIm[index] - bIm[index];
    qreal term = difRe*difRe + difIm*difIm;
    
    // array of each thread's collected term, to be summed
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = term;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

/* computes sqrt(Tr( (a-b) conjTrans(a-b) ) = sqrt( sum of abs vals of (a-b)) */
qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b) {
    
    // we're summing the square of every term in (a-b)
    long long int numValuesToReduce = a.numAmpsPerChunk;
    
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
        
        // spawn threads to sum the probs in each block (store reduction temp values in a's reduction array)
        if (firstTime) {
             densmatr_calcHilbertSchmidtDistanceSquaredKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                 a.deviceStateVec.real, a.deviceStateVec.imag, 
                 b.deviceStateVec.real, b.deviceStateVec.imag, 
                 numValuesToReduce, a.firstLevelReduction);
            firstTime = 0;
            
        // sum the block probs
        } else {
            cudaDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    a.firstLevelReduction, 
                    a.secondLevelReduction, valuesPerCUDABlock); 
            cudaDeviceSynchronize();    
            swapDouble(&(a.firstLevelReduction), &(a.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal trace;
    cudaMemcpy(&trace, a.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    qreal sqrtTrace = sqrt(trace);
    return sqrtTrace;
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
void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
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
void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb) {
    
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

__global__ void densmatr_mixDensityMatrixKernel(Qureg combineQureg, qreal otherProb, Qureg otherQureg, long long int numAmpsToVisit) {
    
    long long int ampInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (ampInd >= numAmpsToVisit) return;
    
    combineQureg.deviceStateVec.real[ampInd] *= 1-otherProb;
    combineQureg.deviceStateVec.imag[ampInd] *= 1-otherProb;
  
    combineQureg.deviceStateVec.real[ampInd] += otherProb*otherQureg.deviceStateVec.real[ampInd];
    combineQureg.deviceStateVec.imag[ampInd] += otherProb*otherQureg.deviceStateVec.imag[ampInd];
}

void densmatr_mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg) {
    
    long long int numAmpsToVisit = combineQureg.numAmpsPerChunk;
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    densmatr_mixDensityMatrixKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        combineQureg, otherProb, otherQureg, numAmpsToVisit
    );
}

/** Called once for every 4 amplitudes in density matrix 
 * Works by establishing the |..0..><..0..| state (for its given index) then 
 * visiting |..1..><..0..| and |..0..><..1..|. Labels |part1 X pa><rt2 NOT(X) part3|
 * From the brain of Simon Benjamin
 */
__global__ void densmatr_mixDephasingKernel(
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


void densmatr_oneQubitDegradeOffDiagonal(Qureg qureg, int targetQubit, qreal dephFac) {
    
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
    densmatr_mixDephasingKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        dephFac, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, colBit, rowBit);
}

void densmatr_mixDephasing(Qureg qureg, int targetQubit, qreal dephase) {
    
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
__global__ void densmatr_mixTwoQubitDephasingKernel(
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
void densmatr_mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal dephase) {
    
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
    densmatr_mixTwoQubitDephasingKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        dephFac, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numBackgroundStates, numAmpsToVisit,
        part1, part2, part3, part4, part5, colBit1, rowBit1, colBit2, rowBit2);
}

/** Works like mixDephasing but modifies every other element, and elements are averaged in pairs */
__global__ void densmatr_mixDepolarisingKernel(
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

/** Works like mixDephasing but modifies every other element, and elements are averaged in pairs */
__global__ void densmatr_mixDampingKernel(
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

void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depolLevel) {
    
    if (depolLevel == 0)
        return;
    
    densmatr_mixDephasing(qureg, targetQubit, depolLevel);
    
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
    densmatr_mixDepolarisingKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        depolLevel, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, bothBits);
}

void densmatr_mixDamping(Qureg qureg, int targetQubit, qreal damping) {
    
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
    densmatr_mixDampingKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        damping, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, bothBits);
}

__global__ void statevec_setWeightedQuregKernel(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out) {

    long long int ampInd = blockIdx.x*blockDim.x + threadIdx.x;
    long long int numAmpsToVisit = qureg1.numAmpsPerChunk;
    if (ampInd >= numAmpsToVisit) return;

    qreal *vecRe1 = qureg1.deviceStateVec.real;
    qreal *vecIm1 = qureg1.deviceStateVec.imag;
    qreal *vecRe2 = qureg2.deviceStateVec.real;
    qreal *vecIm2 = qureg2.deviceStateVec.imag;
    qreal *vecReOut = out.deviceStateVec.real;
    qreal *vecImOut = out.deviceStateVec.imag;

    qreal facRe1 = fac1.real; 
    qreal facIm1 = fac1.imag;
    qreal facRe2 = fac2.real;
    qreal facIm2 = fac2.imag;
    qreal facReOut = facOut.real;
    qreal facImOut = facOut.imag;

    qreal re1,im1, re2,im2, reOut,imOut;
    long long int index = ampInd;

    re1 = vecRe1[index]; im1 = vecIm1[index];
    re2 = vecRe2[index]; im2 = vecIm2[index];
    reOut = vecReOut[index];
    imOut = vecImOut[index];

    vecReOut[index] = (facReOut*reOut - facImOut*imOut) + (facRe1*re1 - facIm1*im1) + (facRe2*re2 - facIm2*im2);
    vecImOut[index] = (facReOut*imOut + facImOut*reOut) + (facRe1*im1 + facIm1*re1) + (facRe2*im2 + facIm2*re2);
}

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out) {

    long long int numAmpsToVisit = qureg1.numAmpsPerChunk;

    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil(numAmpsToVisit / (qreal) threadsPerCUDABlock);
    statevec_setWeightedQuregKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        fac1, qureg1, fac2, qureg2, facOut, out
    );
}

__global__ void statevec_applyDiagonalOpKernel(Qureg qureg, DiagonalOp op) {

    // each thread modifies one value; a wasteful and inefficient strategy
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;

    qreal* stateRe = qureg.deviceStateVec.real;
    qreal* stateIm = qureg.deviceStateVec.imag;
    qreal* opRe = op.deviceOperator.real;
    qreal* opIm = op.deviceOperator.imag;

    qreal a = stateRe[thisTask];
    qreal b = stateIm[thisTask];
    qreal c = opRe[thisTask];
    qreal d = opIm[thisTask];

    // (a + b i)(c + d i) = (a c - b d) + i (a d + b c)
    stateRe[thisTask] = a*c - b*d;
    stateIm[thisTask] = a*d + b*c;
}

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) 
{
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    statevec_applyDiagonalOpKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, op);
}

__global__ void densmatr_applyDiagonalOpKernel(Qureg qureg, DiagonalOp op) {

    // each thread modifies one value; a wasteful and inefficient strategy
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int thisTask = blockIdx.x*blockDim.x + threadIdx.x;
    if (thisTask >= numTasks) return;

    qreal* stateRe = qureg.deviceStateVec.real;
    qreal* stateIm = qureg.deviceStateVec.imag;
    qreal* opRe = op.deviceOperator.real;
    qreal* opIm = op.deviceOperator.imag;

    int opDim = (1 << op.numQubits);
    qreal a = stateRe[thisTask];
    qreal b = stateIm[thisTask];
    qreal c = opRe[thisTask % opDim];
    qreal d = opIm[thisTask % opDim];

    // (a + b i)(c + d i) = (a c - b d) + i (a d + b c)
    stateRe[thisTask] = a*c - b*d;
    stateIm[thisTask] = a*d + b*c;
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op) {
    
    int threadsPerCUDABlock, CUDABlocks;
    threadsPerCUDABlock = 128;
    CUDABlocks = ceil((qreal)(qureg.numAmpsPerChunk)/threadsPerCUDABlock);
    densmatr_applyDiagonalOpKernel<<<CUDABlocks, threadsPerCUDABlock>>>(qureg, op);
}

__global__ void statevec_applySubDiagonalOpKernel(Qureg qureg, int* targets, int numTargets, qreal* opReals, qreal* opImags, int conjFac) {
    
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=qureg.numAmpsPerChunk) return;
    
    long long int v = 0;
    for (int t=0; t<numTargets; t++)
        v |= extractBit(targets[t], index) << t;
        
    qreal elemRe = opReals[v];
    qreal elemIm = opImags[v] * conjFac;
    
    qreal* stateRe = qureg.deviceStateVec.real;
    qreal* stateIm = qureg.deviceStateVec.imag;
    
    qreal ampRe = stateRe[index];
    qreal ampIm = stateIm[index];
    
    // (a + b i)(c + d i) = (a c - b d) + i (a d + b c)
    stateRe[index] = ampRe*elemRe - ampIm*elemIm;
    stateIm[index] = ampRe*elemIm + ampIm*elemRe;
}

void statevec_applySubDiagonalOp(Qureg qureg, int* targets, SubDiagonalOp op, int conj) {
    
    // copy targets to GPU memory
    int* d_targets;
    int numTargets = op.numQubits;
    size_t memTargets = numTargets * sizeof *d_targets;
    cudaMalloc(&d_targets, memTargets);
    cudaMemcpy(d_targets, targets, memTargets, cudaMemcpyHostToDevice);
    
    // copy op to GPU memory
    qreal* d_opReal;
    qreal* d_opImag;
    size_t memOp = op.numElems * sizeof *d_opReal;
    cudaMalloc(&d_opReal, memOp);
    cudaMalloc(&d_opImag, memOp);
    cudaMemcpy(d_opReal, op.real, memOp, cudaMemcpyHostToDevice);
    cudaMemcpy(d_opImag, op.imag, memOp, cudaMemcpyHostToDevice);
    
    // determine factor of imaginary components
    int conjFac = 1;
    if (conj)
        conjFac = -1;
    
    // launch kernels
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil(qureg.numAmpsPerChunk / (qreal) threadsPerCUDABlock);
    statevec_applySubDiagonalOpKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_targets, numTargets, d_opReal, d_opImag, conjFac);
    
    // free temporary GPU memory
    cudaFree(d_targets);
    cudaFree(d_opReal);
    cudaFree(d_opImag);
}

/** computes either a real or imag term of |vec_i|^2 op_i */
__global__ void statevec_calcExpecDiagonalOpKernel(
    int getRealComp,
    qreal* vecReal, qreal* vecImag, qreal* opReal, qreal* opImag, 
    long long int numTermsToSum, qreal* reducedArray) 
{
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >= numTermsToSum) return;
    
    qreal vecAbs = vecReal[index]*vecReal[index] + vecImag[index]*vecImag[index];
    
    // choose whether to calculate the real or imaginary term of the expec term
    qreal expecVal;
    if (getRealComp)
        expecVal = vecAbs * opReal[index];
    else
        expecVal = vecAbs * opImag[index];
    
    // array of each thread's collected sum term, to be summed
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = expecVal;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {
    
    /* @TODO: remove all this reduction boilerplate from QuEST GPU 
     * (e.g. a func which accepts a pointer to do every-value reduction?)
     */

    qreal expecReal, expecImag;
    
    int getRealComp;
    long long int numValuesToReduce;
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    int maxReducedPerLevel;
    int firstTime;
    
    // compute real component of inner product
    getRealComp = 1;
    numValuesToReduce = qureg.numAmpsPerChunk;
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
            statevec_calcExpecDiagonalOpKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                getRealComp,
                qureg.deviceStateVec.real, qureg.deviceStateVec.imag, 
                op.deviceOperator.real, op.deviceOperator.imag, 
                numValuesToReduce, 
                qureg.firstLevelReduction);
            firstTime = 0;
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
    cudaMemcpy(&expecReal, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    // compute imag component of inner product
    getRealComp = 0;
    numValuesToReduce = qureg.numAmpsPerChunk;
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
            statevec_calcExpecDiagonalOpKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                getRealComp,
                qureg.deviceStateVec.real, qureg.deviceStateVec.imag, 
                op.deviceOperator.real, op.deviceOperator.imag, 
                numValuesToReduce, 
                qureg.firstLevelReduction);
            firstTime = 0;
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
    cudaMemcpy(&expecImag, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    // return complex
    Complex expecVal;
    expecVal.real = expecReal;
    expecVal.imag = expecImag;
    return expecVal;
}

__global__ void densmatr_calcExpecDiagonalOpKernel(
    int getRealComp,
    qreal* matReal, qreal* matImag, qreal* opReal, qreal* opImag, 
    int numQubits, long long int numTermsToSum, qreal* reducedArray) 
{
    /** if the thread represents a diagonal op, then it computes either a 
     *  real or imag term of matr_{ii} op_i. Otherwise, it writes a 0 to the 
     * reduction array
     */
    
    // index will identy one of the 2^Q diagonals to be summed
    long long int matInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (matInd >= numTermsToSum) return;
    
    long long int diagSpacing = (1LL << numQubits) + 1LL;
    int isDiag = ((matInd % diagSpacing) == 0);
    
    long long int opInd = matInd / diagSpacing;
    
    qreal val = 0;
    if (isDiag) {
        
        qreal matRe = matReal[matInd];
        qreal matIm = matImag[matInd];
        qreal opRe = opReal[opInd];
        qreal opIm = opImag[opInd];
        
        // (matRe + matIm i)(opRe + opIm i) = 
        //      (matRe opRe - matIm opIm) + i (matRe opIm + matIm opRe)
        if (getRealComp)
            val = matRe * opRe - matIm * opIm;
        else 
            val = matRe * opIm + matIm * opRe;
    }
    
    // array of each thread's collected sum term, to be summed
    extern __shared__ qreal tempReductionArray[];
    tempReductionArray[threadIdx.x] = val;
    __syncthreads();
    
    // every second thread reduces
    if (threadIdx.x<blockDim.x/2)
        reduceBlock(tempReductionArray, reducedArray, blockDim.x);
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {
    
    /* @TODO: remove all this reduction boilerplate from QuEST GPU 
     * (e.g. a func which accepts a pointer to do every-value reduction?)
     */

    qreal expecReal, expecImag;
    
    int getRealComp;
    long long int numValuesToReduce;
    int valuesPerCUDABlock, numCUDABlocks, sharedMemSize;
    int maxReducedPerLevel;
    int firstTime;
    
    // compute real component of inner product
    getRealComp = 1;
    numValuesToReduce = qureg.numAmpsPerChunk;
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
            densmatr_calcExpecDiagonalOpKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                getRealComp,
                qureg.deviceStateVec.real, qureg.deviceStateVec.imag, 
                op.deviceOperator.real, op.deviceOperator.imag, 
                op.numQubits, numValuesToReduce, 
                qureg.firstLevelReduction);
            firstTime = 0;
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
    cudaMemcpy(&expecReal, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    // compute imag component of inner product
    getRealComp = 0;
    numValuesToReduce = qureg.numAmpsPerChunk;
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
            densmatr_calcExpecDiagonalOpKernel<<<numCUDABlocks, valuesPerCUDABlock, sharedMemSize>>>(
                getRealComp,
                qureg.deviceStateVec.real, qureg.deviceStateVec.imag, 
                op.deviceOperator.real, op.deviceOperator.imag, 
                op.numQubits, numValuesToReduce, 
                qureg.firstLevelReduction);
            firstTime = 0;
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
    cudaMemcpy(&expecImag, qureg.firstLevelReduction, sizeof(qreal), cudaMemcpyDeviceToHost);
    
    // return complex
    Complex expecVal;
    expecVal.real = expecReal;
    expecVal.imag = expecImag;
    return expecVal;
}



#ifdef __cplusplus
}
#endif
