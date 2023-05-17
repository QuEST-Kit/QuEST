#include "hip/hip_runtime.h"
// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * An implementation of the backend in ../QuEST_internal.h for a GPU environment.
 *
 * @author Ania Brown 
 * @author Tyson Jones
 */

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "QuEST_internal.h"    // purely to resolve getQuESTDefaultSeedKey
# include "mt19937ar.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

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


void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps) {
    
    hipDeviceSynchronize();
    hipMemcpy(
        qureg.deviceStateVec.real + startInd, 
        reals,
        numAmps * sizeof(*(qureg.deviceStateVec.real)), 
        hipMemcpyHostToDevice);
    hipMemcpy(
        qureg.deviceStateVec.imag + startInd,
        imags,
        numAmps * sizeof(*(qureg.deviceStateVec.imag)), 
        hipMemcpyHostToDevice);
}


/** works for both statevectors and density matrices */
void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg) {
    
    // copy copyQureg's GPU statevec to targetQureg's GPU statevec
    hipDeviceSynchronize();
    hipMemcpy(
        targetQureg.deviceStateVec.real, 
        copyQureg.deviceStateVec.real, 
        targetQureg.numAmpsPerChunk*sizeof(*(targetQureg.deviceStateVec.real)), 
        hipMemcpyDeviceToDevice);
    hipMemcpy(
        targetQureg.deviceStateVec.imag, 
        copyQureg.deviceStateVec.imag, 
        targetQureg.numAmpsPerChunk*sizeof(*(targetQureg.deviceStateVec.imag)), 
        hipMemcpyDeviceToDevice);
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

    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsPerChunk = numAmpsPerRank;
    qureg->numAmpsTotal = numAmps;
    qureg->chunkId = env.rank;
    qureg->numChunks = env.numRanks;
    qureg->isDensityMatrix = 0;

    // check cpu memory allocation was successful
    validateQuregAllocation(qureg, env, __func__);

    // allocate GPU memory
    hipMalloc(&(qureg->deviceStateVec.real), qureg->numAmpsPerChunk*sizeof(*(qureg->deviceStateVec.real)));
    hipMalloc(&(qureg->deviceStateVec.imag), qureg->numAmpsPerChunk*sizeof(*(qureg->deviceStateVec.imag)));
    hipMalloc(&(qureg->firstLevelReduction), ceil(qureg->numAmpsPerChunk/(qreal)REDUCE_SHARED_SIZE)*sizeof(qreal));
    hipMalloc(&(qureg->secondLevelReduction), ceil(qureg->numAmpsPerChunk/(qreal)(REDUCE_SHARED_SIZE*REDUCE_SHARED_SIZE))*
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
    hipFree(qureg.deviceStateVec.real);
    hipFree(qureg.deviceStateVec.imag);
    hipFree(qureg.firstLevelReduction);
    hipFree(qureg.secondLevelReduction);
}

DiagonalOp agnostic_createDiagonalOp(int numQubits, QuESTEnv env) {

    DiagonalOp op;
    op.numQubits = numQubits;
    op.numElemsPerChunk = (1LL << numQubits) / env.numRanks;
    op.chunkId = env.rank;
    op.numChunks = env.numRanks;

    // allocate CPU memory (initialised to zero)
    op.real = (qreal*) calloc(op.numElemsPerChunk, sizeof(qreal));
    op.imag = (qreal*) calloc(op.numElemsPerChunk, sizeof(qreal));
    // @TODO no handling of rank>1 allocation (no distributed GPU)

    // check cpu memory allocation was successful
    validateDiagonalOpAllocation(&op, env, __func__);

    // allocate GPU memory
    size_t arrSize = op.numElemsPerChunk * sizeof(qreal);
    hipMalloc(&(op.deviceOperator.real), arrSize);
    hipMalloc(&(op.deviceOperator.imag), arrSize);

    // check gpu memory allocation was successful
    validateDiagonalOpGPUAllocation(&op, env, __func__);

    // initialise GPU memory to zero
    hipMemset(op.deviceOperator.real, 0, arrSize);
    hipMemset(op.deviceOperator.imag, 0, arrSize);

    return op;
}

void agnostic_destroyDiagonalOp(DiagonalOp op) {
    free(op.real);
    free(op.imag);
    hipFree(op.deviceOperator.real);
    hipFree(op.deviceOperator.imag);
}

void agnostic_syncDiagonalOp(DiagonalOp op) {
    hipDeviceSynchronize();
    size_t mem_elems = op.numElemsPerChunk * sizeof *op.real;
    hipMemcpy(op.deviceOperator.real, op.real, mem_elems, hipMemcpyHostToDevice);
    hipMemcpy(op.deviceOperator.imag, op.imag, mem_elems, hipMemcpyHostToDevice);
}

__global__ void agnostic_initDiagonalOpFromPauliHamilKernel(
    DiagonalOp op, enum pauliOpType* pauliCodes, qreal* termCoeffs, int numSumTerms
) {    
    // each thread processes one diagonal element
    long long int elemInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (elemInd >= op.numElemsPerChunk)
        return;

    qreal elem = 0;
    
    // elem is (+-) every coefficient, with sign determined by parity
    for (int t=0; t<numSumTerms; t++) {
        
        // determine the parity of the Z-targeted qubits in the element's corresponding state
        int isOddNumOnes = 0;
        for (int q=0; q<op.numQubits; q++)
            if (pauliCodes[q + t*op.numQubits] == PAULI_Z)
                if (extractBit(q, elemInd))
                    isOddNumOnes = !isOddNumOnes;
        
        // avoid warp divergence
        int sign = 1 - 2*isOddNumOnes; // (-1 if isOddNumOnes, else +1)
        elem += termCoeffs[t] * sign;
    }
    
    op.deviceOperator.real[elemInd] = elem;
    op.deviceOperator.imag[elemInd] = 0;
}

void agnostic_initDiagonalOpFromPauliHamil(DiagonalOp op, PauliHamil hamil) {
    
    // copy args intop GPU memory
    enum pauliOpType* d_pauliCodes;
    size_t mem_pauliCodes = hamil.numSumTerms * op.numQubits * sizeof *d_pauliCodes;
    hipMalloc(&d_pauliCodes, mem_pauliCodes);
    hipMemcpy(d_pauliCodes, hamil.pauliCodes, mem_pauliCodes, hipMemcpyHostToDevice);
    
    qreal* d_termCoeffs;
    size_t mem_termCoeffs = hamil.numSumTerms * sizeof *d_termCoeffs;
    hipMalloc(&d_termCoeffs, mem_termCoeffs);
    hipMemcpy(d_termCoeffs, hamil.termCoeffs, mem_termCoeffs, hipMemcpyHostToDevice);
    
    int numThreadsPerBlock = 128;
    int numBlocks = ceil(op.numElemsPerChunk / (qreal) numThreadsPerBlock);
    agnostic_initDiagonalOpFromPauliHamilKernel<<<numBlocks, numThreadsPerBlock>>>(
        op, d_pauliCodes, d_termCoeffs, hamil.numSumTerms);
    
    // copy populated operator into to RAM
    hipDeviceSynchronize();
    size_t mem_elems = op.numElemsPerChunk * sizeof *op.real;
    hipMemcpy(op.real, op.deviceOperator.real, mem_elems, hipMemcpyDeviceToHost);
    hipMemcpy(op.imag, op.deviceOperator.imag, mem_elems, hipMemcpyDeviceToHost);

    hipFree(d_pauliCodes);
    hipFree(d_termCoeffs);
}

int GPUExists(void){
    int deviceCount, device;
    int gpuDeviceCount = 0;
    struct hipDeviceProp_t properties;
    hipError_t cudaResultCode = hipGetDeviceCount(&deviceCount);
    if (cudaResultCode != hipSuccess) deviceCount = 0;
    /* machines with no GPUs can still report one emulation device */
    for (device = 0; device < deviceCount; ++device) {
        hipGetDeviceProperties(&properties, device);
        if (properties.major != 9999) { /* 9999 means emulation only */
            ++gpuDeviceCount;
        }
    }
    if (gpuDeviceCount) return 1;
    else return 0;
}

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

void syncQuESTEnv(QuESTEnv env){
    hipDeviceSynchronize();
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);
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
	
void getEnvironmentString(QuESTEnv env, char str[200]){
	
    // OpenMP can be hybridised with GPU in future, so this check is safe and worthwhile
    int ompStatus=0;
    int numThreads=1;
# ifdef _OPENMP
    ompStatus=1;
    numThreads=omp_get_max_threads(); 
# endif
	
    // there is no reporting of CUDA cores/threads/blocks currently (since non-trivial)
    sprintf(str, "CUDA=1 OpenMP=%d MPI=0 threads=%d ranks=1", ompStatus, numThreads);
}

void copyStateToGPU(Qureg qureg)
{
    if (DEBUG) printf("Copying data to GPU\n");
    hipMemcpy(qureg.deviceStateVec.real, qureg.stateVec.real, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.real)), hipMemcpyHostToDevice);
    hipMemcpy(qureg.deviceStateVec.imag, qureg.stateVec.imag, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.imag)), hipMemcpyHostToDevice);
    if (DEBUG) printf("Finished copying data to GPU\n");
}

void copyStateFromGPU(Qureg qureg)
{
    hipDeviceSynchronize();
    if (DEBUG) printf("Copying data from GPU\n");
    hipMemcpy(qureg.stateVec.real, qureg.deviceStateVec.real, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.real)), hipMemcpyDeviceToHost);
    hipMemcpy(qureg.stateVec.imag, qureg.deviceStateVec.imag, 
            qureg.numAmpsPerChunk*sizeof(*(qureg.deviceStateVec.imag)), hipMemcpyDeviceToHost);
    if (DEBUG) printf("Finished copying data from GPU\n");
}

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    if (DEBUG) printf("Copying data to GPU\n");
    hipMemcpy(&(qureg.deviceStateVec.real[startInd]), &(qureg.stateVec.real[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.real)), hipMemcpyHostToDevice);
    hipMemcpy(&(qureg.deviceStateVec.imag[startInd]), &(qureg.stateVec.imag[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.imag)), hipMemcpyHostToDevice);
    if (DEBUG) printf("Finished copying data to GPU\n");
}

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    hipDeviceSynchronize();
    if (DEBUG) printf("Copying data from GPU\n");
    hipMemcpy(&(qureg.stateVec.real[startInd]), &(qureg.deviceStateVec.real[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.real)), hipMemcpyDeviceToHost);
    hipMemcpy(&(qureg.stateVec.imag[startInd]), &(qureg.deviceStateVec.imag[startInd]), 
            numAmps*sizeof(*(qureg.deviceStateVec.imag)), hipMemcpyDeviceToHost);
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
    hipMemcpy(&el, &(qureg.deviceStateVec.real[index]), 
            sizeof(*(qureg.deviceStateVec.real)), hipMemcpyDeviceToHost);
    return el;
}

qreal statevec_getImagAmp(Qureg qureg, long long int index){
    qreal el=0;
    hipMemcpy(&el, &(qureg.deviceStateVec.imag[index]), 
            sizeof(*(qureg.deviceStateVec.imag)), hipMemcpyDeviceToHost);
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
                sscanf(line, REAL_SPECIFIER ", " REAL_SPECIFIER, &(stateVecReal[indexInChunk]),
                        &(stateVecImag[indexInChunk]));
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
    hipMalloc(&d_targs, targMemSize);
    hipMemcpy(d_targs, targs, targMemSize, hipMemcpyHostToDevice);
    
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
    hipMalloc(&d_uRe, uMemSize);
    hipMalloc(&d_uIm, uMemSize);
    hipMemcpy(d_uRe, uReFlat, uMemSize, hipMemcpyHostToDevice);
    hipMemcpy(d_uIm, uImFlat, uMemSize, hipMemcpyHostToDevice);
    
    // allocate device Wspace for thread-local {ampInds}, {reAmps}, {imAmps} (length: 1<<numTargs)
    long long int *d_ampInds;
    qreal *d_reAmps;
    qreal *d_imAmps;
    size_t gridSize = (size_t) threadsPerCUDABlock * CUDABlocks;
    int numTargAmps = uNumRows;
    hipMalloc(&d_ampInds, numTargAmps*gridSize * sizeof *d_ampInds);
    hipMalloc(&d_reAmps,  numTargAmps*gridSize * sizeof *d_reAmps);
    hipMalloc(&d_imAmps,  numTargAmps*gridSize * sizeof *d_imAmps);
    
    // call kernel
    statevec_multiControlledMultiQubitUnitaryKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, ctrlMask, d_targs, numTargs, d_uRe, d_uIm, d_ampInds, d_reAmps, d_imAmps, numTargAmps);
        
    // free kernel memory
    free(uReFlat);
    free(uImFlat);
    hipFree(d_targs);
    hipFree(d_uRe);
    hipFree(d_uIm);
    hipFree(d_ampInds);
    hipFree(d_reAmps);
    hipFree(d_imAmps);
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal zeroProb;
    hipMemcpy(&zeroProb, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    return zeroProb;
}

qreal statevec_findProbabilityOfZero(Qureg qureg, int measureQubit)
{
    qreal stateProb=0;
    
    // 1-qubit edge-case breaks below loop logic
    if (qureg.numQubitsInStateVec == 1) {
        qreal amp;
        hipMemcpy(&amp, qureg.deviceStateVec.real, sizeof(qreal), hipMemcpyDeviceToHost);
        stateProb += amp*amp;
        hipMemcpy(&amp, qureg.deviceStateVec.imag, sizeof(qreal), hipMemcpyDeviceToHost);
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&stateProb, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
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

// atomicAdd on floats/doubles isn't available on <6 CC devices, so we add it ourselves
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
static __inline__ __device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(val + __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

__global__ void statevec_calcProbOfAllOutcomesKernel(
    qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits
) {
    // each thread handles one amplitude (all amplitudes are involved)
    long long int ampInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (ampInd >= qureg.numAmpsTotal) return;
    
    qreal prob = (
        qureg.deviceStateVec.real[ampInd]*qureg.deviceStateVec.real[ampInd] + 
        qureg.deviceStateVec.imag[ampInd]*qureg.deviceStateVec.imag[ampInd]);
    
    // each amplitude contributes to one outcome
    long long int outcomeInd = 0;
    for (int q=0; q<numQubits; q++)
        outcomeInd += extractBit(qubits[q], ampInd) * (1LL << q);
    
    // each thread atomically writes directly to the global output.
    // this beat block-heirarchal atomic reductions in both global and shared memory!
    atomicAdd(&outcomeProbs[outcomeInd], prob);
}

void statevec_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {

    // copy qubits to GPU memory
    int* d_qubits;
    size_t mem_qubits = numQubits * sizeof *d_qubits;
    hipMalloc(&d_qubits, mem_qubits);
    hipMemcpy(d_qubits, qubits, mem_qubits, hipMemcpyHostToDevice);

    // create one thread for every amplitude
    int numThreadsPerBlock = 128;
    int numBlocks = ceil(qureg.numAmpsPerChunk / (qreal) numThreadsPerBlock);
    
    // create global GPU array for outcomeProbs
    qreal* d_outcomeProbs;
    long long int numOutcomes = (1LL << numQubits);
    size_t mem_outcomeProbs = numOutcomes * sizeof *d_outcomeProbs;
    hipMalloc(&d_outcomeProbs, mem_outcomeProbs);
    hipMemset(d_outcomeProbs, 0, mem_outcomeProbs);
    
    // populate per-block subarrays
    statevec_calcProbOfAllOutcomesKernel<<<numBlocks, numThreadsPerBlock>>>(
        d_outcomeProbs, qureg, d_qubits, numQubits);
        
    // copy outcomeProbs from GPU memory
    hipMemcpy(outcomeProbs, d_outcomeProbs, mem_outcomeProbs, hipMemcpyDeviceToHost);
    
    // free GPU memory
    hipFree(d_qubits);
    hipFree(d_outcomeProbs);
}

__global__ void densmatr_calcProbOfAllOutcomesKernel(
    qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits
) {
    // each thread handles one diagonal amplitude
    long long int diagInd = blockIdx.x*blockDim.x + threadIdx.x;
    long long int numDiags = (1LL << qureg.numQubitsRepresented);
    if (diagInd >= numDiags) return;
    
    long long int flatInd = (1 + numDiags)*diagInd;
    qreal prob = qureg.deviceStateVec.real[flatInd];   // im[flatInd] assumed ~ 0
    
    // each diagonal amplitude contributes to one outcome
    long long int outcomeInd = 0;
    for (int q=0; q<numQubits; q++)
        outcomeInd += extractBit(qubits[q], diagInd) * (1LL << q);
    
    // each thread atomically writes directly to the global output.
    // this beat block-heirarchal atomic reductions in both global and shared memory!
    atomicAdd(&outcomeProbs[outcomeInd], prob);
}

void densmatr_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {

    // copy qubits to GPU memory
    int* d_qubits;
    size_t mem_qubits = numQubits * sizeof *d_qubits;
    hipMalloc(&d_qubits, mem_qubits);
    hipMemcpy(d_qubits, qubits, mem_qubits, hipMemcpyHostToDevice);
    
    // create global array, with per-block subarrays
    int numThreadsPerBlock = 128;
    int numDiags = (1LL << qureg.numQubitsRepresented);
    int numBlocks = ceil(numDiags / (qreal) numThreadsPerBlock);
        
    // create global GPU array for outcomeProbs
    qreal* d_outcomeProbs;
    long long int numOutcomes = (1LL << numQubits);
    size_t mem_outcomeProbs = numOutcomes * sizeof *d_outcomeProbs;
    hipMalloc(&d_outcomeProbs, mem_outcomeProbs);
    hipMemset(d_outcomeProbs, 0, mem_outcomeProbs);
    
    // populate per-block subarrays
    densmatr_calcProbOfAllOutcomesKernel<<<numBlocks, numThreadsPerBlock>>>(
        d_outcomeProbs, qureg, d_qubits, numQubits);
        
    // copy outcomeProbs from GPU memory
    hipMemcpy(outcomeProbs, d_outcomeProbs, mem_outcomeProbs, hipMemcpyDeviceToHost);
    
    // free GPU memory
    hipFree(d_qubits);
    hipFree(d_outcomeProbs);
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    b.firstLevelReduction, 
                    b.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(b.firstLevelReduction), &(b.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal innerprod;
    hipMemcpy(&innerprod, b.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    bra.firstLevelReduction, 
                    bra.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(bra.firstLevelReduction), &(bra.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&innerProdReal, bra.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    bra.firstLevelReduction, 
                    bra.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(bra.firstLevelReduction), &(bra.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&innerProdImag, bra.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    pureState.firstLevelReduction, 
                    pureState.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(pureState.firstLevelReduction), &(pureState.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal fidelity;
    hipMemcpy(&fidelity, pureState.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    a.firstLevelReduction, 
                    a.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(a.firstLevelReduction), &(a.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal trace;
    hipMemcpy(&trace, a.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    
    qreal traceDensSquared;
    hipMemcpy(&traceDensSquared, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
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

/** Called once for every 16 amplitudes */
__global__ void densmatr_mixTwoQubitDepolarisingKernel(
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

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel) {
    
    if (depolLevel == 0)
        return;
    
    // assumes qubit2 > qubit1
    
    densmatr_mixTwoQubitDephasing(qureg, qubit1, qubit2, depolLevel);
    
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
    densmatr_mixTwoQubitDepolarisingKernel<<<CUDABlocks, threadsPerCUDABlock>>>(
        depolLevel, qureg.deviceStateVec.real, qureg.deviceStateVec.imag, numAmpsToVisit,
        part1, part2, part3, part4, part5, rowCol1, rowCol2);
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&expecReal, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&expecImag, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&expecReal, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
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
            hipDeviceSynchronize();    
            copySharedReduceBlock<<<numCUDABlocks, valuesPerCUDABlock/2, sharedMemSize>>>(
                    qureg.firstLevelReduction, 
                    qureg.secondLevelReduction, valuesPerCUDABlock); 
            hipDeviceSynchronize();    
            swapDouble(&(qureg.firstLevelReduction), &(qureg.secondLevelReduction));
        }
        numValuesToReduce = numValuesToReduce/maxReducedPerLevel;
    }
    hipMemcpy(&expecImag, qureg.firstLevelReduction, sizeof(qreal), hipMemcpyDeviceToHost);
    
    // return complex
    Complex expecVal;
    expecVal.real = expecReal;
    expecVal.imag = expecImag;
    return expecVal;
}

void agnostic_setDiagonalOpElems(DiagonalOp op, long long int startInd, qreal* real, qreal* imag, long long int numElems) {

    // update both RAM and VRAM, for consistency
    memcpy(&op.real[startInd], real, numElems * sizeof(qreal));
    memcpy(&op.imag[startInd], imag, numElems * sizeof(qreal));

    hipDeviceSynchronize();
    hipMemcpy(
        op.deviceOperator.real + startInd, 
        real,
        numElems * sizeof(*(op.deviceOperator.real)), 
        hipMemcpyHostToDevice);
    hipMemcpy(
        op.deviceOperator.imag + startInd,
        imag,
        numElems * sizeof(*(op.deviceOperator.imag)), 
        hipMemcpyHostToDevice);
}

__global__ void statevec_applyPhaseFuncOverridesKernel(
    Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int numTerms, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides, 
    int conj
) {
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=qureg.numAmpsPerChunk) return;

    // determine global amplitude index (non-distributed, so it's just local index)
    long long int globalAmpInd = index;

    // determine phase index of {qubits}
    long long int phaseInd = 0LL;
    if (encoding == UNSIGNED) {
        for (int q=0; q<numQubits; q++)
            phaseInd += (1LL << q) * extractBit(qubits[q], globalAmpInd);
    }
    else if (encoding == TWOS_COMPLEMENT) {
        for (int q=0; q<numQubits-1; q++) // use final qubit to indicate sign 
            phaseInd += (1LL << q) * extractBit(qubits[q], globalAmpInd);
        if (extractBit(qubits[numQubits-1], globalAmpInd) == 1)
            phaseInd -= (1LL << (numQubits-1));
    }

    // determine if this phase index has an overriden value (i < numOverrides)
    int i;
    for (i=0; i<numOverrides; i++)
        if (phaseInd == overrideInds[i])
            break;

    // determine phase from {coeffs}, {exponents} (unless overriden)
    qreal phase = 0;
    if (i < numOverrides)
        phase = overridePhases[i];
    else
        for (int t=0; t<numTerms; t++)
            phase += coeffs[t] * pow((qreal) phaseInd, (qreal) exponents[t]);
            
    // negate phase to conjugate operator 
    if (conj)
        phase *= -1;

    // modify amp to amp * exp(i phase) 
    qreal c = cos(phase);
    qreal s = sin(phase);
    qreal re = qureg.deviceStateVec.real[index];
    qreal im = qureg.deviceStateVec.imag[index];

    // = {re[amp] cos(phase) - im[amp] sin(phase)} + i {re[amp] sin(phase) + im[amp] cos(phase)}
    qureg.deviceStateVec.real[index] = re*c - im*s;
    qureg.deviceStateVec.imag[index] = re*s + im*c;
}

 void statevec_applyPhaseFuncOverrides(
     Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding,
     qreal* coeffs, qreal* exponents, int numTerms, 
     long long int* overrideInds, qreal* overridePhases, int numOverrides,
     int conj
 ) {
    // allocate device space for global list of {qubits}, {coeffs}, {exponents}, {overrideInds} and {overridePhases}
    int* d_qubits;                          size_t mem_qubits = numQubits * sizeof *d_qubits;
    qreal* d_coeffs;                        size_t mem_terms = numTerms * sizeof *d_coeffs;
    qreal* d_exponents;                 
    long long int* d_overrideInds;          size_t mem_inds = numOverrides * sizeof *d_overrideInds;
    qreal* d_overridePhases;                size_t mem_phas = numOverrides * sizeof *d_overridePhases;
    hipMalloc(&d_qubits, mem_qubits);      hipMemcpy(d_qubits, qubits, mem_qubits, hipMemcpyHostToDevice);
    hipMalloc(&d_coeffs, mem_terms);       hipMemcpy(d_coeffs, coeffs, mem_terms, hipMemcpyHostToDevice);
    hipMalloc(&d_exponents, mem_terms);    hipMemcpy(d_exponents, exponents, mem_terms, hipMemcpyHostToDevice);
    hipMalloc(&d_overrideInds, mem_inds);  hipMemcpy(d_overrideInds, overrideInds, mem_inds, hipMemcpyHostToDevice);
    hipMalloc(&d_overridePhases,mem_phas); hipMemcpy(d_overridePhases, overridePhases, mem_phas, hipMemcpyHostToDevice);

    // call kernel
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal) qureg.numAmpsPerChunk / threadsPerCUDABlock);
    statevec_applyPhaseFuncOverridesKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_qubits, numQubits, encoding, 
        d_coeffs, d_exponents, numTerms, 
        d_overrideInds, d_overridePhases, numOverrides,
        conj);

    // cleanup device memory 
    hipFree(d_qubits);
    hipFree(d_coeffs);
    hipFree(d_exponents);
    hipFree(d_overrideInds);
    hipFree(d_overridePhases);
}

__global__ void statevec_applyMultiVarPhaseFuncOverridesKernel(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int* numTermsPerReg, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    long long int *phaseInds,
    int conj
) {
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=qureg.numAmpsPerChunk) return;

    // determine global amplitude index (non-distributed, so it's just local index)
    long long int globalAmpInd = index;

    /*
     * each thread needs to write to a local:
     *      long long int phaseInds[numRegs];
     * but instead has access to shared array phaseInds, with below stride and offset
    */
    size_t stride = gridDim.x*blockDim.x;
    size_t offset = blockIdx.x*blockDim.x + threadIdx.x;

    // determine phase indices
    int flatInd = 0;
    if (encoding == UNSIGNED) {
        for (int r=0; r<numRegs; r++) {
            phaseInds[r*stride+offset] = 0LL;
            for (int q=0; q<numQubitsPerReg[r]; q++)
                phaseInds[r*stride+offset] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);
        }
    }
    else if  (encoding == TWOS_COMPLEMENT) {
        for (int r=0; r<numRegs; r++) {
            phaseInds[r*stride+offset] = 0LL;
            for (int q=0; q<numQubitsPerReg[r]-1; q++)  
                phaseInds[r*stride+offset] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);
            // use final qubit to indicate sign
            if (extractBit(qubits[flatInd++], globalAmpInd) == 1)
                phaseInds[r*stride+offset] -= (1LL << (numQubitsPerReg[r]-1)); 
        }
    }

    // determine if this phase index has an overriden value (i < numOverrides)
    int i;
    for (i=0; i<numOverrides; i++) {
        int found = 1;
        for (int r=0; r<numRegs; r++) {
            if (phaseInds[r*stride+offset] != overrideInds[i*numRegs+r]) {
                found = 0;
                break;
            }
        }
        if (found)
            break;
    }

    // compute the phase (unless overriden)
    qreal phase = 0;
    if (i < numOverrides)
        phase = overridePhases[i];
    else {
        flatInd = 0;
        for (int r=0; r<numRegs; r++) {
            for (int t=0; t<numTermsPerReg[r]; t++) {
                phase += coeffs[flatInd] * pow((qreal) phaseInds[r*stride+offset], (qreal) exponents[flatInd]);
                flatInd++;
            }
        }
    }
    
    // negate phase to conjugate operator 
    if (conj)
        phase *= -1;

    // modify amp to amp * exp(i phase) 
    qreal c = cos(phase);
    qreal s = sin(phase);
    qreal re = qureg.deviceStateVec.real[index];
    qreal im = qureg.deviceStateVec.imag[index];

    // = {re[amp] cos(phase) - im[amp] sin(phase)} + i {re[amp] sin(phase) + im[amp] cos(phase)}
    qureg.deviceStateVec.real[index] = re*c - im*s;
    qureg.deviceStateVec.imag[index] = re*s + im*c;
}

void statevec_applyMultiVarPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int* numTermsPerReg, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj
) {
    // determine size of arrays, for cloning into GPU memory
    size_t mem_numQubitsPerReg = numRegs * sizeof *numQubitsPerReg;
    size_t mem_numTermsPerReg = numRegs * sizeof *numTermsPerReg;
    size_t mem_overridePhases = numOverrides * sizeof *overridePhases;
    size_t mem_overrideInds = numOverrides * numRegs * sizeof *overrideInds;
    size_t mem_qubits = 0;
    size_t mem_coeffs = 0;  
    size_t mem_exponents = 0;
    for (int r=0; r<numRegs; r++) {
        mem_qubits += numQubitsPerReg[r] * sizeof *qubits;
        mem_coeffs += numTermsPerReg[r] * sizeof *coeffs;
        mem_exponents += numTermsPerReg[r] * sizeof *exponents;
    }

    // allocate global GPU memory
    int* d_qubits;                  hipMalloc(&d_qubits,           mem_qubits);
    qreal* d_coeffs;                hipMalloc(&d_coeffs,           mem_coeffs);
    qreal* d_exponents;             hipMalloc(&d_exponents,        mem_exponents);
    int* d_numQubitsPerReg;         hipMalloc(&d_numQubitsPerReg,  mem_numQubitsPerReg);
    int* d_numTermsPerReg;          hipMalloc(&d_numTermsPerReg,   mem_numTermsPerReg);
    long long int* d_overrideInds;  hipMalloc(&d_overrideInds,     mem_overrideInds);
    qreal* d_overridePhases;        hipMalloc(&d_overridePhases,   mem_overridePhases);

    // copy function args into GPU memory
    hipMemcpy(d_qubits, qubits,                    mem_qubits,             hipMemcpyHostToDevice);
    hipMemcpy(d_coeffs, coeffs,                    mem_coeffs,             hipMemcpyHostToDevice);
    hipMemcpy(d_exponents, exponents,              mem_exponents,          hipMemcpyHostToDevice);
    hipMemcpy(d_numQubitsPerReg, numQubitsPerReg,  mem_numQubitsPerReg,    hipMemcpyHostToDevice);
    hipMemcpy(d_numTermsPerReg, numTermsPerReg,    mem_numTermsPerReg,     hipMemcpyHostToDevice);
    hipMemcpy(d_overrideInds, overrideInds,        mem_overrideInds,       hipMemcpyHostToDevice);
    hipMemcpy(d_overridePhases, overridePhases,    mem_overridePhases,     hipMemcpyHostToDevice);

    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal) qureg.numAmpsPerChunk / threadsPerCUDABlock);

    // allocate thread-local working space {phaseInds}
    long long int *d_phaseInds;
    size_t gridSize = (size_t) threadsPerCUDABlock * CUDABlocks;
    hipMalloc(&d_phaseInds, numRegs*gridSize * sizeof *d_phaseInds);

    // call kernel
    statevec_applyMultiVarPhaseFuncOverridesKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_qubits, d_numQubitsPerReg, numRegs, encoding,
        d_coeffs, d_exponents, d_numTermsPerReg, 
        d_overrideInds, d_overridePhases, numOverrides,
        d_phaseInds, 
        conj);

    // free device memory
    hipFree(d_qubits);
    hipFree(d_coeffs);
    hipFree(d_exponents);
    hipFree(d_numQubitsPerReg);
    hipFree(d_numTermsPerReg);
    hipFree(d_overrideInds);
    hipFree(d_overridePhases);
    hipFree(d_phaseInds);
}

__global__ void statevec_applyParamNamedPhaseFuncOverridesKernel(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    enum phaseFunc phaseFuncName, qreal* params, int numParams,
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    long long int* phaseInds,
    int conj
) {
    long long int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index>=qureg.numAmpsPerChunk) return;

    // determine global amplitude index (non-distributed, so it's just local index)
    long long int globalAmpInd = index;

    /*
     * each thread needs to write to a local:
     *      long long int phaseInds[numRegs];
     * but instead has access to shared array phaseInds, with below stride and offset
    */
    size_t stride = gridDim.x*blockDim.x;
    size_t offset = blockIdx.x*blockDim.x + threadIdx.x;

    // determine phase indices
    if (encoding == UNSIGNED) {
        int flatInd = 0;
        for (int r=0; r<numRegs; r++) {
            phaseInds[r*stride+offset] = 0LL;
            for (int q=0; q<numQubitsPerReg[r]; q++)
                phaseInds[r*stride+offset] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);
        }
    }
    else if  (encoding == TWOS_COMPLEMENT) {
        int flatInd = 0;
        for (int r=0; r<numRegs; r++) {
            phaseInds[r*stride+offset] = 0LL;
            for (int q=0; q<numQubitsPerReg[r]-1; q++)  
                phaseInds[r*stride+offset] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);
            // use final qubit to indicate sign
            if (extractBit(qubits[flatInd++], globalAmpInd) == 1)
                phaseInds[r*stride+offset] -= (1LL << (numQubitsPerReg[r]-1));
        }
    }

    // determine if this phase index has an overriden value (i < numOverrides)
    int i;
    for (i=0; i<numOverrides; i++) {
        int found = 1;
        for (int r=0; r<numRegs; r++) {
            if (phaseInds[r*stride+offset] != overrideInds[i*numRegs+r]) {
                found = 0;
                break;
            }
        }
        if (found)
            break;
    }

    // compute the phase (unless overriden)
    qreal phase = 0;
    if (i < numOverrides)
        phase = overridePhases[i];
    else {
        // compute norm related phases
        if (phaseFuncName == NORM || phaseFuncName == INVERSE_NORM ||
            phaseFuncName == SCALED_NORM || phaseFuncName == SCALED_INVERSE_NORM ||
            phaseFuncName == SCALED_INVERSE_SHIFTED_NORM) {
            qreal norm = 0;
            if (phaseFuncName == SCALED_INVERSE_SHIFTED_NORM) {
                for (int r=0; r<numRegs; r++) {
                    qreal dif = phaseInds[r*stride+offset] - params[2+r];
                    norm += dif*dif;
                }
            }
            else
                for (int r=0; r<numRegs; r++)
                    norm += phaseInds[r*stride+offset]*phaseInds[r*stride+offset];
            norm = sqrt(norm);

            if (phaseFuncName == NORM)
                phase = norm;
            else if (phaseFuncName == INVERSE_NORM)
                phase = (norm == 0.)? params[0] : 1/norm; // smallest non-zero norm is 1
            else if (phaseFuncName == SCALED_NORM)
                phase = params[0] * norm;
            else if (phaseFuncName == SCALED_INVERSE_NORM || phaseFuncName == SCALED_INVERSE_SHIFTED_NORM)
                phase = (norm <= REAL_EPS)? params[1] : params[0] / norm; // unless shifted closer to zero
        }
        // compute product related phases
        else if (phaseFuncName == PRODUCT || phaseFuncName == INVERSE_PRODUCT ||
                 phaseFuncName == SCALED_PRODUCT || phaseFuncName == SCALED_INVERSE_PRODUCT) {

            qreal prod = 1;
            for (int r=0; r<numRegs; r++)
                prod *= phaseInds[r*stride+offset];

            if (phaseFuncName == PRODUCT)
                phase = prod;
            else if (phaseFuncName == INVERSE_PRODUCT)
                phase = (prod == 0.)? params[0] : 1/prod; // smallest non-zero prod is +- 1
            else if (phaseFuncName == SCALED_PRODUCT)
                phase = params[0] * prod;
            else if (phaseFuncName == SCALED_INVERSE_PRODUCT)
                phase = (prod == 0.)? params[1] : params[0] / prod;
        }
        // compute Euclidean distance related phases 
        else if (phaseFuncName == DISTANCE || phaseFuncName == INVERSE_DISTANCE ||
                 phaseFuncName == SCALED_DISTANCE || phaseFuncName == SCALED_INVERSE_DISTANCE ||
                 phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE) {

            qreal dist = 0;
            if (phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE) {
                for (int r=0; r<numRegs; r+=2) {
                    qreal dif = (phaseInds[r*stride+offset] - phaseInds[(r+1)*stride+offset] - params[2+r/2]);
                    dist += dif*dif;
                }
            }
            else
                for (int r=0; r<numRegs; r+=2) {
                    qreal dif = (phaseInds[(r+1)*stride+offset] - phaseInds[r*stride+offset]);
                    dist += dif*dif;
                }
            dist = sqrt(dist);

            if (phaseFuncName == DISTANCE)
                phase = dist;
            else if (phaseFuncName == INVERSE_DISTANCE)
                phase = (dist == 0.)? params[0] : 1/dist; // smallest non-zero dist is 1
            else if (phaseFuncName == SCALED_DISTANCE)
                phase = params[0] * dist;
            else if (phaseFuncName == SCALED_INVERSE_DISTANCE || phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE)
                phase = (dist <= REAL_EPS)? params[1] : params[0] / dist; // unless shifted closer
        }
    }
    
    
    // negate phase to conjugate operator 
    if (conj)
        phase *= -1;

    // modify amp to amp * exp(i phase) 
    qreal c = cos(phase);
    qreal s = sin(phase);
    qreal re = qureg.deviceStateVec.real[index];
    qreal im = qureg.deviceStateVec.imag[index];

    // = {re[amp] cos(phase) - im[amp] sin(phase)} + i {re[amp] sin(phase) + im[amp] cos(phase)}
    qureg.deviceStateVec.real[index] = re*c - im*s;
    qureg.deviceStateVec.imag[index] = re*s + im*c;
}

void statevec_applyParamNamedPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    enum phaseFunc phaseFuncName, qreal* params, int numParams,
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj 
) {
    // determine size of arrays, for cloning into GPU memory
    size_t mem_numQubitsPerReg = numRegs * sizeof *numQubitsPerReg;
    size_t mem_overridePhases = numOverrides * sizeof *overridePhases;
    size_t mem_overrideInds = numOverrides * numRegs * sizeof *overrideInds;
    size_t mem_params = numParams * sizeof *params;
    size_t mem_qubits = 0;
    for (int r=0; r<numRegs; r++)
        mem_qubits += numQubitsPerReg[r] * sizeof *qubits;

    // allocate global GPU memory
    int* d_qubits;                  hipMalloc(&d_qubits,           mem_qubits);
    int* d_numQubitsPerReg;         hipMalloc(&d_numQubitsPerReg,  mem_numQubitsPerReg);
    long long int* d_overrideInds;  hipMalloc(&d_overrideInds,     mem_overrideInds);
    qreal* d_overridePhases;        hipMalloc(&d_overridePhases,   mem_overridePhases);
    qreal* d_params = NULL;         if (numParams > 0) hipMalloc(&d_params, mem_params);

    // copy function args into GPU memory
    hipMemcpy(d_qubits, qubits,                    mem_qubits,             hipMemcpyHostToDevice);
    hipMemcpy(d_numQubitsPerReg, numQubitsPerReg,  mem_numQubitsPerReg,    hipMemcpyHostToDevice);
    hipMemcpy(d_overrideInds, overrideInds,        mem_overrideInds,       hipMemcpyHostToDevice);
    hipMemcpy(d_overridePhases, overridePhases,    mem_overridePhases,     hipMemcpyHostToDevice);
    if (numParams > 0)
        hipMemcpy(d_params, params, mem_params, hipMemcpyHostToDevice);

    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal) qureg.numAmpsPerChunk / threadsPerCUDABlock);

    // allocate thread-local working space {phaseInds}
    long long int *d_phaseInds;
    size_t gridSize = (size_t) threadsPerCUDABlock * CUDABlocks;
    hipMalloc(&d_phaseInds, numRegs*gridSize * sizeof *d_phaseInds);

    // call kernel
    statevec_applyParamNamedPhaseFuncOverridesKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_qubits, d_numQubitsPerReg, numRegs, encoding,
        phaseFuncName, d_params, numParams,
        d_overrideInds, d_overridePhases, numOverrides,
        d_phaseInds,
        conj);

    // free device memory
    hipFree(d_qubits);
    hipFree(d_numQubitsPerReg);
    hipFree(d_overrideInds);
    hipFree(d_overridePhases);
    hipFree(d_phaseInds);
    if (numParams > 0)
        hipFree(d_params);
}

void seedQuEST(QuESTEnv *env, unsigned long int *seedArray, int numSeeds) {

    // free existing seed array, if exists
    if (env->seeds != NULL)
        free(env->seeds);
        
    // record keys in permanent heap
    env->seeds = (unsigned long int*) malloc(numSeeds * sizeof *(env->seeds));
    for (int i=0; i<numSeeds; i++)
        (env->seeds)[i] = seedArray[i];
    env->numSeeds = numSeeds;
    
    // pass keys to Mersenne Twister seeder
    init_by_array(seedArray, numSeeds); 
}



#ifdef __cplusplus
}
#endif
