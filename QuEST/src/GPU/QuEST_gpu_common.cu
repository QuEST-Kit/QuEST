// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * GPU routines which are agnostic to the cuQuantum backend or the custom kernels
 *
 * @author Tyson Jones
 * @author Ania Brown (defined some original functions moved here)
 */


# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "QuEST_internal.h"
# include "mt19937ar.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

#ifdef USE_HIP
// Translate CUDA calls into HIP calls 
#include "cuda_to_hip.h"
#endif


#ifdef __cplusplus
extern "C" {
#endif



/*
 * GPU CONFIG
 */

const int NUM_THREADS_PER_BLOCK = 128;



/*
 * BACKEND AGNOSTIC COMPLEX ARITHMETIC
 */

#define PROD_REAL(aRe, aIm, bRe, bIm) \
    ((aRe)*(bRe) - (aIm)*(bIm))
#define PROD_IMAG(aRe, aIm, bRe, bIm) \
    ((aRe)*(bIm) + (aIm)*(bRe))

#ifdef USE_CUQUANTUM
    #define GET_AMP(aRe, aIm, qureg, i) \
        aRe = qureg.deviceCuStateVec[i].x; \
        aIm = qureg.deviceCuStateVec[i].y;
#else
    #define GET_AMP(aRe, aIm, qureg, i) \
        aRe = qureg.deviceStateVec.real[i]; \
        aIm = qureg.deviceStateVec.imag[i];
#endif

#ifdef USE_CUQUANTUM
    #define SET_AMP(qureg, i, aRe, aIm) \
        qureg.deviceCuStateVec[i].x = aRe; \
        qureg.deviceCuStateVec[i].y = aIm;
#else
    #define SET_AMP(qureg, i, aRe, aIm) \
        qureg.deviceStateVec.real[i] = aRe; \
        qureg.deviceStateVec.imag[i] = aIm;
#endif



/*
 * ENVIRONMENT MANAGEMENT
 */

int GPUExists(void) {
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

void syncQuESTEnv(QuESTEnv env){
    cudaDeviceSynchronize();
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void reportQuESTEnv(QuESTEnv env){
    printf("EXECUTION ENVIRONMENT:\n");
    printf("Running locally on one node with GPU ");
# ifdef USE_CUQUANTUM
    printf("via cuQuantum\n");
# else 
    printf("via custom kernels\n");
# endif
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

    // ascertain whether we're using cuQuantum or bespoke kernels
    int cuQuantumStatus=0;
# ifdef USE_CUQUANTUM
    cuQuantumStatus=1;
# endif
	
    // there is no reporting of CUDA cores/threads/blocks currently (since non-trivial)
    sprintf(str, "CUDA=1 cuQuantum=%d OpenMP=%d MPI=0 threads=%d ranks=1", cuQuantumStatus, ompStatus, numThreads);
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



/*
 * STATE INITIALISATION
 */

__global__ void densmatr_initPureStateKernel(Qureg dens, Qureg pure)
{
    long long int pureInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (pureInd >= pure.numAmpsPerChunk) return;

    qreal aRe, aIm;
    GET_AMP( aRe, aIm, pure, pureInd );

    for (long long int col=0; col<pure.numAmpsPerChunk; col++) {
        long long int densInd = col*pure.numAmpsPerChunk + pureInd;

        qreal bRe, bIm;
        GET_AMP( bRe, bIm, pure, col );
        bIm *= -1; // conjugated
        
        qreal cRe = PROD_REAL( aRe, aIm, bRe, bIm );
        qreal cIm = PROD_IMAG( aRe, aIm, bRe, bIm );

        SET_AMP( dens, densInd, cRe, cIm );
    }
}

void densmatr_initPureState(Qureg dens, Qureg pure)
{
    int CUDABlocks = ceil(pure.numAmpsPerChunk/ (qreal) NUM_THREADS_PER_BLOCK);
    densmatr_initPureStateKernel<<<CUDABlocks, NUM_THREADS_PER_BLOCK>>>(dens, pure);
}

__global__ void densmatr_setQuregToPauliHamilKernel(
    Qureg qureg, enum pauliOpType* pauliCodes, qreal* termCoeffs, int numSumTerms
) {
    long long int n = blockIdx.x*blockDim.x + threadIdx.x;
    if (n>=qureg.numAmpsPerChunk) return;
    
    // flattened {I,X,Y,Z} matrix elements, where [k] = [p][i][j]
    const int pauliRealElems[] = {   1,0, 0,1,   0,1, 1,0,   0,0, 0,0,   1,0, 0,-1  };
    const int pauliImagElems[] = {   0,0, 0,0,   0,0, 0,0,   0,-1,1,0,   0,0, 0,0   };

    // |n> = |c>|r>
    const int numQubits = qureg.numQubitsRepresented;
    const long long int r = n & ((1LL << numQubits) - 1);
    const long long int c = n >> numQubits;

    // new amplitude of |n>
    qreal elemRe = 0;
    qreal elemIm = 0;
    
    for (long long int t=0; t<numSumTerms; t++) {
        
        // pauliKronecker[r][c] = prod_q Pauli[q][q-th bit of r and c]
        int kronRe = 1;
        int kronIm = 0;
        long long int pInd = t * numQubits;
        
        for (int q=0; q<numQubits; q++) {
            
            // get element of Pauli matrix
            int i = (r >> q) & 1;
            int j = (c >> q) & 1;
            int p = (int) pauliCodes[pInd++];
            int k = (p<<2) + (i<<1) + j;
            int pauliRe = pauliRealElems[k]; 
            int pauliIm = pauliImagElems[k];
            
            // kron *= pauli
            int tmp = (pauliRe*kronRe) - (pauliIm*kronIm);
            kronIm = (pauliRe*kronIm) + (pauliIm*kronRe);
            kronRe = tmp;
        }
        
        // elem = sum_t coeffs[t] pauliKronecker[r][c]
        elemRe += termCoeffs[t] * kronRe;
        elemIm += termCoeffs[t] * kronIm;
    }
    
    // overwrite the density matrix entry
    SET_AMP( qureg, n, elemRe, elemIm );
}

void densmatr_setQuregToPauliHamil(Qureg qureg, PauliHamil hamil) {
    
    // copy hamil into GPU memory
    enum pauliOpType* d_pauliCodes;
    size_t mem_pauliCodes = hamil.numSumTerms * hamil.numQubits * sizeof *d_pauliCodes;
    cudaMalloc(&d_pauliCodes, mem_pauliCodes);
    cudaMemcpy(d_pauliCodes, hamil.pauliCodes, mem_pauliCodes, cudaMemcpyHostToDevice);
    
    qreal* d_termCoeffs;
    size_t mem_termCoeffs = hamil.numSumTerms * sizeof *d_termCoeffs;
    cudaMalloc(&d_termCoeffs, mem_termCoeffs);
    cudaMemcpy(d_termCoeffs, hamil.termCoeffs, mem_termCoeffs, cudaMemcpyHostToDevice);
    
    int numBlocks = ceil(qureg.numAmpsPerChunk / (qreal) NUM_THREADS_PER_BLOCK);
    densmatr_setQuregToPauliHamilKernel<<<numBlocks, NUM_THREADS_PER_BLOCK>>>(
        qureg, d_pauliCodes, d_termCoeffs, hamil.numSumTerms);

    // free tmp GPU memory
    cudaFree(d_pauliCodes);
    cudaFree(d_termCoeffs);
}



/*
 * DECOHERENCE
 */

__global__ void densmatr_mixTwoQubitDepolarisingKernel(
    Qureg qureg, qreal depolLevel, long long int part1, long long int part2, 
    long long int part3, long long int part4, long long int part5,
    long long int rowCol1, long long int rowCol2)
{
    long long int scanInd = blockIdx.x*blockDim.x + threadIdx.x;
    if (scanInd >= qureg.numAmpsPerChunk) return;
    
    long long int ind00 = (scanInd&part1) + ((scanInd&part2)<<1) + ((scanInd&part3)<<2) + ((scanInd&part4)<<3) + ((scanInd&part5)<<4);
    long long int ind01 = ind00 + rowCol1;
    long long int ind10 = ind00 + rowCol2;
    long long int ind11 = ind00 + rowCol1 + rowCol2;

    qreal re00, re01, re10, re11;
    qreal im00, im01, im10, im11;
    GET_AMP( re00, im00, qureg, ind00 );
    GET_AMP( re01, im01, qureg, ind01 );
    GET_AMP( re10, im10, qureg, ind10 );
    GET_AMP( re11, im11, qureg, ind11 );
    
    qreal realAvDepol = depolLevel * 0.25 * (re00 + re01 + re10 + re11);
    qreal imagAvDepol = depolLevel * 0.25 * (im00 + im01 + im10 + im11);
    qreal retain = 1 - depolLevel;

    SET_AMP( qureg, ind00, retain*re00 + realAvDepol, retain*im00 + imagAvDepol );
    SET_AMP( qureg, ind01, retain*re01 + realAvDepol, retain*im01 + imagAvDepol );
    SET_AMP( qureg, ind10, retain*re10 + realAvDepol, retain*im10 + imagAvDepol );
    SET_AMP( qureg, ind11, retain*re11 + realAvDepol, retain*im11 + imagAvDepol );
}

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel) {
    
    if (depolLevel == 0)
        return;

    densmatr_mixTwoQubitDephasing(qureg, qubit1, qubit2, depolLevel);
    
    // this code is painfully and unnecessarily verbose; it computes index offsets 
    // using bitwise logic, only to (within the kernel) effect those offsets with
    // integer arithmetic. Instead, one might as well obtain the ultimate indices
    // directly using bitwise logic within the kernel, passing no offsets at all.
    // We defer this to when the insertBits (and other bit twiddles) have been moved
    // into QuEST_gpu_common.cu.

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
    
    int CUDABlocks = ceil(numAmpsToVisit / (qreal) NUM_THREADS_PER_BLOCK);
    densmatr_mixTwoQubitDepolarisingKernel<<<CUDABlocks, NUM_THREADS_PER_BLOCK>>>(
        qureg, depolLevel, part1, part2, part3, part4, part5, rowCol1, rowCol2);
}



/*
 * MANAGING NON-QUREG TYPES 
 */

DiagonalOp agnostic_createDiagonalOp(int numQubits, QuESTEnv env) {

    DiagonalOp op;
    op.numQubits = numQubits;
    op.numElemsPerChunk = (1LL << numQubits) / env.numRanks;
    op.chunkId = env.rank;
    op.numChunks = env.numRanks;

    // allocate CPU memory (initialised to zero)
    op.real = (qreal*) calloc(op.numElemsPerChunk, sizeof(qreal));
    op.imag = (qreal*) calloc(op.numElemsPerChunk, sizeof(qreal));

    // check cpu memory allocation was successful
    validateDiagonalOpAllocation(&op, env, __func__);

    // allocate GPU memory
    size_t arrSize = op.numElemsPerChunk * sizeof(qreal);
    cudaMalloc(&(op.deviceOperator.real), arrSize);
    cudaMalloc(&(op.deviceOperator.imag), arrSize);

    // check gpu memory allocation was successful
    validateDiagonalOpGPUAllocation(&op, env, __func__);

    // initialise GPU memory to zero
    cudaMemset(op.deviceOperator.real, 0, arrSize);
    cudaMemset(op.deviceOperator.imag, 0, arrSize);

    return op;
}

void agnostic_destroyDiagonalOp(DiagonalOp op) {
    free(op.real);
    free(op.imag);
    cudaFree(op.deviceOperator.real);
    cudaFree(op.deviceOperator.imag);
}

void agnostic_syncDiagonalOp(DiagonalOp op) {
    cudaDeviceSynchronize();
    size_t mem_elems = op.numElemsPerChunk * sizeof *op.real;
    cudaMemcpy(op.deviceOperator.real, op.real, mem_elems, cudaMemcpyHostToDevice);
    cudaMemcpy(op.deviceOperator.imag, op.imag, mem_elems, cudaMemcpyHostToDevice);
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
                if ((elemInd >> q) & 1)
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
    cudaMalloc(&d_pauliCodes, mem_pauliCodes);
    cudaMemcpy(d_pauliCodes, hamil.pauliCodes, mem_pauliCodes, cudaMemcpyHostToDevice);
    
    qreal* d_termCoeffs;
    size_t mem_termCoeffs = hamil.numSumTerms * sizeof *d_termCoeffs;
    cudaMalloc(&d_termCoeffs, mem_termCoeffs);
    cudaMemcpy(d_termCoeffs, hamil.termCoeffs, mem_termCoeffs, cudaMemcpyHostToDevice);
    
    int numBlocks = ceil(op.numElemsPerChunk / (qreal) NUM_THREADS_PER_BLOCK);
    agnostic_initDiagonalOpFromPauliHamilKernel<<<numBlocks, NUM_THREADS_PER_BLOCK>>>(
        op, d_pauliCodes, d_termCoeffs, hamil.numSumTerms);
    
    // copy populated operator into to RAM
    cudaDeviceSynchronize();
    size_t mem_elems = op.numElemsPerChunk * sizeof *op.real;
    cudaMemcpy(op.real, op.deviceOperator.real, mem_elems, cudaMemcpyDeviceToHost);
    cudaMemcpy(op.imag, op.deviceOperator.imag, mem_elems, cudaMemcpyDeviceToHost);

    cudaFree(d_pauliCodes);
    cudaFree(d_termCoeffs);
}

void agnostic_setDiagonalOpElems(DiagonalOp op, long long int startInd, qreal* real, qreal* imag, long long int numElems) {

    // update both RAM and VRAM, for consistency
    memcpy(&op.real[startInd], real, numElems * sizeof(qreal));
    memcpy(&op.imag[startInd], imag, numElems * sizeof(qreal));

    cudaDeviceSynchronize();
    cudaMemcpy(
        op.deviceOperator.real + startInd, 
        real,
        numElems * sizeof(*(op.deviceOperator.real)), 
        cudaMemcpyHostToDevice);
    cudaMemcpy(
        op.deviceOperator.imag + startInd,
        imag,
        numElems * sizeof(*(op.deviceOperator.imag)), 
        cudaMemcpyHostToDevice);
}



#ifdef __cplusplus
}
#endif
