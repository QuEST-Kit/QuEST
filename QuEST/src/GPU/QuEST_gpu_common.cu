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
    cudaDeviceSynchronize();
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);
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
    
    int numThreadsPerBlock = 128;
    int numBlocks = ceil(op.numElemsPerChunk / (qreal) numThreadsPerBlock);
    agnostic_initDiagonalOpFromPauliHamilKernel<<<numBlocks, numThreadsPerBlock>>>(
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
