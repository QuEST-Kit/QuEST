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

__forceinline__ __device__ long long int getThreadInd() {
    return blockIdx.x*blockDim.x + threadIdx.x;
}



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
 * DEBUG 
 */

void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank){

    if (qureg.numQubitsInStateVec > 5) {
        printf("State reporting disabled for >5 qubits");
        return;
    }

    copyStateFromGPU(qureg); 

    // distributed GPU not yet supported
    if (reportRank != 0)
        return;

    printf("Reporting state from rank %d [\n", qureg.chunkId);
    printf("real, imag\n");

    for(long long int i=0; i<qureg.numAmpsPerChunk; i++)
        printf(
            REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", 
            qureg.stateVec.real[i], qureg.stateVec.imag[i]);

    printf("]\n");
}



/*
 * STATE INITIALISATION
 */

__global__ void densmatr_initPureStateKernel(Qureg dens, Qureg pure)
{
    long long int pureInd = getThreadInd();
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
    long long int n = getThreadInd();
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
 * CALCULATIONS
 */

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
    long long int ampInd = getThreadInd();
    if (ampInd >= qureg.numAmpsTotal) return;

    qreal ampRe, ampIm;
    GET_AMP( ampRe, ampIm, qureg, ampInd );
    qreal prob = ampRe*ampRe + ampIm*ampIm;
    
    // each amplitude contributes to one outcome
    long long int outcomeInd = 0;
    for (int q=0; q<numQubits; q++)
        outcomeInd += ((ampInd >> qubits[q]) & 1) << q;
    
    // each thread atomically writes directly to the global output.
    // this beat block-heirarchal atomic reductions in both global and shared memory!
    atomicAdd(&outcomeProbs[outcomeInd], prob);
}

void statevec_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {

    // copy qubits to GPU memory
    int* d_qubits;
    size_t mem_qubits = numQubits * sizeof *d_qubits;
    cudaMalloc(&d_qubits, mem_qubits);
    cudaMemcpy(d_qubits, qubits, mem_qubits, cudaMemcpyHostToDevice);

    // create one thread for every amplitude
    int numThreadsPerBlock = 128;
    int numBlocks = ceil(qureg.numAmpsPerChunk / (qreal) numThreadsPerBlock);
    
    // create global GPU array for outcomeProbs
    qreal* d_outcomeProbs;
    long long int numOutcomes = (1LL << numQubits);
    size_t mem_outcomeProbs = numOutcomes * sizeof *d_outcomeProbs;
    cudaMalloc(&d_outcomeProbs, mem_outcomeProbs);
    cudaMemset(d_outcomeProbs, 0, mem_outcomeProbs);
    
    // populate per-block subarrays
    statevec_calcProbOfAllOutcomesKernel<<<numBlocks, numThreadsPerBlock>>>(
        d_outcomeProbs, qureg, d_qubits, numQubits);
        
    // copy outcomeProbs from GPU memory
    cudaMemcpy(outcomeProbs, d_outcomeProbs, mem_outcomeProbs, cudaMemcpyDeviceToHost);
    
    // free GPU memory
    cudaFree(d_qubits);
    cudaFree(d_outcomeProbs);
}

__global__ void densmatr_calcProbOfAllOutcomesKernel(
    qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits
) {
    // each thread handles one diagonal amplitude
    long long int diagInd = getThreadInd();
    long long int numDiags = (1LL << qureg.numQubitsRepresented);
    if (diagInd >= numDiags) return;
    
    long long int flatInd = (1 + numDiags)*diagInd;

    qreal ampRe, ampIm;
    GET_AMP( ampRe, ampIm, qureg, flatInd );
    qreal prob = ampRe + 0*ampIm;   // ampIm assumed ~ 0, silence unused warning
    
    // each diagonal amplitude contributes to one outcome
    long long int outcomeInd = 0;
    for (int q=0; q<numQubits; q++)
        outcomeInd += ((diagInd >> qubits[q]) & 1) << q;
    
    // each thread atomically writes directly to the global output.
    // this beat block-heirarchal atomic reductions in both global and shared memory!
    atomicAdd(&outcomeProbs[outcomeInd], prob);
}

void densmatr_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {

    // copy qubits to GPU memory
    int* d_qubits;
    size_t mem_qubits = numQubits * sizeof *d_qubits;
    cudaMalloc(&d_qubits, mem_qubits);
    cudaMemcpy(d_qubits, qubits, mem_qubits, cudaMemcpyHostToDevice);
    
    // create global array, with per-block subarrays
    int numThreadsPerBlock = 128;
    int numDiags = (1LL << qureg.numQubitsRepresented);
    int numBlocks = ceil(numDiags / (qreal) numThreadsPerBlock);
        
    // create global GPU array for outcomeProbs
    qreal* d_outcomeProbs;
    long long int numOutcomes = (1LL << numQubits);
    size_t mem_outcomeProbs = numOutcomes * sizeof *d_outcomeProbs;
    cudaMalloc(&d_outcomeProbs, mem_outcomeProbs);
    cudaMemset(d_outcomeProbs, 0, mem_outcomeProbs);
    
    // populate per-block subarrays
    densmatr_calcProbOfAllOutcomesKernel<<<numBlocks, numThreadsPerBlock>>>(
        d_outcomeProbs, qureg, d_qubits, numQubits);
        
    // copy outcomeProbs from GPU memory
    cudaMemcpy(outcomeProbs, d_outcomeProbs, mem_outcomeProbs, cudaMemcpyDeviceToHost);
    
    // free GPU memory
    cudaFree(d_qubits);
    cudaFree(d_outcomeProbs);
}



/*
 * DECOHERENCE
 */

__global__ void densmatr_mixTwoQubitDepolarisingKernel(
    Qureg qureg, qreal depolLevel, long long int part1, long long int part2, 
    long long int part3, long long int part4, long long int part5,
    long long int rowCol1, long long int rowCol2)
{
    long long int scanInd = getThreadInd();
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
    long long int elemInd = getThreadInd();
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



/*
 * PHASE FUNCTIONS
 */

__forceinline__ __device__ int getBit(long long int num, int ind) {
    return (num >> ind) & 1;
}

__forceinline__ __device__ void applyPhaseToAmp(Qureg qureg, long long int index, qreal phase, int conj) 
{
    // negate phase to conjugate operator
    phase *= (1 - 2*conj);
    qreal c = cos(phase);
    qreal s = sin(phase);

    // modify amp to amp * exp(i phase) 
    qreal re, im;
    GET_AMP( re, im, qureg, index );
    SET_AMP( qureg, index, re*c - im*s, re*s + im*c );
}

__forceinline__ __device__ void getMultiRegPhaseIndOffsets(size_t* stride, size_t* offset)
{
    // determine the phaseInds tuple relevant for this thread
    *stride = gridDim.x*blockDim.x;
    *offset = blockIdx.x*blockDim.x + threadIdx.x;
}

__forceinline__ __device__ long long int getPhaseInd(long long int globalAmpInd, int* qubits, int numQubits, enum bitEncoding encoding)
{
    long long int phaseInd = 0LL;

    if (encoding == UNSIGNED) {
        for (int q=0; q<numQubits; q++)
            phaseInd += (1LL << q) * getBit(globalAmpInd, qubits[q]);
    }
    else if (encoding == TWOS_COMPLEMENT) {
        for (int q=0; q<numQubits-1; q++) // use final qubit to indicate sign 
            phaseInd += (1LL << q) * getBit(globalAmpInd, qubits[q]);
        if (getBit(globalAmpInd, qubits[numQubits-1]) == 1)
            phaseInd -= (1LL << (numQubits-1));
    }

    return phaseInd;
}

__forceinline__ __device__ void setMultiRegPhaseInds(
    long long int *phaseInds, long long int globalAmpInd, 
    int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding
) {
    // determine the phaseInds tuple relevant for this thread
    size_t stride, offset;
    getMultiRegPhaseIndOffsets(&stride, &offset);

    // determine phase indices
    if (encoding == UNSIGNED) {
        int flatInd = 0;
        for (int r=0; r<numRegs; r++) {
            phaseInds[r*stride+offset] = 0LL;
            for (int q=0; q<numQubitsPerReg[r]; q++)
                phaseInds[r*stride+offset] += (1LL << q) * getBit(globalAmpInd, qubits[flatInd++]);
        }
    }
    else if  (encoding == TWOS_COMPLEMENT) {
        int flatInd = 0;
        for (int r=0; r<numRegs; r++) {
            phaseInds[r*stride+offset] = 0LL;
            for (int q=0; q<numQubitsPerReg[r]-1; q++)  
                phaseInds[r*stride+offset] += (1LL << q) * getBit(globalAmpInd, qubits[flatInd++]);
            // use final qubit to indicate sign
            if (getBit(globalAmpInd, qubits[flatInd++]) == 1)
                phaseInds[r*stride+offset] -= (1LL << (numQubitsPerReg[r]-1));
        }
    }
}

__forceinline__ __device__ long long int getIndOfPhaseOverride(long long int phaseInd, long long int* overrideInds, int numOverrides)
{
    // determine if this phase index has an overriden value (i < numOverrides).
    // this is painfully inefficient and ill-suited to kernel lockstep, but we
    // assume there are very few (<5) overrides, so optimisation is unnecessary

    int i;
    for (i=0; i<numOverrides; i++)
        if (phaseInd == overrideInds[i])
            break;

    return i;
}

__forceinline__ __device__ long long int getIndOfMultiRegPhaseOverride(
    long long int* phaseInds, int numRegs, long long int* overrideInds, int numOverrides
) {
    // determine the phaseInds tuple relevant for this thread
    size_t stride, offset;
    getMultiRegPhaseIndOffsets(&stride, &offset);

    // determine if this phase index has an overriden value (i < numOverrides).
    // this is painfully inefficient and ill-suited to kernel lockstep, but we
    // assume there are very few (<5) overrides, so optimisation is unnecessary

    int i;
    for (i=0; i<numOverrides; i++) {

        // to be a match, every ind in the index-tuple must match
        int found = 1;
        for (int r=0; r<numRegs; r++) {
            if (phaseInds[r*stride+offset] != overrideInds[i*numRegs+r]) {
                found = 0;
                break;
            }
        }

        // breaking is mostly cosmetic; if even a single thread encounters no overrides,
        // it must iterate all overrides and all other lock-step threads must wait
        if (found)
            break;
    }

    // (i<numOverrides) indicates an override was found
    return i;
}

__forceinline__ __device__ qreal evalNormPhaseFunc(
    long long int* phaseInds, size_t stride, size_t offset,
    int numRegs, enum phaseFunc phaseFuncName, qreal* params, int numParams
) {
    // determine norm
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

    // determine phase via phase function

    if (phaseFuncName == NORM)
        return norm;

    if (phaseFuncName == INVERSE_NORM)
        return (norm == 0.)? params[0] : 1/norm; // smallest non-zero norm is 1

    if (phaseFuncName == SCALED_NORM)
        return params[0] * norm;

    if (
        phaseFuncName == SCALED_INVERSE_NORM || 
        phaseFuncName == SCALED_INVERSE_SHIFTED_NORM
    )
        return (norm <= REAL_EPS)? params[1] : params[0] / norm; // unless shifted closer to zero
}

__forceinline__ __device__ qreal evalProductPhaseFunc(
    long long int* phaseInds, size_t stride, size_t offset,
    int numRegs, enum phaseFunc phaseFuncName, qreal* params, int numParams
) {
    // determine product of phase indices
    qreal prod = 1;
    for (int r=0; r<numRegs; r++)
        prod *= phaseInds[r*stride+offset];

    // determine phase via phase function
    if (phaseFuncName == PRODUCT)
        return prod;

    if (phaseFuncName == INVERSE_PRODUCT)
        return (prod == 0.)? params[0] : 1/prod; // smallest non-zero prod is +- 1
    
    if (phaseFuncName == SCALED_PRODUCT)
        return params[0] * prod;

    if (phaseFuncName == SCALED_INVERSE_PRODUCT)
        return (prod == 0.)? params[1] : params[0] / prod;
}

__forceinline__ __device__ qreal evalDistancePhaseFunc(
    long long int* phaseInds, size_t stride, size_t offset,
    int numRegs, enum phaseFunc phaseFuncName, qreal* params, int numParams
) {
    // evaluate distance (depends on phase function)
    qreal dist = 0;
    if (phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE) {
        for (int r=0; r<numRegs; r+=2) {
            qreal dif = (phaseInds[r*stride+offset] - phaseInds[(r+1)*stride+offset] - params[2+r/2]);
            dist += dif*dif;
        }
    }
    else if (phaseFuncName == SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE) {
        for (int r=0; r<numRegs; r+=2) {
            qreal dif = (phaseInds[r*stride+offset] - phaseInds[(r+1)*stride+offset] - params[2+r+1]);
            dist += params[2+r] * dif*dif;
        }
    }
    else
        for (int r=0; r<numRegs; r+=2) {
            qreal dif = (phaseInds[(r+1)*stride+offset] - phaseInds[r*stride+offset]);
            dist += dif*dif;
        }

    // if sqrt() arg of distance would be negative, set it to zero, to subsequently be set to the divergence param
    if (dist < 0)
        dist = 0;

    dist = sqrt(dist);

    if (phaseFuncName == DISTANCE)
        return dist;

    if (phaseFuncName == INVERSE_DISTANCE)
        return (dist == 0.)? params[0] : 1/dist; // smallest non-zero dist is 1
    
    if (phaseFuncName == SCALED_DISTANCE)
        return params[0] * dist;
    
    if (
        phaseFuncName == SCALED_INVERSE_DISTANCE || 
        phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE || 
        phaseFuncName == SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE
    )
        return (dist <= REAL_EPS)? params[1] : params[0] / dist; // unless shifted closer
}

__forceinline__ __device__ qreal getPhaseFromParamNamedFunc(
    long long int* phaseInds, int numRegs, enum phaseFunc phaseFuncName, qreal* params, int numParams
) {
    // determine the phaseInds tuple relevant for the calling thread
    size_t stride, offset;
    getMultiRegPhaseIndOffsets(&stride, &offset);

    // compute norm related phases
    if (
        phaseFuncName == NORM || 
        phaseFuncName == INVERSE_NORM ||
        phaseFuncName == SCALED_NORM || 
        phaseFuncName == SCALED_INVERSE_NORM ||
        phaseFuncName == SCALED_INVERSE_SHIFTED_NORM
    )
        return evalNormPhaseFunc(phaseInds, stride, offset, numRegs, phaseFuncName, params, numParams);

    // compute product related phases
    if (
        phaseFuncName == PRODUCT || 
        phaseFuncName == INVERSE_PRODUCT ||
        phaseFuncName == SCALED_PRODUCT || 
        phaseFuncName == SCALED_INVERSE_PRODUCT
    )
        return evalProductPhaseFunc(phaseInds, stride, offset, numRegs, phaseFuncName, params, numParams);

    // compute Euclidean distance related phases 
    if (
        phaseFuncName == DISTANCE || 
        phaseFuncName == INVERSE_DISTANCE ||
        phaseFuncName == SCALED_DISTANCE || 
        phaseFuncName == SCALED_INVERSE_DISTANCE ||
        phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE || 
        phaseFuncName == SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE
    )
        return evalDistancePhaseFunc(phaseInds, stride, offset, numRegs, phaseFuncName, params, numParams);
}

__global__ void statevec_applyPhaseFuncOverridesKernel(
    Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int numTerms, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides, 
    int conj
) {
    long long int index = getThreadInd();
    if (index>=qureg.numAmpsPerChunk) return;

    // determine global amplitude index (non-distributed, so it's just local index)
    long long int globalAmpInd = index;

    // determine phaseInd from the qubit encoding
    long long int phaseInd = getPhaseInd(globalAmpInd, qubits, numQubits, encoding);

    // determine if this phase index has an overriden value (i < numOverrides)
    int overInd = getIndOfPhaseOverride(phaseInd, overrideInds, numOverrides);

    // determine phase from {coeffs}, {exponents} (unless overriden)
    qreal phase = 0;
    if (overInd < numOverrides)
        phase = overridePhases[overInd];
    else
        for (int t=0; t<numTerms; t++)
            phase += coeffs[t] * pow((qreal) phaseInd, (qreal) exponents[t]);
            
    // modify amp to amp * exp(i phase) 
    applyPhaseToAmp(qureg, index, phase, conj);
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
    cudaMalloc(&d_qubits, mem_qubits);      cudaMemcpy(d_qubits, qubits, mem_qubits, cudaMemcpyHostToDevice);
    cudaMalloc(&d_coeffs, mem_terms);       cudaMemcpy(d_coeffs, coeffs, mem_terms, cudaMemcpyHostToDevice);
    cudaMalloc(&d_exponents, mem_terms);    cudaMemcpy(d_exponents, exponents, mem_terms, cudaMemcpyHostToDevice);
    cudaMalloc(&d_overrideInds, mem_inds);  cudaMemcpy(d_overrideInds, overrideInds, mem_inds, cudaMemcpyHostToDevice);
    cudaMalloc(&d_overridePhases,mem_phas); cudaMemcpy(d_overridePhases, overridePhases, mem_phas, cudaMemcpyHostToDevice);

    // call kernel
    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal) qureg.numAmpsPerChunk / threadsPerCUDABlock);
    statevec_applyPhaseFuncOverridesKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_qubits, numQubits, encoding, 
        d_coeffs, d_exponents, numTerms, 
        d_overrideInds, d_overridePhases, numOverrides,
        conj);

    // cleanup device memory 
    cudaFree(d_qubits);
    cudaFree(d_coeffs);
    cudaFree(d_exponents);
    cudaFree(d_overrideInds);
    cudaFree(d_overridePhases);
}

__global__ void statevec_applyMultiVarPhaseFuncOverridesKernel(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int* numTermsPerReg, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    long long int *phaseInds,
    int conj
) {
    long long int index = getThreadInd();
    if (index>=qureg.numAmpsPerChunk) return;

    // determine global amplitude index (non-distributed, so it's just local index)
    long long int globalAmpInd = index;

    // determine phase indices (each thread has phaseInds[numRegs] sub-array)
    setMultiRegPhaseInds(phaseInds, globalAmpInd, qubits, numQubitsPerReg, numRegs, encoding);

    // determine if this phase index has an overriden value
    long long int overInd = getIndOfMultiRegPhaseOverride(phaseInds, numRegs, overrideInds, numOverrides);

    // either use the overriden phase...
    qreal phase = 0;
    if (overInd < numOverrides)
        phase = overridePhases[overInd];

    else {
        // else determine the phaseInds tuple relevant for this thread
        size_t stride, offset;
        getMultiRegPhaseIndOffsets(&stride, &offset);

        // and compute the phase from coeffs and exponents
        long long int flatInd = 0;
        for (int r=0; r<numRegs; r++) {
            for (int t=0; t<numTermsPerReg[r]; t++) {
                phase += coeffs[flatInd] * pow((qreal) phaseInds[r*stride+offset], (qreal) exponents[flatInd]);
                flatInd++;
            }
        }
    }
    
    // modify amp to amp * exp(i phase) 
    applyPhaseToAmp(qureg, index, phase, conj);
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
    int* d_qubits;                  cudaMalloc(&d_qubits,           mem_qubits);
    qreal* d_coeffs;                cudaMalloc(&d_coeffs,           mem_coeffs);
    qreal* d_exponents;             cudaMalloc(&d_exponents,        mem_exponents);
    int* d_numQubitsPerReg;         cudaMalloc(&d_numQubitsPerReg,  mem_numQubitsPerReg);
    int* d_numTermsPerReg;          cudaMalloc(&d_numTermsPerReg,   mem_numTermsPerReg);
    long long int* d_overrideInds;  cudaMalloc(&d_overrideInds,     mem_overrideInds);
    qreal* d_overridePhases;        cudaMalloc(&d_overridePhases,   mem_overridePhases);

    // copy function args into GPU memory
    cudaMemcpy(d_qubits, qubits,                    mem_qubits,             cudaMemcpyHostToDevice);
    cudaMemcpy(d_coeffs, coeffs,                    mem_coeffs,             cudaMemcpyHostToDevice);
    cudaMemcpy(d_exponents, exponents,              mem_exponents,          cudaMemcpyHostToDevice);
    cudaMemcpy(d_numQubitsPerReg, numQubitsPerReg,  mem_numQubitsPerReg,    cudaMemcpyHostToDevice);
    cudaMemcpy(d_numTermsPerReg, numTermsPerReg,    mem_numTermsPerReg,     cudaMemcpyHostToDevice);
    cudaMemcpy(d_overrideInds, overrideInds,        mem_overrideInds,       cudaMemcpyHostToDevice);
    cudaMemcpy(d_overridePhases, overridePhases,    mem_overridePhases,     cudaMemcpyHostToDevice);

    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal) qureg.numAmpsPerChunk / threadsPerCUDABlock);

    // allocate thread-local working space {phaseInds}
    long long int *d_phaseInds;
    size_t gridSize = (size_t) threadsPerCUDABlock * CUDABlocks;
    cudaMalloc(&d_phaseInds, numRegs*gridSize * sizeof *d_phaseInds);

    // call kernel
    statevec_applyMultiVarPhaseFuncOverridesKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_qubits, d_numQubitsPerReg, numRegs, encoding,
        d_coeffs, d_exponents, d_numTermsPerReg, 
        d_overrideInds, d_overridePhases, numOverrides,
        d_phaseInds, 
        conj);

    // free device memory
    cudaFree(d_qubits);
    cudaFree(d_coeffs);
    cudaFree(d_exponents);
    cudaFree(d_numQubitsPerReg);
    cudaFree(d_numTermsPerReg);
    cudaFree(d_overrideInds);
    cudaFree(d_overridePhases);
    cudaFree(d_phaseInds);
}

__global__ void statevec_applyParamNamedPhaseFuncOverridesKernel(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    enum phaseFunc phaseFuncName, qreal* params, int numParams,
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    long long int* phaseInds,
    int conj
) {
    long long int index = getThreadInd();
    if (index>=qureg.numAmpsPerChunk) return;

    // determine global amplitude index (non-distributed, so it's just local index)
    long long int globalAmpInd = index;

    // determine phase indices (each thread has phaseInds[numRegs] sub-array)
    setMultiRegPhaseInds(phaseInds, globalAmpInd, qubits, numQubitsPerReg, numRegs, encoding);

    // determine if this phase index has an overriden value
    long long int overInd = getIndOfMultiRegPhaseOverride(phaseInds, numRegs, overrideInds, numOverrides);

    // determine the phase, or the overriden one
    qreal phase = 0;
    if (overInd < numOverrides)
        phase = overridePhases[overInd];
    else
        phase = getPhaseFromParamNamedFunc(phaseInds, numRegs, phaseFuncName, params, numParams);
    
    // modify amp to amp * exp(i phase) 
    applyPhaseToAmp(qureg, index, phase, conj);
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
    int* d_qubits;                  cudaMalloc(&d_qubits,           mem_qubits);
    int* d_numQubitsPerReg;         cudaMalloc(&d_numQubitsPerReg,  mem_numQubitsPerReg);
    long long int* d_overrideInds;  cudaMalloc(&d_overrideInds,     mem_overrideInds);
    qreal* d_overridePhases;        cudaMalloc(&d_overridePhases,   mem_overridePhases);
    qreal* d_params = NULL;         if (numParams > 0) cudaMalloc(&d_params, mem_params);

    // copy function args into GPU memory
    cudaMemcpy(d_qubits, qubits,                    mem_qubits,             cudaMemcpyHostToDevice);
    cudaMemcpy(d_numQubitsPerReg, numQubitsPerReg,  mem_numQubitsPerReg,    cudaMemcpyHostToDevice);
    cudaMemcpy(d_overrideInds, overrideInds,        mem_overrideInds,       cudaMemcpyHostToDevice);
    cudaMemcpy(d_overridePhases, overridePhases,    mem_overridePhases,     cudaMemcpyHostToDevice);
    if (numParams > 0)
        cudaMemcpy(d_params, params, mem_params, cudaMemcpyHostToDevice);

    int threadsPerCUDABlock = 128;
    int CUDABlocks = ceil((qreal) qureg.numAmpsPerChunk / threadsPerCUDABlock);

    // allocate thread-local working space {phaseInds}
    long long int *d_phaseInds;
    size_t gridSize = (size_t) threadsPerCUDABlock * CUDABlocks;
    cudaMalloc(&d_phaseInds, numRegs*gridSize * sizeof *d_phaseInds);

    // call kernel
    statevec_applyParamNamedPhaseFuncOverridesKernel<<<CUDABlocks,threadsPerCUDABlock>>>(
        qureg, d_qubits, d_numQubitsPerReg, numRegs, encoding,
        phaseFuncName, d_params, numParams,
        d_overrideInds, d_overridePhases, numOverrides,
        d_phaseInds,
        conj);

    // free device memory
    cudaFree(d_qubits);
    cudaFree(d_numQubitsPerReg);
    cudaFree(d_overrideInds);
    cudaFree(d_overridePhases);
    cudaFree(d_phaseInds);
    if (numParams > 0)
        cudaFree(d_params);
}



#ifdef __cplusplus
}
#endif
