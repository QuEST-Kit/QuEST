// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * An implementation of QuEST's backend (../QuEST_internal.h) using NVIDIA's cuQuantum library.
 * This makes no use of the ComplexArray qureg.deviceStateVec, used by the bespoke GPU kernels,
 * which is not malloc'd in this deployment. Instead, this cuQuantum backend mallocs and uses
 * two dedicated arrays of 'cuAmp' complex primitives; qureg.cuStateVec (CPU memory) and
 * qureg.deviceCuStateVec (GPU memory). Note that some API functions are implemented directly
 * in QuEST_gpu_common.cu, in a way agnostic to cuQuantum vs bespoke kernels.
 *
 * @author Tyson Jones
 */

# include "QuEST.h"
# include "QuEST_gpu_common.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"

# include <cuComplex.h>
# include <vector>
# include <custatevec.h>
# include <thrust/device_ptr.h>
# include <thrust/sequence.h>
# include <thrust/iterator/zip_iterator.h>
# include <thrust/iterator/counting_iterator.h>
# include <thrust/iterator/transform_iterator.h>
# include <thrust/for_each.h>
# include <thrust/inner_product.h>



/*
 * TYPES AND ADAPTERS
 */

// precision-agnostic conversions between cuAmp, qreal and Complex
# if QuEST_PREC==1
    # define TO_CU_AMP(re, im) make_cuFloatComplex(re, im)
    # define cuAmpReal(amp) cuCrealf(amp)
    # define cuAmpImag(amp) cuCimagf(amp)
    # define cuAmpConj(amp) cuConjf(amp)
    # define CU_AMP_IN_STATE_PREC CUDA_C_F32
    # define CU_AMP_IN_MATRIX_PREC CUDA_C_64F
# elif QuEST_PREC==2
    # define TO_CU_AMP(re, im) make_cuDoubleComplex(re, im)
    # define cuAmpReal(amp) cuCreal(amp)
    # define cuAmpImag(amp) cuCimag(amp)
    # define cuAmpConj(amp) cuConj(amp)
    # define CU_AMP_IN_STATE_PREC CUDA_C_64F
    # define CU_AMP_IN_MATRIX_PREC CUDA_C_64F
# elif QuEST_PREC==4
    # define TO_CU_AMP(re, im) -1 // invalid precision config
    # define cuAmpReal(amp) -1
    # define cuAmpImag(amp) -1
    # define CU_AMP_IN_STATE_PREC void // invalid
    # define CU_AMP_IN_MATRIX_PREC void // invalid
#endif

// convenient operator overloads for cuAmp, for doing complex artihmetic.
// some of these are defined to be used by Thrust's backend, because we
// avoided Thrust's complex<qreal> (see QuEST_precision.h for explanation).
// notice we are manually performing the arithmetic using qreals and 
// re-packing the result into TO_CU_AMP, rather than using cuComplex.h's 
// functions like cuCadd(). This is because such functions are precision
// specific (grr) and do exactly the same thing themselves!
__host__ __device__ inline cuAmp operator - (const cuAmp& a) {
    return TO_CU_AMP(-cuAmpReal(a), -cuAmpImag(a));
}
__host__ __device__ inline cuAmp operator * (const cuAmp& a, const std::size_t n) {
    return TO_CU_AMP(n*cuAmpReal(a), n*cuAmpImag(a));
}
__host__ __device__ inline cuAmp operator + (const cuAmp& a, const cuAmp& b) {
    return TO_CU_AMP(cuAmpReal(a) + cuAmpReal(b), cuAmpImag(a) + cuAmpImag(b));
}
__host__ __device__ inline cuAmp operator * (const cuAmp& a, const cuAmp& b) {
    qreal aRe = cuAmpReal(a);
    qreal aIm = cuAmpImag(a);
    qreal bRe = cuAmpReal(b);
    qreal bIm = cuAmpImag(b);
    return TO_CU_AMP(aRe*bRe - aIm*bIm, aRe*bIm + aIm*bRe);
}

// convert user-facing Complex to cuQuantum-facing cuAmp
cuAmp toCuAmp(Complex c) {
    return TO_CU_AMP(c.real, c.imag);
}

// concise alias for row-wise flattened complex matrix
typedef std::vector<cuAmp> cuMatr;

// flatten ComplexMatrixN mIn to a cuMatr mOut
#define GET_cuMatr_FROM_ComplexMatrix( mOut, mIn, nQubits ) \
    long long int dim = (1LL << nQubits); \
    cuMatr mOut(dim*dim); \
    long long int i=0; \
    for (long long int r=0; r<(dim); r++) \
        for (long long int c=0; c<(dim); c++) \
            mOut[i++] = TO_CU_AMP(mIn.real[r][c], mIn.imag[r][c]);

// convert user-facing ComplexMatrixN to cuQuantum-facing cuMatr
cuMatr toCuMatr(ComplexMatrix2 mIn) {
    GET_cuMatr_FROM_ComplexMatrix(mOut, mIn, 1);
    return mOut;
}
cuMatr toCuMatr(ComplexMatrix4 mIn) {
    GET_cuMatr_FROM_ComplexMatrix(mOut, mIn, 2);
    return mOut;
}
cuMatr toCuMatr(ComplexMatrixN mIn) {
    GET_cuMatr_FROM_ComplexMatrix(mOut, mIn, mIn.numQubits);
    return mOut;
}

// convert QuEST backend masks back into user-input qubit lists (needed by cuQuantum)
std::vector<int> getIndsFromMask(long long int mask, int numBits) {
    std::vector<int> inds;
    for (int i=0; i<numBits; i++)
        if (mask & (1LL<<i))
            inds.push_back(i);
    return inds;
}



/*
 * CUQUANTUM MEMORY MANAGEMENT
 */

int GPUSupportsMemPools() {

    // consult only the first device (garuanteed already to exist)
    int deviceId;
    cudaGetDevice(&deviceId);
    int supports;
    cudaDeviceGetAttribute(&supports, cudaDevAttrMemoryPoolsSupported, deviceId);
    return supports;
}

int memPoolAlloc(void* ctx, void** ptr, size_t size, cudaStream_t stream) {
    cudaMemPool_t& pool = *static_cast<cudaMemPool_t*>(ctx);
    return cudaMallocFromPoolAsync(ptr, size, pool, stream); 
}
int memPoolFree(void* ctx, void* ptr, size_t size, cudaStream_t stream) {
    return cudaFreeAsync(ptr, stream); 
}

CuQuantumConfig* createCuConfig() {

    // create cuQuantumConfig in heap memory
    CuQuantumConfig* config = (CuQuantumConfig*) malloc(sizeof(CuQuantumConfig));

    // bind existing memory pool (does not need later manual freeing)
    int deviceId;
    cudaGetDevice(&deviceId);
    cudaDeviceGetMemPool(&(config->cuMemPool), deviceId);

    // create new custatevecHandle_t
    custatevecCreate(&(config->cuQuantumHandle));

    // create new cudaStream_t
    cudaStreamCreate(&(config->cuStream));

    // custatevecDeviceMemHandler_t needs no explicit creation

    // return config's heap pointer
    return config;
}

void initCuConfig(CuQuantumConfig* config) {

    // get existing memPool threshold, above which memory gets freed at every stream synch
    size_t currMaxMem;
    cudaMemPoolGetAttribute(config->cuMemPool, cudaMemPoolAttrReleaseThreshold, &currMaxMem); 

    // if memPool threshold smaller than 1 MiB = 16 qubits, extend it
    size_t desiredMaxMem = 16*(1<<15);
    if (currMaxMem < desiredMaxMem)
        cudaMemPoolSetAttribute(config->cuMemPool, cudaMemPoolAttrReleaseThreshold, &desiredMaxMem); 

    // bind mempool to deviceMemHandler
    config->cuMemHandler.ctx = &(config->cuMemPool);
    config->cuMemHandler.device_alloc = memPoolAlloc;
    config->cuMemHandler.device_free = memPoolFree;
    strcpy(config->cuMemHandler.name, "mempool");

    // bind deviceMemHandler to cuQuantum
    custatevecSetDeviceMemHandler(config->cuQuantumHandle, &(config->cuMemHandler));

    // bind stream to cuQuantum
    custatevecSetStream(config->cuQuantumHandle, config->cuStream);
}

void destroyCuConfig(CuQuantumConfig* config) {

    // free config's heap attributes
    cudaStreamDestroy(config->cuStream);
    custatevecDestroy(config->cuQuantumHandle);

    // don't need to free cuMemPool; it already existed
    // don't need to free cuMemHandler; it's a struct included in config's heap memory

    // free config's heap memory
    free(config);
}



/*
 * CUQUANTUM WRAPPERS (to reduce boilerplate)
 */

void custatevec_applyMatrix(Qureg qureg, std::vector<int> ctrls, std::vector<int> targs, cuMatr matr) {

    // do not adjoint matrix
    int adj = 0;

    // condition all ctrls on =1 state
    int* ctrlVals = nullptr;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    custatevecApplyMatrix(
        qureg.cuConfig->cuQuantumHandle, 
        qureg.deviceCuStateVec, CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        matr.data(), CU_AMP_IN_MATRIX_PREC, CUSTATEVEC_MATRIX_LAYOUT_ROW, adj, 
        targs.data(), targs.size(), 
        ctrls.data(), ctrlVals, ctrls.size(), 
        CUSTATEVEC_COMPUTE_DEFAULT,
        work, workSize);
}

void custatevec_applyDiagonal(Qureg qureg, std::vector<int> ctrls, std::vector<int> targs, cuAmp* elems) {

    // apply no permutation matrix
    custatevecIndex_t *perm = nullptr;

    // do not adjoint elems
    int adj = 0;

    // condition all ctrls on =1 state
    int* ctrlVals = nullptr;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    custatevecApplyGeneralizedPermutationMatrix(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec,
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec,
        perm, elems, CU_AMP_IN_MATRIX_PREC, adj, 
        targs.data(), targs.size(), 
        ctrls.data(), ctrlVals, ctrls.size(),
        work, workSize);
}



/*
 * THRUST WRAPPERS AND FUNCTORS (to reduce boilerplate, and for custom modification of statevectors)
 */

struct functor_prodOfConj : public thrust::binary_function<cuAmp,cuAmp,cuAmp>
{
    __host__ __device__ cuAmp operator()(cuAmp braAmp, cuAmp ketAmp) { 
        return  ketAmp * cuAmpConj(braAmp);
    }
};

struct functor_absOfDif : public thrust::binary_function<cuAmp,cuAmp,cuAmp>
{
    __host__ __device__ cuAmp operator()(cuAmp a, cuAmp b) { 
        qreal difRe = cuAmpReal(a) - cuAmpReal(b);
        qreal difIm = cuAmpImag(a) - cuAmpImag(b);
        qreal difAbs = difRe*difRe + difIm*difIm;
        return TO_CU_AMP(difAbs, 0);
    }
};

struct functor_absSquared : public thrust::unary_function<cuAmp,qreal>
{
    __host__ __device__ qreal operator()(cuAmp amp) {
        qreal re = cuAmpReal(amp);
        qreal im = cuAmpImag(amp);
        qreal absSq = re*re + im*im;
        return absSq;
    }
};

struct functor_prodOfAbsSquared : public thrust::binary_function<cuAmp,cuAmp,cuAmp>
{
    __host__ __device__ cuAmp operator()(cuAmp probAmp, cuAmp rawAmp) { 
        qreal re = cuAmpReal(probAmp);
        qreal im = cuAmpImag(probAmp);
        qreal absSq = re*re + im*im;
        cuAmp prod = rawAmp * TO_CU_AMP(absSq, 0); 
        return  prod;
    }
};

struct functor_scaleInd : public thrust::unary_function<long long int,long long int>
{
    const long long int factor;
    functor_scaleInd(long long int _factor) : factor(_factor) {}

    __host__ __device__ long long int operator()(long long int ind) {
        return factor * ind;
    }
};

struct functor_mapIndToAlternatingBlockScaledInd : public thrust::unary_function<long long int,long long int>
{
    long long int blockSize;
    long long int factor;
    int useOffsetBlocks;

    functor_mapIndToAlternatingBlockScaledInd(long long int _blockSize, long long int _factor, int _useOffsetBlocks) :
        blockSize(_blockSize), factor(_factor), useOffsetBlocks(_useOffsetBlocks) {}

    __host__ __device__ long long int operator()(long long int rawInd) {

        long long int blockNum = rawInd / blockSize;
        long long int subBlockInd = rawInd % blockSize;

        long long int blockStartInd = (blockNum * 2 * blockSize) + (useOffsetBlocks * blockSize);
        long long int blockifiedInd = blockStartInd + subBlockInd;

        long long int finalInd = factor * blockifiedInd;

        return finalInd;
    }
};

struct functor_setWeightedQureg
{
    // functor requires 3 complex scalars
    const cuAmp fac1;
    const cuAmp fac2;
    const cuAmp facOut;
    functor_setWeightedQureg(cuAmp _fac1, cuAmp _fac2, cuAmp _facOut) : 
        fac1(_fac1), fac2(_fac2), facOut(_facOut) {}

    // and modifies out to (facOut out + fac1 qureg1 + fac2 qureg2)
    template <typename Tuple> __host__ __device__ void operator()(Tuple t) {
        thrust::get<2>(t) = facOut*thrust::get<2>(t) + fac1*thrust::get<0>(t) + fac2*thrust::get<1>(t);
    }
};

struct functor_matrixColumnDotVector : public thrust::unary_function<long long int,cuAmp>
{
    cuAmp* matrix; // flattened column-wise
    cuAmp* vector;
    long long int dim;
    functor_matrixColumnDotVector(cuAmp* _matrix, cuAmp* _vector, long long int _dim) :
        matrix(_matrix), vector(_vector), dim(_dim) {}

    __host__ __device__ cuAmp operator()(long long int columnInd) {

        // safe to iterate serially, since #columnInd is exponentially growing
        cuAmp sum = TO_CU_AMP(0, 0);
        for (long long int rowInd=0; rowInd<dim; rowInd++)
            sum = sum + vector[rowInd] * cuAmpConj(matrix[columnInd*dim + rowInd]);
        
        return sum;
    }
};

struct functor_getDiagonalOpAmp : public thrust::unary_function<long long int,cuAmp>
{
    DiagonalOp op;
    functor_getDiagonalOpAmp(DiagonalOp _op) : op(_op) {}

    __host__ __device__ cuAmp operator()(long long int columnInd) {

        return TO_CU_AMP(
            op.deviceOperator.real[columnInd], 
            op.deviceOperator.imag[columnInd]);
    }
};

struct functor_multDiagOntoDensMatr
{
    cuAmp* matrix;
    DiagonalOp op;
    functor_multDiagOntoDensMatr(cuAmp* _matrix, DiagonalOp _op) : matrix(_matrix), op(_op) {}

    __host__ __device__ void operator()(long long int matrixInd) {

        long long int opDim = 1LL << op.numQubits;
        long long int opInd = matrixInd % opDim;
        cuAmp opAmp = TO_CU_AMP(
            op.deviceOperator.real[opInd], 
            op.deviceOperator.imag[opInd]);

        matrix[matrixInd] = opAmp * matrix[matrixInd];
    }
};

struct functor_collapseDensMatrToOutcome
{
    cuAmp* matrix;
    int numQubits;
    int outcome;
    qreal outcomeProb;
    int target;
    functor_collapseDensMatrToOutcome(cuAmp* _matrix, int _numQubits, int _outcome, qreal _outcomeProb, int _target) : 
        matrix(_matrix), numQubits(_numQubits), outcome(_outcome), outcomeProb(_outcomeProb), target(_target) {}
    
    __host__ __device__ void operator()(long long int ind) {

        // obtain bits |...b2...><...b1...|
        int b1 = (ind >> target) & 1;
        int b2 = (ind >> target >> numQubits) & 1;

        int f1 = !(b1 ^ b2);        // true if b1 == b2
        int f2 = !(b1 ^ outcome);   // true if b1 == outcome
        qreal fac = f1 * f2 / outcomeProb;  // 0 or 1/prob

        matrix[ind] = TO_CU_AMP(fac, 0) * matrix[ind];
    }
};

thrust::device_ptr<cuAmp> getStartPtr(Qureg qureg) {

    return thrust::device_pointer_cast(qureg.deviceCuStateVec);
}

thrust::device_ptr<cuAmp> getEndPtr(Qureg qureg) {

    return getStartPtr(qureg) + qureg.numAmpsTotal;
}

auto iter_strided(Qureg qureg, long long int stride, long long int strideIndInit) {

    // iterates qureg[i * stride] over i, from i = strideIndInit, until exceeding qureg
    auto rawIndIter = thrust::make_counting_iterator(strideIndInit);
    auto stridedIndIter = thrust::make_transform_iterator(rawIndIter, functor_scaleInd(stride));
    auto stridedAmpIter = thrust::make_permutation_iterator(getStartPtr(qureg), stridedIndIter);
    return stridedAmpIter;
}

auto getStartStridedAmpIter(Qureg qureg, long long int stride) {

    return iter_strided(qureg, stride, 0);
}

auto getEndStridedAmpIter(Qureg qureg, long long int stride) {

    long long int numIters = ceil(qureg.numAmpsTotal / (qreal) stride);
    return iter_strided(qureg, stride, numIters);
}

auto iter_blockStrided(Qureg qureg, long long int blockSize, int useOffsetBlocks, long long int stride, long long int rawIndInit) {

    // iterates qureg[(i mapped to alternating blocks) * stride] over i, from i = rawIndInit, until exceeding qureg
    auto functor = functor_mapIndToAlternatingBlockScaledInd(blockSize, stride, useOffsetBlocks);
    auto rawIndIter = thrust::make_counting_iterator(rawIndInit);
    auto stridedBlockIndIter = thrust::make_transform_iterator(rawIndIter, functor);
    auto stridedBlockAmpIter = thrust::make_permutation_iterator(getStartPtr(qureg), stridedBlockIndIter);
    return stridedBlockAmpIter;
}

auto getStartBlockStridedAmpIter(Qureg qureg, long long int blockSize, int useOffsetBlocks, long long int stride) {

    return iter_blockStrided(qureg, blockSize, useOffsetBlocks, stride, 0);
}

auto getEndBlockStridedAmpIter(Qureg qureg, long long int blockSize, int useOffsetBlocks, long long int stride) {

    long long int numAltBlockAmps = qureg.numAmpsTotal / 2;
    long long int numIters = ceil(numAltBlockAmps / (qreal) stride);
    return iter_blockStrided(qureg, blockSize, useOffsetBlocks, stride, numIters);
}

auto iter_diagonalOp(DiagonalOp op, long long int initInd) {

    auto rawIndIter = thrust::make_counting_iterator(initInd);
    auto diagAmpIter = thrust::make_transform_iterator(rawIndIter, functor_getDiagonalOpAmp(op));
    return diagAmpIter;
}

auto getStartDiagonalOpAmpIter(DiagonalOp op) {

    return iter_diagonalOp(op, 0);
}

auto getEndDiagonalOpAmpIter(DiagonalOp op) {

    return iter_diagonalOp(op, op.numElemsPerChunk);
}



/*
 * Start C-linkage, so that below functions may only use C signatures (not C++)
 */
#ifdef __cplusplus
extern "C" {
#endif



/* 
 * ENVIRONMENT MANAGEMENT
 */

QuESTEnv createQuESTEnv(void) {
    validateGPUExists(GPUExists(), __func__);
    validateGPUIsCuQuantumCompatible(GPUSupportsMemPools(),__func__);
    
    QuESTEnv env;
    env.rank=0;
    env.numRanks=1;
    
    env.seeds = NULL;
    env.numSeeds = 0;
    seedQuESTDefault(&env);

    // prepare cuQuantum with automatic workspaces
    env.cuConfig = createCuConfig();
    initCuConfig(env.cuConfig);

    return env;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);

    // finalise cuQuantum
    destroyCuConfig(env.cuConfig);
}



/* 
 * QUREG CREATION AND AMP SET/GET
 */

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env)
{   
    // set standard fields
    long long int numAmps = 1LL << numQubits;
    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsPerChunk = numAmps;
    qureg->numAmpsTotal = numAmps;
    qureg->chunkId = 0;
    qureg->numChunks = 1;
    qureg->isDensityMatrix = 0;

    // copy env's cuQuantum config handle
    qureg->cuConfig = env.cuConfig;

    // allocate user-facing CPU memory
    qureg->stateVec.real = (qreal*) malloc(numAmps * sizeof(qureg->stateVec.real));
    qureg->stateVec.imag = (qreal*) malloc(numAmps * sizeof(qureg->stateVec.imag));
    validateQuregAllocation(qureg, env, __func__);

    // allocate cuQuantum GPU memory (unvalidated)
    cudaMalloc( &(qureg->deviceCuStateVec), numAmps * sizeof(*(qureg->deviceCuStateVec)) );

    // allocate private cuQuantum CPU memory (for exchanging with GPU memory)
    qureg->cuStateVec = (cuAmp*) malloc(numAmps * sizeof(*(qureg->cuStateVec)));
}

void statevec_destroyQureg(Qureg qureg, QuESTEnv env)
{
    // free user-facing CPU memory 
    free(qureg.stateVec.real);
    free(qureg.stateVec.imag);
    
    // free private cuQuantum CPU memory
    free(qureg.cuStateVec);

    // free cuQuantum GPU memory
    cudaFree(qureg.deviceCuStateVec);
}

void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps)
{
    // slowly manually overwrite subset of private cuQuantum CPU memory
    for (long long int i=0; i<numAmps; i++)
        qureg.cuStateVec[i+startInd] = TO_CU_AMP(reals[i], imags[i]);

    // cuda-copy subset to GPU memory subset
    cudaDeviceSynchronize();
    cudaMemcpy(
        &(qureg.deviceCuStateVec[startInd]), 
        &(qureg.cuStateVec[startInd]), 
        numAmps * sizeof(cuAmp), cudaMemcpyHostToDevice);
}

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    statevec_setAmps(qureg, startInd, &(qureg.stateVec.real[startInd]), &(qureg.stateVec.imag[startInd]), numAmps);
}

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps)
{
    // cuda-copy subset of GPU memory to private cuQuantum CPU memory
    cudaDeviceSynchronize();
    cudaMemcpy(
        &(qureg.cuStateVec[startInd]), 
        &(qureg.deviceCuStateVec[startInd]), 
        numAmps * sizeof(*(qureg.cuStateVec)), 
        cudaMemcpyDeviceToHost);

    // slowly manually overwrite public CPU memory from private
    for (long long int i=startInd; i<(startInd+numAmps); i++) {
        qureg.stateVec.real[i] = cuAmpReal(qureg.cuStateVec[i]);
        qureg.stateVec.imag[i] = cuAmpImag(qureg.cuStateVec[i]);
    }
}

void copyStateToGPU(Qureg qureg)
{
    statevec_copySubstateToGPU(qureg, 0, qureg.numAmpsTotal);
}

void copyStateFromGPU(Qureg qureg)
{
    statevec_copySubstateFromGPU(qureg, 0, qureg.numAmpsTotal);
}

void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg)
{
    // directly cuda-copy the GPU memory 
    cudaDeviceSynchronize();
    cudaMemcpy(
        targetQureg.deviceCuStateVec,
        copyQureg.deviceCuStateVec,
        copyQureg.numAmpsTotal * sizeof(cuAmp),
        cudaMemcpyDeviceToDevice);
}

qreal statevec_getRealAmp(Qureg qureg, long long int index)
{
    cuAmp amp;
    cudaDeviceSynchronize();
    cudaMemcpy(&amp, &(qureg.deviceCuStateVec[index]), sizeof(cuAmp), cudaMemcpyDeviceToHost);
    return cuAmpReal(amp);
}

qreal statevec_getImagAmp(Qureg qureg, long long int index)
{
    cuAmp amp;
    cudaDeviceSynchronize();
    cudaMemcpy(&amp, &(qureg.deviceCuStateVec[index]), sizeof(cuAmp), cudaMemcpyDeviceToHost);
    return cuAmpImag(amp);
}



/*
 * STATE INITIALISATION
 */

void densmatr_initPlusState(Qureg qureg)
{
    qreal val = 1./(1LL << qureg.numQubitsRepresented);
    thrust::fill(getStartPtr(qureg), getEndPtr(qureg), TO_CU_AMP(val, 0));
}

void statevec_initZeroState(Qureg qureg)
{
    custatevecInitializeStateVector(
        qureg.cuConfig->cuQuantumHandle, 
        qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, 
        qureg.numQubitsInStateVec, 
        CUSTATEVEC_STATE_VECTOR_TYPE_ZERO);
}

void statevec_initBlankState(Qureg qureg)
{
    // init to |0> = {1, 0, 0, ...}
    statevec_initZeroState(qureg);

    // overwrite amps[0] to 0
    qreal re[] = {0};
    qreal im[] = {0};
    statevec_setAmps(qureg, 0, re, im, 1);
}

void statevec_initPlusState(Qureg qureg)
{
    custatevecInitializeStateVector(
        qureg.cuConfig->cuQuantumHandle,
        qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, 
        qureg.numQubitsInStateVec, 
        CUSTATEVEC_STATE_VECTOR_TYPE_UNIFORM);
}

void statevec_initClassicalState(Qureg qureg, long long int stateInd)
{
    // init to |null> = {0, 0, 0, ...}
    statevec_initBlankState(qureg);

    // overwrite amps[stateInd] to 1
    qreal re[] = {1};
    qreal im[] = {0};
    statevec_setAmps(qureg, stateInd, re, im, 1);
}

void densmatr_initClassicalState(Qureg qureg, long long int stateInd)
{
    // init to |null> = {0, 0, 0, ...}
    statevec_initBlankState(qureg);

    // index of the desired state in the flat density matrix
    long long int densityDim = 1LL << qureg.numQubitsRepresented;
    long long int densityInd = (densityDim + 1)*stateInd;

    // overwrite amps[densityInd] to 1
    qreal re[] = {1};
    qreal im[] = {0};
    statevec_setAmps(qureg, densityInd, re, im, 1);
}

void statevec_initDebugState(Qureg qureg)
{
    // |n> -> (.2n + (.2n+.1)i) |n>
    cuAmp init = TO_CU_AMP(0,  .1);
    cuAmp step = TO_CU_AMP(.2, .2);
    thrust::sequence(getStartPtr(qureg), getEndPtr(qureg), init, step);
}

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out)
{
    // out ->  facOut + fac1 qureg1 + fac2 qureg2
    thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(getStartPtr(qureg1), getStartPtr(qureg2), getStartPtr(out))),
        thrust::make_zip_iterator(thrust::make_tuple(getEndPtr(qureg1),   getEndPtr(qureg2),   getEndPtr(out))),
        functor_setWeightedQureg(toCuAmp(fac1), toCuAmp(fac2), toCuAmp(facOut)));
}



/*
 * OPERATORS
 */

void statevec_compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta) 
{
    cuAmp a = toCuAmp(alpha);
    cuAmp b = toCuAmp(beta);
    cuMatr matrix{
        a, -cuAmpConj(b),
        b,  cuAmpConj(a)
    };
    custatevec_applyMatrix(qureg, {}, {targetQubit}, matrix);
}

void statevec_controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta) 
{
    cuAmp a = toCuAmp(alpha);
    cuAmp b = toCuAmp(beta);
    cuMatr matrix{
        a, -cuAmpConj(b),
        b,  cuAmpConj(a)
    };
    custatevec_applyMatrix(qureg, {controlQubit}, {targetQubit}, matrix);
}

void statevec_unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u)
{
    custatevec_applyMatrix(qureg, {}, {targetQubit}, toCuMatr(u));
}

void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u)
{
    std::vector<int> c = getIndsFromMask(ctrlMask,qureg.numQubitsInStateVec);
    std::vector<int> t(targs,targs+numTargs); 
    custatevec_applyMatrix(qureg, c, t, toCuMatr(u));
}

void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int q1, int q2, ComplexMatrix4 u)
{
    std::vector<int> c = getIndsFromMask(ctrlMask,qureg.numQubitsInStateVec);
    custatevec_applyMatrix(qureg, c, {q1,q2}, toCuMatr(u));
}

void statevec_controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u)
{
    custatevec_applyMatrix(qureg, {controlQubit}, {targetQubit}, toCuMatr(u));
}

void statevec_multiControlledUnitary(Qureg qureg, long long int ctrlQubitsMask, long long int ctrlFlipMask, int targetQubit, ComplexMatrix2 u)
{
    int targs[] = {targetQubit};
    std::vector<int> ctrlInds = getIndsFromMask(ctrlQubitsMask,qureg.numQubitsInStateVec);
    std::vector<int> ctrlVals(ctrlInds.size());
    for (size_t i=0; i<ctrlInds.size(); i++)
        ctrlVals[i] = !(ctrlFlipMask & (1LL<<ctrlInds[i]));

    custatevecApplyMatrix(
        qureg.cuConfig->cuQuantumHandle,
        qureg.deviceCuStateVec, CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        toCuMatr(u).data(), CU_AMP_IN_MATRIX_PREC, CUSTATEVEC_MATRIX_LAYOUT_ROW, 0, 
        targs, 1, ctrlInds.data(), ctrlVals.data(), ctrlInds.size(), 
        CUSTATEVEC_COMPUTE_DEFAULT, nullptr, 0);
}

void statevec_pauliX(Qureg qureg, int targetQubit) 
{
    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuMatr matrix{
        a0, a1,
        a1, a0
    };
    custatevec_applyMatrix(qureg, {}, {targetQubit}, matrix);
}

void statevec_pauliY(Qureg qureg, int targetQubit) 
{
    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp aI = TO_CU_AMP(0, 1);
    cuMatr matrix{
        a0, -aI,
        aI,  a0
    };
    custatevec_applyMatrix(qureg, {}, {targetQubit}, matrix);
}

void statevec_pauliYConj(Qureg qureg, int targetQubit) 
{
    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp aI = TO_CU_AMP(0, 1);
    cuMatr matrix{
         a0, aI,
        -aI, a0
    };
    custatevec_applyMatrix(qureg, {}, {targetQubit}, matrix);
}

void statevec_controlledPauliY(Qureg qureg, int controlQubit, int targetQubit)
{
    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp aI = TO_CU_AMP(0, 1);
    cuMatr matrix{
        a0, -aI,
        aI,  a0
    };
    custatevec_applyMatrix(qureg, {controlQubit}, {targetQubit}, matrix);
}

void statevec_controlledPauliYConj(Qureg qureg, int controlQubit, int targetQubit)
{
    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp aI = TO_CU_AMP(0, 1);
    cuMatr matrix{
         a0, aI,
        -aI, a0
    };
    custatevec_applyMatrix(qureg, {controlQubit}, {targetQubit}, matrix);
}

void statevec_phaseShiftByTerm(Qureg qureg, int targetQubit, Complex term)
{   
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuAmp a2 = toCuAmp(term);
    cuAmp elems[] = {a1, a2};

    custatevec_applyDiagonal(qureg, {}, {targetQubit}, elems);
}

void statevec_controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle)
{
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuAmp a2 = TO_CU_AMP(cos(angle), sin(angle));
    cuAmp elems[] = {a1, a2};

    custatevec_applyDiagonal(qureg, {idQubit1}, {idQubit2}, elems);
}

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle)
{   
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuAmp a2 = TO_CU_AMP(cos(angle), sin(angle));
    cuAmp elems[] = {a1, a2};

    std::vector<int> targs{controlQubits[0]};
    std::vector<int> ctrls(controlQubits + 1, controlQubits + numControlQubits); 
    custatevec_applyDiagonal(qureg, ctrls, targs, elems);
}

void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle)
{   
    qreal theta = - angle/2.;
    std::vector<int> targs = getIndsFromMask(mask, qureg.numQubitsInStateVec);
    std::vector<custatevecPauli_t> paulis(targs.size(), CUSTATEVEC_PAULI_Z);

    custatevecApplyPauliRotation(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        theta, paulis.data(), targs.data(), targs.size(),
        nullptr, nullptr, 0);
}

void statevec_multiControlledMultiRotateZ(Qureg qureg, long long int ctrlMask, long long int targMask, qreal angle)
{   
    qreal theta = - angle/2.;
    std::vector<int> ctrls = getIndsFromMask(ctrlMask, qureg.numQubitsInStateVec);
    std::vector<int> targs = getIndsFromMask(targMask, qureg.numQubitsInStateVec);
    std::vector<custatevecPauli_t> paulis(targs.size(), CUSTATEVEC_PAULI_Z);

    custatevecApplyPauliRotation(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        theta, paulis.data(), targs.data(), targs.size(),
        ctrls.data(), nullptr, ctrls.size());
}

void statevec_controlledPhaseFlip(Qureg qureg, int idQubit1, int idQubit2)
{
    cuAmp a1 = TO_CU_AMP( 1, 0);
    cuAmp a2 = TO_CU_AMP(-1, 0);
    cuAmp elems[] = {a1, a2};

    custatevec_applyDiagonal(qureg, {idQubit1}, {idQubit2}, elems);
}

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits)
{
    cuAmp a1 = TO_CU_AMP( 1, 0);
    cuAmp a2 = TO_CU_AMP(-1, 0);
    cuAmp elems[] = {a1, a2};

    std::vector<int> targs{controlQubits[0]};
    std::vector<int> ctrls(controlQubits + 1, controlQubits + numControlQubits); 
    custatevec_applyDiagonal(qureg, ctrls, targs, elems);
}

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) 
{
    int2 targPairs[] = {{qb1, qb2}}; 
    int numPairs = 1;

    custatevecSwapIndexBits(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        targPairs, numPairs,
        nullptr, nullptr, 0);
}

void statevec_hadamard(Qureg qureg, int targetQubit) 
{
    cuAmp a = TO_CU_AMP(1/sqrt(2.), 0);
    cuMatr matrix{
        a,  a,
        a, -a
    };
    custatevec_applyMatrix(qureg, {}, {targetQubit}, matrix);
}

void statevec_controlledNot(Qureg qureg, int controlQubit, int targetQubit)
{
    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuMatr matrix{
        a0, a1,
        a1, a0
    };
    custatevec_applyMatrix(qureg, {controlQubit}, {targetQubit}, matrix);
}

void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask)
{
    // this operator can be effected in one-shot using a custom kernel, but we here
    // instead resort to slowly (by at most a factor #targs) effect it as a sequence
    // of one-target multi-ctrl NOT gates.

    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuMatr matrix{
        a0, a1,
        a1, a0
    };
    std::vector<int> ctrls = getIndsFromMask(ctrlMask, qureg.numQubitsInStateVec);
    std::vector<int> targs = getIndsFromMask(targMask, qureg.numQubitsInStateVec);
    for (int targ : targs)
        custatevec_applyMatrix(qureg, ctrls, {targ}, matrix);
}

void statevec_applySubDiagonalOp(Qureg qureg, int* targs, SubDiagonalOp op, int conj)
{
    // sneakily leverage the CPU cuQuantum memory in order to convert op
    // (as separate arrays op.real and op.imag) into cuAmp*
    cuAmp* elems = qureg.cuStateVec;
    for (long long int i=0; i<op.numElems; i++)
        elems[i] = TO_CU_AMP(op.real[i], ((conj)? -1 : 1) * op.imag[i]);

    std::vector<int> t(targs, targs + op.numQubits); 
    custatevec_applyDiagonal(qureg, {}, t, elems);
}

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) 
{
    thrust::transform(
        getStartDiagonalOpAmpIter(op), getEndDiagonalOpAmpIter(op),
        getStartPtr(qureg), getStartPtr(qureg), // both are begin iters
        thrust::multiplies<cuAmp>());
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op)
{
    auto startIndIter = thrust::make_counting_iterator(0);
    auto endIndIter = startIndIter + qureg.numAmpsTotal;
    auto functor = functor_multDiagOntoDensMatr(qureg.deviceCuStateVec, op);
    thrust::for_each(startIndIter, endIndIter, functor);
}



/*
 * DECOHERENCE
 */

void densmatr_mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg)
{
    Complex facOut   = {.real=(1-otherProb), .imag=0};
    Complex facOther = {.real=otherProb,     .imag=0};
    Complex zero     = {.real=0,             .imag=0};

    statevec_setWeightedQureg(facOther, otherQureg, zero, otherQureg, facOut, combineQureg);
}

void densmatr_mixDephasing(Qureg qureg, int targetQubit, qreal dephase) 
{
    cuAmp a = TO_CU_AMP(1, 0);
    cuAmp b = TO_CU_AMP(1-dephase, 0); // dephase = 2*prob
    cuAmp elems[] = {a, b, b, a};

    std::vector<int> targs{targetQubit, targetQubit + qureg.numQubitsRepresented};
    custatevec_applyDiagonal(qureg, {}, targs, elems);
}

void densmatr_mixTwoQubitDephasing(Qureg qureg, int qb1, int qb2, qreal dephase)
{
    // this function effects the two-qubit dephasing on a density matrix, via the
    // four-qubit diagonal superoperator on a Choi vector. The 16 elements of the
    // diagonal have only two unique entries; 1 and 1-2*dephase/3. It's conceivable
    // that a bespoke kernel could be faster, but likely by little.

    cuAmp a = TO_CU_AMP(1, 0);
    cuAmp b = TO_CU_AMP(1 - dephase, 0); // dephase = 4*prob/3
    cuAmp elems[] = {a, b, b, b,   b, a, b, b,   b, b, a, b,   b, b, b, a};

    int shift = qureg.numQubitsRepresented;
    std::vector<int> targs{qb1, qb2, qb1 + shift, qb2 + shift};
    custatevec_applyDiagonal(qureg, {}, targs, elems);
}

void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depol)
{
    // this function effects depolarising as a dense two-qubit superoperator 
    // on a Choi vector, where only 6 of the 16 elements are non-zero. This is
    // potentially wasteful, and a bespoke kernel could be faster, leveraging
    // QuEST's existing GPU code (or the optimised algorithm in the "distributed"
    // manuscript).

    // depol = (4*prob)/3.0
    cuAmp a = TO_CU_AMP(depol/2., 0);     // 2*prob/3
    cuAmp b = TO_CU_AMP(1 - depol/2., 0); // 1-2*prob/3
    cuAmp c = TO_CU_AMP(1 - depol, 0);    // 1-4*prob/3
    cuAmp z = TO_CU_AMP(0, 0);            // 0

    cuMatr matr{
        b, z, z, a,
        z, c, z, z,
        z, z, c, z,
        a, z, z, b
    };
    std::vector<int> targs{ targetQubit, targetQubit + qureg.numQubitsRepresented };
    custatevec_applyMatrix(qureg, {}, targs, matr);
}

void densmatr_mixDamping(Qureg qureg, int qb, qreal prob)
{
    // this function effects damping as a dense two-qubit superoperator 
    // on a Choi vector, where only 5 of the 16 elements are non-zero. This is
    // potentially wasteful, and a bespoke kernel could be faster, leveraging
    // QuEST's existing GPU code (or the optimised algorithm in the "distributed"
    // manuscript).

    cuAmp w = TO_CU_AMP(1, 0);
    cuAmp z = TO_CU_AMP(0, 0);
    cuAmp p = TO_CU_AMP(prob, 0);
    cuAmp a = TO_CU_AMP(sqrt(1 - prob), 0);
    cuAmp b = TO_CU_AMP(1-prob, 0);

    cuMatr matr{
        w, z, z, p,
        z, a, z, z,
        z, z, a, z,
        z, z, z, b
    };
    std::vector<int> targs{qb, qb + qureg.numQubitsRepresented};
    custatevec_applyMatrix(qureg, {}, targs, matr);
}



/*
 * CALCULATIONS
 */

qreal densmatr_calcTotalProb(Qureg qureg)
{
    long long int diagStride = 1 + (1LL << qureg.numQubitsRepresented);

    cuAmp sum = thrust::reduce(
        getStartStridedAmpIter(qureg, diagStride),
        getEndStridedAmpIter(qureg, diagStride),
        TO_CU_AMP(0, 0));

    return cuAmpReal(sum);
}

qreal statevec_calcTotalProb(Qureg qureg)
{
    qreal abs2sum0;
    qreal abs2sum1;
    int basisBits[] = {0};
    int numBasisBits = 1;

    custatevecAbs2SumOnZBasis(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec,
        &abs2sum0, &abs2sum1, basisBits, numBasisBits);

    return abs2sum0 + abs2sum1;
}

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    // we only ever compute prob(outcome=0) as per the QuEST API's limitation
    qreal prob0;
    qreal* prob1 = nullptr;

    int basisBits[] = {measureQubit};
    int numBasisBits = 1;

    custatevecAbs2SumOnZBasis(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec,
        &prob0, prob1, basisBits, numBasisBits);

    return (outcome == 0)? prob0 : (1-prob0);
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    long long int blockSize = (1LL << measureQubit);
    long long int diagStride = (1LL << qureg.numQubitsRepresented) + 1;

    // we could set this to outcome to sum the outcome=1 amps directly, but the
    // QuEST API specifies that the outcome=0 amps are always summed instead.
    int useOffsetBlocks = 0;

    cuAmp sum = thrust::reduce(
        getStartBlockStridedAmpIter(qureg, blockSize, useOffsetBlocks, diagStride),
        getEndBlockStridedAmpIter(qureg, blockSize, useOffsetBlocks, diagStride),
        TO_CU_AMP(0, 0));

    qreal prob0 = cuAmpReal(sum);
    return (outcome == 0)? prob0 : (1-prob0);
}

qreal densmatr_calcInnerProduct(Qureg a, Qureg b)
{
    cuAmp prod = thrust::inner_product(
        getStartPtr(a), getEndPtr(a), getStartPtr(b), 
        TO_CU_AMP(0,0), thrust::plus<cuAmp>(), functor_prodOfConj());

    return cuAmpReal(prod);
}

Complex statevec_calcInnerProduct(Qureg bra, Qureg ket)
{
    cuAmp prod = thrust::inner_product(
        getStartPtr(bra), getEndPtr(bra), getStartPtr(ket), 
        TO_CU_AMP(0,0), thrust::plus<cuAmp>(), functor_prodOfConj());

    return (Complex) {.real=cuAmpReal(prod), .imag=cuAmpImag(prod)};
}

qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState)
{
    // create iterator f(i) = sum_j qureg_ij pureState_j
    auto functor = functor_matrixColumnDotVector(
        qureg.deviceCuStateVec, pureState.deviceCuStateVec, // both init iters
        pureState.numAmpsTotal);
    auto rawIndIter = thrust::make_counting_iterator(0);
    auto prodAmpIter = thrust::make_transform_iterator(rawIndIter, functor);

    // compute sum_i conj(pureState_i) * f(i)
    cuAmp prod = thrust::inner_product(
        getStartPtr(pureState), getEndPtr(pureState), prodAmpIter, 
        TO_CU_AMP(0,0), thrust::plus<cuAmp>(), functor_prodOfConj());

    qreal prodRe = cuAmpReal(prod);
    return prodRe;
}

qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b)
{
    cuAmp trace = thrust::inner_product(
        getStartPtr(a), getEndPtr(a), getStartPtr(b), 
        TO_CU_AMP(0,0), thrust::plus<cuAmp>(), functor_absOfDif());

    qreal dist = sqrt(cuAmpReal(trace));
    return dist;
}

qreal densmatr_calcPurity(Qureg qureg)
{
    qreal sumOfAbsSquared = thrust::transform_reduce(
        getStartPtr(qureg), getEndPtr(qureg), functor_absSquared(),
        0., thrust::plus<qreal>());

    return sumOfAbsSquared;
}

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    cuAmp prod = thrust::inner_product(
        getStartPtr(qureg), getEndPtr(qureg), getStartDiagonalOpAmpIter(op), 
        TO_CU_AMP(0,0), thrust::plus<cuAmp>(), functor_prodOfAbsSquared());

    return (Complex) {.real=cuAmpReal(prod), .imag=cuAmpImag(prod)};
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
   long long int diagStride = 1 + (1LL << qureg.numQubitsRepresented);

    cuAmp prod = thrust::inner_product(
        getStartStridedAmpIter(qureg, diagStride),
        getEndStridedAmpIter(qureg, diagStride),
        getStartDiagonalOpAmpIter(op),
        TO_CU_AMP(0,0), thrust::plus<cuAmp>(), thrust::multiplies<cuAmp>());

    return (Complex) {.real=cuAmpReal(prod), .imag=cuAmpImag(prod)};
}



/*
 * REDUCTIONS
 */

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{
    int basisBits[] = {measureQubit};

    custatevecCollapseOnZBasis(
        qureg.cuConfig->cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        outcome, basisBits, 1, outcomeProb);
}

void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{
    auto functor = functor_collapseDensMatrToOutcome(
        qureg.deviceCuStateVec, qureg.numQubitsRepresented, 
        outcome, outcomeProb, measureQubit);

    auto startIndIter = thrust::make_counting_iterator(0);
    auto endIndIter = startIndIter + qureg.numAmpsTotal;
    thrust::for_each(startIndIter, endIndIter, functor);
}



#ifdef __cplusplus
}
#endif
