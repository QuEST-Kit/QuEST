// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * An implementation of QuEST's backend (../QuEST_internal.h) using NVIDIA's cuQuantum library.
 * This makes no use of the ComplexArray qureg.deviceStateVec, used by the bespoke GPU kernels,
 * which is not malloc'd in this deployment. Instead, this cuQuantum backend mallocs and uses
 * two dedicated arrays of 'cuAmp' complex primitives; qureg.cuStateVec (CPU memory) and
 * qureg.deviceCuStateVec (GPU memory)
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

// convenient operator overloads for cuAmp, for doing complex artihmetic
cuAmp operator - (const cuAmp& a) {
    return TO_CU_AMP(-cuAmpReal(a), -cuAmpImag(a));
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



#ifdef __cplusplus
extern "C" {
#endif



/*
 * CUQUANTUM WRAPPERS (to reduce boilerplate)
 */

void custatevec_applyMatrix(Qureg qureg, std::vector<int> ctrls, std::vector<int> targs, cuMatr matr) {

    // do not adjoint matrix
    int adj = 0;

    // condition all ctrls on =1 state
    int* ctrlBits = nullptr;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    custatevecApplyMatrix(
        qureg.cuQuantumHandle, 
        qureg.deviceCuStateVec, CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        matr.data(), CU_AMP_IN_MATRIX_PREC, CUSTATEVEC_MATRIX_LAYOUT_ROW, adj, 
        targs.data(), targs.size(), 
        ctrls.data(), ctrlBits, ctrls.size(), 
        CUSTATEVEC_COMPUTE_DEFAULT,
        work, workSize);
}



/* 
 * ENVIRONMENT MANAGEMENT
 */

int GPUSupportsMemPools() {

    // consult only the first device (garuanteed already to exist)
    int device = 0;

    int supports;
    cudaDeviceGetAttribute(&supports, cudaDevAttrMemoryPoolsSupported, device);
    return supports;
}

int memPoolAlloc(void* ctx, void** ptr, size_t size, cudaStream_t stream) {
    cudaMemPool_t pool = *static_cast<cudaMemPool_t*>(ctx);
    return cudaMallocFromPoolAsync(ptr, size, pool, stream); 
}
int memPoolFree(void* ctx, void* ptr, size_t size, cudaStream_t stream) {
    return cudaFreeAsync(ptr, stream); 
}

void setupAutoWorkspaces(custatevecHandle_t cuQuantumHandle) {

    // get the current (device's default) stream-ordered memory pool (assuming single GPU)
    int deviceId;
    cudaGetDevice(&deviceId);
    cudaMemPool_t memPool;
    cudaDeviceGetMemPool(&memPool, deviceId);

    // get its current memory threshold, above which memory gets freed at every stream synch
    size_t currMaxMem;
    cudaMemPoolGetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &currMaxMem); 

    // if it's smaller than 1 MiB = 16 qubits, extend it
    size_t desiredMaxMem = 16*(1<<15);
    if (currMaxMem < desiredMaxMem)
        cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &desiredMaxMem); 

    // create a mem handler around the mem pool
    custatevecDeviceMemHandler_t memHandler;
    memHandler.ctx = &memPool;
    memHandler.device_alloc = memPoolAlloc;
    memHandler.device_free = memPoolFree;
    strcpy(memHandler.name, "mempool");
    printf("created mem handler\n");

    // set cuQuantum to use this handler and pool, to automate workspace memory management
    custatevecSetDeviceMemHandler(cuQuantumHandle, &memHandler);
}

QuESTEnv createQuESTEnv(void) {
    validateGPUExists(GPUExists(), __func__);
    validateGPUIsCuQuantumCompatible(GPUSupportsMemPools(),__func__);
    
    QuESTEnv env;
    env.rank=0;
    env.numRanks=1;
    
    env.seeds = NULL;
    env.numSeeds = 0;
    seedQuESTDefault(&env);

    // prepare cuQuantum
    custatevecCreate(&env.cuQuantumHandle);
    cudaStreamCreate(&env.cuStream);    
    custatevecSetStream(env.cuQuantumHandle, env.cuStream);
    setupAutoWorkspaces(env.cuQuantumHandle);
    
    return env;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);

    // finalise cuQuantum
    custatevecDestroy(env.cuQuantumHandle);
    cudaStreamDestroy(env.cuStream);
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

    // copy env's cuQuantum handles
    qureg->cuQuantumHandle = env.cuQuantumHandle;
    qureg->cuStream = env.cuStream;

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

void densmatr_initPureState(Qureg targetQureg, Qureg copyQureg)
{
}

void densmatr_initPlusState(Qureg qureg)
{
}

void densmatr_initClassicalState(Qureg qureg, long long int stateInd)
{
}

void statevec_initBlankState(Qureg qureg)
{
}

void statevec_initZeroState(Qureg qureg)
{
}

void statevec_initPlusState(Qureg qureg)
{
}

void statevec_initClassicalState(Qureg qureg, long long int stateInd)
{
}

void statevec_initDebugState(Qureg qureg)
{
}

void statevec_initStateOfSingleQubit(Qureg *qureg, int qubitId, int outcome)
{
}

int statevec_initStateFromSingleFile(Qureg *qureg, char filename[200], QuESTEnv env)
{
    return -1;
}

void densmatr_setQuregToPauliHamil(Qureg qureg, PauliHamil hamil)
{
}

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out)
{
}



/*
 * DEBUG
 */

void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank)
{
}

int statevec_compareStates(Qureg mq1, Qureg mq2, qreal precision)
{
    return -1;
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
        qureg.cuQuantumHandle, 
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
    // this diagonal operator, otherwise embarrasingly parallel with unit stride, 
    // is here treated as a generic one-target unitary, wastefully inducing non-unit 
    // stride and unnecessary memory reads, and potentially unnecessary communication
    // in multi-GPU mode. 

    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuAmp aE = toCuAmp(term);
    cuMatr matrix{
        a1, a0,
        a0, aE
    };
    custatevec_applyMatrix(qureg, {}, {targetQubit}, matrix);
}

void statevec_controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle)
{
    // this diagonal operator, otherwise embarrasingly parallel with unit stride, 
    // is here treated as a generic one-target unitary, wastefully inducing non-unit 
    // stride and unnecessary memory reads, and potentially unnecessary communication
    // in multi-GPU mode. 

    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuAmp aE = TO_CU_AMP(cos(angle), sin(angle));
    cuMatr matrix{
        a1, a0,
        a0, aE
    };
    custatevec_applyMatrix(qureg, {idQubit1}, {idQubit2}, matrix);
}

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle)
{   
    // this diagonal operator, otherwise embarrasingly parallel with unit stride, 
    // is here treated as a generic one-target unitary, wastefully inducing non-unit 
    // stride and unnecessary memory reads, and potentially unnecessary communication
    // in multi-GPU mode. 

    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuAmp aE = TO_CU_AMP(cos(angle), sin(angle));
    cuMatr matrix{
        a1, a0,
        a0, aE
    };
    std::vector<int> targs{controlQubits[0]};
    std::vector<int> ctrls(controlQubits + 1, controlQubits + numControlQubits); 
    custatevec_applyMatrix(qureg, ctrls, targs, matrix);
}

void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle)
{   
    qreal theta = - angle/2.;
    std::vector<int> targs = getIndsFromMask(mask, qureg.numQubitsInStateVec);
    std::vector<custatevecPauli_t> paulis(targs.size(), CUSTATEVEC_PAULI_Z);

    custatevecApplyPauliRotation(
        qureg.cuQuantumHandle, qureg.deviceCuStateVec, 
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
        qureg.cuQuantumHandle, qureg.deviceCuStateVec, 
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec, 
        theta, paulis.data(), targs.data(), targs.size(),
        ctrls.data(), nullptr, ctrls.size());
}

void statevec_controlledPhaseFlip(Qureg qureg, int idQubit1, int idQubit2)
{
    // this diagonal operator, otherwise embarrasingly parallel with unit stride, 
    // is here treated as a generic one-target unitary, wastefully inducing non-unit 
    // stride and unnecessary memory reads, and potentially unnecessary communication
    // in multi-GPU mode. 

    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuMatr matrix{
        a1,  a0,
        a0, -a1
    };
    custatevec_applyMatrix(qureg, {idQubit1}, {idQubit2}, matrix);
}

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits)
{
    // this diagonal operator, otherwise embarrasingly parallel with unit stride, 
    // is here treated as a generic one-target unitary, wastefully inducing non-unit 
    // stride and unnecessary memory reads, and potentially unnecessary communication
    // in multi-GPU mode. 

    cuAmp a0 = TO_CU_AMP(0, 0);
    cuAmp a1 = TO_CU_AMP(1, 0);
    cuMatr matrix{
        a1,  a0,
        a0, -a1
    };
    std::vector<int> targs{controlQubits[0]};
    std::vector<int> ctrls(controlQubits + 1, controlQubits + numControlQubits); 
    custatevec_applyMatrix(qureg, ctrls, targs, matrix);
}

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) 
{
    int2 targPairs[] = {{qb1, qb2}}; 
    int numPairs = 1;

    custatevecSwapIndexBits(
        qureg.cuQuantumHandle, qureg.deviceCuStateVec, 
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
    // isntead resort to slowly (by at most a factor #targs) effect it as a sequence
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

void statevec_applySubDiagonalOp(Qureg qureg, int* targets, SubDiagonalOp op, int conj)
{
    // sneakily leverage the CPU cuQuantum memory in order to convert op
    // (as separate arrays op.real and op.imag) into cuAmp*
    cuAmp* diagonals = qureg.cuStateVec;
    for (long long int i=0; i<op.numElems; i++)
        diagonals[i] = TO_CU_AMP(op.real[i], op.imag[i]);

    custatevecApplyGeneralizedPermutationMatrix(
        qureg.cuQuantumHandle, qureg.deviceCuStateVec,
        CU_AMP_IN_STATE_PREC, qureg.numQubitsInStateVec,
        nullptr, diagonals, CU_AMP_IN_MATRIX_PREC, 0, 
        targets, op.numQubits, nullptr, nullptr, 0,
        nullptr, 0);
}

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) 
{
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op)
{
}

void statevec_applyPhaseFuncOverrides(
    Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int numTerms, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
}

void statevec_applyMultiVarPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int* numTermsPerReg, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
}

void statevec_applyParamNamedPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    enum phaseFunc phaseFuncName, qreal* params, int numParams,
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
}



/*
 * DECOHERENCE
 */

void densmatr_mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg)
{
}

void densmatr_mixDephasing(Qureg qureg, int targetQubit, qreal dephase) 
{
}

void densmatr_mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal dephase)
{
}

void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depolLevel)
{
}

void densmatr_mixDamping(Qureg qureg, int targetQubit, qreal damping)
{
}

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel)
{
}



/*
 * CALCULATIONS
 */

qreal densmatr_calcTotalProb(Qureg qureg)
{
    return -1;
}

qreal statevec_calcTotalProb(Qureg qureg)
{
    return -1;
}

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    return -1;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    return -1;
}

void statevec_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits)
{
}

void densmatr_calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits)
{
}

qreal densmatr_calcInnerProduct(Qureg a, Qureg b)
{
    return -1;
}

Complex statevec_calcInnerProduct(Qureg bra, Qureg ket)
{
    return (Complex) {.real=-1, .imag=-1};
}

qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState)
{
    return -1;
}

qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b)
{
    return -1;
}

qreal densmatr_calcPurity(Qureg qureg)
{
    return -1;
}

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    return (Complex) {.real=-1, .imag=-1};
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op)
{
    return (Complex) {.real=-1, .imag=-1};
}



/*
 * REDUCTIONS
 */

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{        
}

void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal outcomeProb)
{
}



#ifdef __cplusplus
}
#endif
