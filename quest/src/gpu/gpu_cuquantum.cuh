/** @file
 * Subroutines which invoke cuStateVec, which are alternatives to the
 * kernels defined in gpu_kernels.cuh, as invoked by gpu_subroutines.cpp
 * 
 * This file is only ever included when COMPILE_CUQUANTUM=1 and COMPILE_CUDA=1
 * so it can safely invoke CUDA signatures without guards. Note that many of 
 * the statevector functions herein will be re-leveraged by QuEST's density
 * matrix simulation, so it important we do not pass Qureg.numQubits to the 
 * cuStateVec API, and instead pass qureg.logNumAmpsPerNode. 
 * 
 * This file is a (CUDA) header since only ever included by gpu_subroutines.cpp.
 * 
 * @author Tyson Jones
 */


// because this file uses a global instance of CuQuantumConfig (not inlined),
// it must never be included by multiple translation units; it can only ever
// be included by gpu_subroutines.cpp
#ifdef GPU_CUQUANTUM_HPP
    #error "File gpu_cuquantum.hpp was erroneously included by multiple source files."
#endif

#ifndef GPU_CUQUANTUM_HPP
#define GPU_CUQUANTUM_HPP


// check preprocessors and compilers are valid before #includes to avoid compile errors

#if ! COMPILE_CUQUANTUM
    #error "A file being compiled somehow included gpu_cuquantum.hpp despite QuEST not being compiled in cuQuantum mode."
#endif

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_cuquantum.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#ifndef __NVCC__
    #error "A file which included gpu_cuquantum.hpp was attemptedly compiled with a non-CUDA compiler."
#endif


#include "quest/include/precision.h"

#include "quest/src/core/utilities.hpp"
#include "quest/src/gpu/gpu_config.hpp"
#include "quest/src/gpu/gpu_types.cuh"

#include <custatevec.h>
#include <vector>
#include <new>

using std::vector;



/*
 * PRECISION FLAG
 *
 * which we will use for both the state precision AND gate matrix precision,
 * because QuEST uses only a single qcomp type for both in the API.
 */

#if (FLOAT_PRECISION == 1)
    #define CUQUANTUM_QCOMP CUDA_C_32F

#elif (FLOAT_PRECISION == 2)
    #define CUQUANTUM_QCOMP CUDA_C_64F

#else
    #error "Build bug; precision.h should have prevented non-float non-double qcomp precision on GPU (and cuQuantum)."
    
#endif



/*
 * ENVIRONMENT MANAGEMENT
 */


// sets the size at which the CUDA memory pool will 
// automatically deallocate temporary memory. Below this
// size, temporary memory structures (like a CompMatr)
// will persist in GPU memory to save time. This is
// only relevant to GPU-mode with cuQuantum enabled,
// and is effected at createQuESTEnv().
size_t CUQUANTUM_MEM_POOL_BYTES = 16*(1<<16); // 1 MiB ~ 8 qubit complex<double> matrix


struct CuQuantumConfig {
    cudaMemPool_t mempool;
    cudaStream_t stream;
    custatevecHandle_t handle;
    custatevecDeviceMemHandler_t memhandler;
};

// singleton handle to cuQuantum env needed for applying gates and finalizing env
CuQuantumConfig config;


cudaMemPool_t getExistingMemPool() {

    // validation gaurantees memory pool already exists
    cudaMemPool_t memPool;
    int deviceId;
    CUDA_CHECK( cudaGetDevice(&deviceId) );
    CUDA_CHECK( cudaDeviceGetMemPool(&memPool, deviceId) );
    return memPool;
}


void adjustMemPoolSize(cudaMemPool_t memPool) {

    // find existing memPool threshold (above which memory gets freed at every stream synch)
    size_t currMaxMem;
    CUDA_CHECK( cudaMemPoolGetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &currMaxMem) ); 

    // optionally increase memPool threshold
    if (currMaxMem < CUQUANTUM_MEM_POOL_BYTES)
        CUDA_CHECK( cudaMemPoolSetAttribute(memPool, cudaMemPoolAttrReleaseThreshold, &CUQUANTUM_MEM_POOL_BYTES) ); 
}

int allocMemInPool(void* ctx, void** ptr, size_t size, cudaStream_t stream) {
    cudaMemPool_t& pool = *static_cast<cudaMemPool_t*>(ctx);
    return cudaMallocFromPoolAsync(ptr, size, pool, stream); 
}

int deallocMemInPool(void* ctx, void* ptr, size_t size, cudaStream_t stream) {
    return cudaFreeAsync(ptr, stream); 
}


void gpu_initCuQuantum() {

    // the cuStateVec docs say custatevecCreate() should be called
    // once per physical GPU, though oversubscribing MPI processes
    // while setting PERMIT_NODES_TO_SHARE_GPU=1 worked fine in our
    // testing - we will treat it as tolerable but undefined behaviour

    // create new stream and cuQuantum handle, binding to global config
    CUDA_CHECK( custatevecCreate(&config.handle) );
    CUDA_CHECK( cudaStreamCreate(&config.stream) );

    // get and configure existing memory pool (for later automatic alloc/dealloc of gate matrices)
    config.mempool = getExistingMemPool();
    adjustMemPoolSize(config.mempool);

    // create a temporary cuQuantum memory handler
    config.memhandler.ctx = &config.mempool;
    config.memhandler.device_alloc = allocMemInPool;
    config.memhandler.device_free = deallocMemInPool;
    strcpy(config.memhandler.name, "mempool");

    // bind memory handler and stream to cuQuantum handle
    CUDA_CHECK( custatevecSetDeviceMemHandler(config.handle, &config.memhandler) );
    CUDA_CHECK( custatevecSetStream(config.handle, config.stream) );
}


void gpu_finalizeCuQuantum() {

    CUDA_CHECK( cudaStreamDestroy(config.stream) );
    CUDA_CHECK( custatevecDestroy(config.handle) );
}



/*
 * GATES
 */


void cuquantum_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    // our SWAP targets are bundled into pairs
    int2 targPairs[] = {{targ1, targ2}};;
    int numTargPairs = 1;

    CUDA_CHECK( custatevecSwapIndexBits(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        targPairs, numTargPairs,

        // swap mask params seem to be in the reverse order to the remainder of the cuStateVec API
        ctrlStates.data(), ctrls.data(), ctrls.size()) );
}


// there is no bespoke cuquantum_statevec_anyCtrlSwap_subB()



/*
 * MATRICES
 */


void cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, cu_qcomp* flatMatrElems) {

    // this funciton is called 'subA' instead of just 'sub', because it is also called in 
    // the one-target case whereby it is strictly the embarrassingly parallel _subA scenario

    // do not adjoint matrix
    int adj = 0;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    CUDA_CHECK( custatevecApplyMatrix(
        config.handle, 
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode, 
        flatMatrElems, CUQUANTUM_QCOMP, CUSTATEVEC_MATRIX_LAYOUT_ROW, adj, 
        targs.data(), targs.size(),
        ctrls.data(), ctrlStates.data(), ctrls.size(), 
        CUSTATEVEC_COMPUTE_DEFAULT,
        work, workSize) );
}


// there is no bespoke cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subB()


void cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, cu_qcomp* flatMatrElems, bool conj) {

    // beware that despite diagonal matrices being embarrassingly parallel,
    // the target qubits must still all be suffix-only to avoid a cuStateVec error

    // apply no permutation matrix
    custatevecIndex_t *perm = nullptr;

    // effect conjugate by adjointing, which is equivalent for a diagonal matrix
    int adj = conj;

    // use automatic workspace management
    void* work = nullptr;
    size_t workSize = 0;

    CUDA_CHECK( custatevecApplyGeneralizedPermutationMatrix(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        perm, flatMatrElems, CUQUANTUM_QCOMP, adj, 
        targs.data(), targs.size(), 
        ctrls.data(), ctrlStates.data(), ctrls.size(),
        work, workSize) );
}



/*
 * DECOHERENCE
 *
 * which is mostly unsupported by the cuStateVec backend, although we repurpose
 * the templated diagonal-matrix functions to effect diagonal superoperators
 */


void cuquantum_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob) {

    // effect the superoperator as a two-qubit diagonal upon a statevector suffix state
    cu_qcomp a = {1,        0};
    cu_qcomp b = {1-2*prob, 0};
    cu_qcomp elems[] = {a, b, b, a};
    vector<int> targs {qubit, qubit + qureg.numQubits};

    bool conj = false;
    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, {}, {}, targs, elems, conj);
}


void cuquantum_densmatr_oneQubitDephasing_subB(Qureg qureg, int ketQubit, qreal prob) {

    // we need to merely scale every amp where ketQubit differs from braBit, which is
    // equivalent to a state-controlled global phase upon a statevector, which is 
    // itself a same-element one-qubit diagonal applied to any target
    int braBit = getBit(qureg.rank, ketQubit - qureg.logNumColsPerNode);
    cu_qcomp fac = {1 - 2*prob, 0};
    cu_qcomp elems[] = {fac, fac};

    // we choose to target the largest possible qubit, expecting best cuStateVec performance;
    // note it must still be a suffix qubit since cuQuantum does not know qureg is distributed
    int targ = qureg.logNumAmpsPerNode - 1; // leftmost suffix bra qubit

    bool conj = false;
    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, {ketQubit}, {!braBit}, {targ}, elems, conj);
}


void cuquantum_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubitA, int qubitB, qreal prob) {

    /// @todo
    /// only 75% of the amps are changed, each of which is multiplied by the same scalar,
    /// but our below method multiplies all amps with 16 separate scalars - can we accel?
    /// we are applying a prefactor b to all states except where the ket & bra qubits
    /// are the same, i.e. we skip |00><00|, |01><01|, |10><10|, |11><11|

    // effect the superoperator as a four-qubit diagonal upon a statevector suffix state
    cu_qcomp a = {1,          0};
    cu_qcomp b = {1-4*prob/3, 0};
    cu_qcomp elems[] = {a,b,b,b, b,a,b,b, b,b,a,b, b,b,b,a};
    vector<int> targs {qubitA, qubitB, qubitA + qureg.numQubits, qubitB + qureg.numQubits};

    bool conj = false;
    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, {}, {}, targs, elems, conj);
}


// there is no bespoke cuquantum_densmatr_twoQubitDephasing_subB()



/*
 * PROBABILITIES
 */


qreal cuquantum_statevec_calcTotalProb_sub(Qureg qureg) {

    // cuQuantum probabilities are always double (not qreal)
    double prob0;
    double prob1;

    // we can find the probability via any qubit, though we target
    // the leftmost so that the reduction is contiguous which we
    // expect has the best performance; beware though that cuQuantum
    // does not know the state is distributed, so we must use the
    // leftmost suffix qubit
    int qubit = qureg.logNumAmpsPerNode - 1;
    int numQubits = 1;

    CUDA_CHECK( custatevecAbs2SumOnZBasis(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        &prob0, &prob1, &qubit, numQubits ) );

    qreal total = prob0 + prob1;
    return total;
}


qreal cuquantum_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    // cuQuantum probabilities are always double
    double prob;

    CUDA_CHECK( custatevecAbs2SumArray(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        &prob, nullptr, 0, outcomes.data(), qubits.data(), qubits.size()) );

    return static_cast<qreal>(prob);
}


void cuquantum_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    // cuQuantum can accept a host-pointer (like outProbs), but only
    // double precision; if qreal != double, we use temporary memory
    #if (FLOAT_PRECISION == 2)
        double* outPtr = outProbs;
    #else
        vector<double> tmpProbs;
        util_tryAllocVector(tmpProbs, powerOf2(qubits.size()), error_cuQuantumTempCpuAllocFailed);
        double* outPtr = tmpProbs.data();
    #endif

    CUDA_CHECK( custatevecAbs2SumArray(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        outPtr, qubits.data(), qubits.size(), nullptr, nullptr, 0) );

    // serially cast and copy output probabilities, if necessary
    #if (FLOAT_PRECISION != 2)
        for (size_t i=0; i<tmpProbs.size(); i++)
            outProbs[i] = static_cast<qreal>(tmpProbs[i]);
    #endif
}



/*
 * EXPECTATION VALUES
 */


qreal cuquantum_statevec_calcExpecPauliStr_subA(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    // prepare term (XX...YY...ZZ...)
    size_t numPaulis = x.size() + y.size() + z.size();
    vector<custatevecPauli_t> paulis; 
    vector<int32_t> targs; 
    
    paulis.reserve(numPaulis);
    targs.reserve(numPaulis);

    for (int t : x) { paulis.push_back(CUSTATEVEC_PAULI_X); targs.push_back(t); }
    for (int t : y) { paulis.push_back(CUSTATEVEC_PAULI_Y); targs.push_back(t); }
    for (int t : z) { paulis.push_back(CUSTATEVEC_PAULI_Z); targs.push_back(t); }

    // prepare terms = {term}
    const custatevecPauli_t* termPaulis[] = {paulis.data()};
    const int32_t* termTargets[] = {targs.data()};
    const uint32_t numPaulisPerTerm[] = { (uint32_t) paulis.size()};
    uint32_t numTerms = 1;

    // cuStateVec output is always double
    double value = 0;

    CUDA_CHECK( custatevecComputeExpectationsOnPauliBasis(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        &value, termPaulis, numTerms, termTargets, numPaulisPerTerm) );

    return static_cast<qreal>(value);
}


qreal cuquantum_statevec_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

    return cuquantum_statevec_calcExpecPauliStr_subA(qureg, {}, {}, targs);
}



/*
 * PROJECTORS
 */


void cuquantum_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {

    CUDA_CHECK( custatevecCollapseByBitString(
        config.handle,
        toCuQcomps(qureg.gpuAmps), CUQUANTUM_QCOMP, qureg.logNumAmpsPerNode,
        outcomes.data(), qubits.data(), qubits.size(), prob) );
}


#endif // GPU_CUQUANTUM_HPP