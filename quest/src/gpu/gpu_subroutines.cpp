/** @file
 * CUDA GPU-accelerated definitions of the subroutines called by
 * accelerator.cpp. This file is always compiled, even when GPU
 * acceleration is disabled and when parsed by a non-CUDA compiler,
 * and so uses precompiler guards to disable CUDA-only code. 
 * This file contains host definitions and associated memory and
 * thread management, and invokes custom kernels defined in
 * gpu_kernels.hpp which is never parsed by non-CUDA compilers. 
 * This file also invokes Thrust and cuQuantum routines, defined in
 * gpu_thrust.hpp and gpu_cuquantum.hpp respectively, which are
 * also never parsed by non-CUDA compilers.
 * 
 * Note that some custom kernels are templated in order to apply
 * compile-time optimisations like automatic loop unrolling. So
 * too are their calling host definitions in this file, which are
 * called by accelerator.cpp which chooses the template parameter.
 * This unnecessarily duplicates other parts of the host functions
 * responsible for dispatching to Thrust or cuQuantum, bloating the
 * compiled binary size; we accept this design wart over having
 * non-templated host functions because this requires duplicating
 * the template-dispatch logic (which would then also be defined in
 * cpu_subroutines.cpp) and moving it out of the aptly-named
 * accelerator.cpp file.
 * 
 * Despite COMPILE_CUDA=1 whenever COMPILE_CUQUANTUM=1, we will
 * still use superfluous (COMPILE_CUDA || COMPILE_CUQUANTUM) guards
 * to communicate when there is no bespoke cuQuantum routine.
 *
 * When compiling for AMD GPUs, the CUDA symbols invoked herein are
 * mapped to HIP symbols by cuda_to_hip.h 
 * 
 * @author Tyson Jones
 */

#if (COMPILE_CUQUANTUM && ! COMPILE_CUDA)
    #error "Cannot define COMPILE_CUQUANTUM=1 without simultaneously defining COMPILE_CUDA=1"
#endif

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/comm/comm_indices.hpp"
#include "quest/src/gpu/gpu_config.hpp"
#include "quest/src/gpu/gpu_subroutines.hpp"

#if COMPILE_CUDA
    #include "quest/src/gpu/gpu_types.cuh"
    #include "quest/src/gpu/gpu_kernels.cuh"
    #include "quest/src/gpu/gpu_thrust.cuh"
#endif

#if COMPILE_CUQUANTUM
    #include "quest/src/gpu/gpu_cuquantum.cuh"
#endif

#include <vector>
using std::vector;



/*
 * GETTERS
 */


qcomp gpu_statevec_getAmp_sub(Qureg qureg, qindex ind) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    // this bespoke function exists (in lieu of caller 
    // just calling copyGpuToCpu() directly) mostly for
    // consistency with the CPU equivalent (which has a
    // performance motivation), and so that we can one
    // day update this function if single-scalar random
    // access of GPU memory can be accelerated.

    qcomp amp;

    // compiler guards harmlessly duplicated therein
    gpu_copyGpuToCpu(&qureg.gpuAmps[ind], &amp, 1);

    return amp;

#else
    error_gpuCopyButGpuNotCompiled();
    return -1;
#endif
}



/*
 * SETTERS
 */


void gpu_densmatr_setAmpsToPauliStrSum_sub(Qureg qureg, PauliStrSum sum) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_densmatr_setAmpsToPauliStrSum_sub(qureg, sum);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_fullstatediagmatr_setElemsToPauliStrSum(out, in);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * COMMUNICATION BUFFER PACKING
 */


template <int NumQubits>
qindex gpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates) {

    assert_numQubitsMatchesQubitStatesAndTemplateParam(qubits.size(), qubitStates.size(), NumQubits);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(qubits.size());
    qindex numBlocks = getNumBlocks(numThreads);
    qindex sendInd = getSubBufferSendInd(qureg);

    devints sortedQubits = util_getSorted(qubits);
    qindex qubitStateMask  = util_getBitMask(qubits, qubitStates);

    kernel_statevec_packAmpsIntoBuffer <NumQubits> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[sendInd], numThreads, 
        getPtr(sortedQubits), qubits.size(), qubitStateMask
    );

    // return the number of packed amps
    return numThreads;

#else
    error_gpuSimButGpuNotCompiled();
    return 0;
#endif
}


qindex gpu_statevec_packPairSummedAmpsIntoBuffer(Qureg qureg, int qubit1, int qubit2, int qubit3, int bit2) {

    assert_bufferPackerGivenIncreasingQubits(qubit1, qubit2, qubit3);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 8;
    qindex numBlocks = getNumBlocks(numThreads);
    qindex sendInd = getSubBufferSendInd(qureg);

    kernel_statevec_packPairSummedAmpsIntoBuffer <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[sendInd], numThreads, 
        qubit1, qubit2, qubit3, bit2
    );

    // return the number of packed amps
    return numThreads;

#else
    error_gpuSimButGpuNotCompiled();
    return 0;
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qindex, gpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )



/*
 * SWAPS
 */


template <int NumCtrls> 
void gpu_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUQUANTUM

    cuquantum_statevec_anyCtrlSwap_subA(qureg, ctrls, ctrlStates, targ1, targ2);

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(2 + ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devints sortedQubits = util_getSorted(ctrls, {targ2, targ1});
    qindex qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ2, targ1}, {0, 1});

    kernel_statevec_anyCtrlSwap_subA <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, 
        getPtr(sortedQubits), ctrls.size(), qubitStateMask, targ1, targ2
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <int NumCtrls> 
void gpu_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    devints sortedCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    kernel_statevec_anyCtrlSwap_subB <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        getPtr(sortedCtrls), ctrls.size(), ctrlStateMask
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <int NumCtrls> 
void gpu_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(1 + ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    devints sortedQubits = util_getSorted(ctrls, {targ});
    qindex qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ}, {targState});

    kernel_statevec_anyCtrlSwap_subC <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        getPtr(sortedQubits), ctrls.size(), qubitStateMask
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlSwap_subA, (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlSwap_subB, (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlSwap_subC, (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) )



/*
 * ONE-TARGET DENSE MATRIX
 */


template <int NumCtrls>
void gpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUQUANTUM

    auto arr = unpackMatrixToCuQcomps(matr);
    cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(qureg, ctrls, ctrlStates, {targ}, arr.data());

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size() + 1);
    qindex numBlocks = getNumBlocks(numThreads);

    devints sortedQubits = util_getSorted(ctrls, {targ});
    qindex qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ}, {0});

    auto [m00, m01, m10, m11] = unpackMatrixToCuQcomps(matr);

    kernel_statevec_anyCtrlOneTargDenseMatr_subA <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, 
        getPtr(sortedQubits), ctrls.size(), qubitStateMask, targ, 
        m00, m01, m10, m11
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <int NumCtrls>
void gpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    devints sortedCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    kernel_statevec_anyCtrlOneTargDenseMatr_subB <NumCtrls> <<<numBlocks,NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        getPtr(sortedCtrls), ctrls.size(), ctrlStateMask, 
        toCuQcomp(fac0), toCuQcomp(fac1)
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * TWO-TARGET DENSE MATRIX
 */


template <int NumCtrls> 
void gpu_statevec_anyCtrlTwoTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUQUANTUM

    auto arr = unpackMatrixToCuQcomps(matr);
    cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(qureg, ctrls, ctrlStates, {targ1, targ2}, arr.data());

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size() + 2);
    qindex numBlocks = getNumBlocks(numThreads);

    devints sortedQubits = util_getSorted(ctrls, {targ1,targ2});
    qindex qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ1,targ2}, {0,0});

    // unpack matrix elems which are more efficiently accessed by kernels as args than shared mem (... maybe...)
    auto m = unpackMatrixToCuQcomps(matr);

    kernel_statevec_anyCtrlTwoTargDenseMatr_sub <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, 
        getPtr(sortedQubits), ctrls.size(), qubitStateMask, targ1, targ2,
        m[0], m[1], m[2],  m[3],  m[4],  m[5],  m[6],  m[7],
        m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]
    );
 
#else
    error_gpuSimButGpuNotCompiled();
#endif
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlTwoTargDenseMatr_sub, (Qureg, vector<int>, vector<int>, int, int, CompMatr2) )



/*
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj>
void gpu_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

#if COMPILE_CUQUANTUM

    auto matrElemsPtr = toCuQcomps(matr.gpuElemsFlat);
    auto matrElemsLen = matr.numRows * matr.numRows;

    // conjugate every matrix element if necessary (cuStateVec cannot conj for us; only adjoint)
    if (ApplyConj)
        thrust_setElemsToConjugate(matrElemsPtr, matrElemsLen);

    cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(qureg, ctrls, ctrlStates, targs, matrElemsPtr);

    // undo conjugation (which is only not done if cuQuantum encounters a non-recoverable internal error)
    if (ApplyConj)
        thrust_setElemsToConjugate(matrElemsPtr, matrElemsLen);

#elif COMPILE_CUDA

    // a 'batch' refers to 2^N amps which become mixed by the matrix,
    // distinguished in this kernel from 'numThreads' since we may
    // task each thread with processing more than a single batch
    qindex numBatches = qureg.numAmpsPerNode / powerOf2(ctrls.size() + targs.size());

    devints deviceTargs = targs;
    devints deviceQubits = util_getSorted(ctrls, targs);
    qindex qubitStateMask = util_getBitMask(ctrls, ctrlStates, targs, vector<int>(targs.size(),0));

    // unpacking args (to better distinguish below signatures)
    auto ampsPtr   = toCuQcomps(qureg.gpuAmps);
    auto matrPtr   = toCuQcomps(matr.gpuElemsFlat);
    auto qubitsPtr = getPtr(deviceQubits);
    auto targsPtr  = getPtr(deviceTargs);
    auto nCtrls    = ctrls.size();

    // this function updates amplitudes in batches of 2^NumTargs, where each is
    // determined by distinct mixtures of the existing 2^NumTargs values, which
    // must ergo be cached. As such, each thread needs private memory, which is
    // provided either by fast registers (when NumTargs is compile-time known)
    // or by slow global memory (necessary when blockDim * 2^NumTargs exceeds the
    // total shared memory) which is strided for coalesced access by warp threads
    
    if constexpr (NumTargs != -1) {

        // when NumTargs <= 5, each thread has a private array stored in the registers,
        // enabling rapid IO. Given NUM_THREADS_PER_BLOCK = 128, the maximum size of 
        // this array per-block is 16 * 128 * 2^5 B = 64 KiB which exceeds shared
        // memory capacity, but does NOT exceed maximum register capacity.

        /// @todo
        /// We should really check the above claims, otherwise the thread-private arrays could
        /// silently "spill" from registers into "local memory" (which is really slow,
        /// global memory) and greatly sabotage performance on some GPUs.

        qindex numThreads = numBatches;
        qindex numBlocks = getNumBlocks(numThreads);

        kernel_statevec_anyCtrlFewTargDenseMatr<NumCtrls, NumTargs, ApplyConj> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
            ampsPtr, numThreads, 
            qubitsPtr, nCtrls, qubitStateMask, 
            targsPtr, matrPtr
        );

    } else {

        // when NumTargs > 5, we must use global memory to give each thread a private
        // workspace of 2^NumTargs elements. Alas, we should not simply allocate this 
        // space per-thread, since all threads being potentially concurrent means we 
        // would allocate a total cache equal to the Qureg size (when nctrls=0)! 
        // Instead, we change the parallelisation granularity, giving each thread more
        // batches of 2^NumTargs amps to modify, re-using its private cache, so that the
        // number of potentially-concurrent threads is reduced, as is the total cache.
        // We choose the granularity by upper-bounding the number of concurrent threads,
        // where we assign one-block per multiprocessor because we are anyway memory-
        // bandwidth bound (so we don't expect many interweaved blocks per MP).
        qindex numThreads = gpu_getMaxNumConcurrentThreads();
        
        // use strictly 2^# threads to maintain precondition of all kernels
        if (!isPowerOf2(numThreads))
            numThreads = util_getNextPowerOf2(numThreads);

        // no point in dispatching more threads than batches
        if (numThreads > numBatches)
            numThreads = numBatches;

        // evenly distribute the batches between threads, and the threads unevenly between blocks
        qindex numBatchesPerThread = numBatches / numThreads; // divides evenly
        qindex numBlocks = getNumBlocks(numThreads);

        // expand the cache if necessary
        qindex numKernelInvocations = numBlocks * NUM_THREADS_PER_BLOCK;
        qcomp* cache = gpu_getCacheOfSize(powerOf2(targs.size()), numKernelInvocations);

        kernel_statevec_anyCtrlManyTargDenseMatr <NumCtrls, ApplyConj> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
            toCuQcomps(cache),
            ampsPtr, numThreads, numBatchesPerThread, 
            qubitsPtr, nCtrls, qubitStateMask, 
            targsPtr, targs.size(), powerOf2(targs.size()), matrPtr
        );
    }

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, gpu_statevec_anyCtrlAnyTargDenseMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, CompMatr) )



/*
 * ONE-TARGET DIAG MATRIX
 */


template <int NumCtrls> 
void gpu_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // diagonal matrices are always embarrassingly parallel, but alas
    // cuQuantum's API cannot handle when the targeted qubits are prefix,
    // since they appear larger than the local statevector (unknown by it
    // to be distributed). cuQuantum would require we gratuitously swap
    // such qubits into the suffix qubits before invocation, which is unnece-
    // ssary. As such, we will only invoke cuQuantum when the targeted qubits
    // (in this function, only one) are within the suffix substate, otherwise
    // we fall back to using our custom kernels which never require comm.

#if COMPILE_CUQUANTUM

    if (util_isQubitInSuffix(targ, qureg)) {

        // we never conjugate DiagMatr1 at this level; the caller will have already conjugated
        bool conj = false;

        // we can pass 1D CPU .elems array directly to cuQuantum which will recognise host pointers
        cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, ctrls, ctrlStates, {targ}, toCuQcomps(matr.elems), conj);
        
        // explicitly return to avoid re-simulation below
        return;
    }

#endif

// note preprocessors are not exclusive
#if COMPILE_CUDA

    /// @todo
    /// when NumCtrls==0, a Thrust functor would be undoubtedly more
    /// efficient (because of improved parallelisation granularity) 

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devints deviceCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    auto elems = unpackMatrixToCuQcomps(matr);

    kernel_statevec_anyCtrlOneTargDiagMatr_sub <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qureg.rank, qureg.logNumAmpsPerNode,
        getPtr(deviceCtrls), ctrls.size(), ctrlStateMask, targ, elems[0], elems[1]
    );

    // explicitly return to avoid runtime error below
    return;

#endif

    // only reachable when nothing above simulated
    error_gpuSimButGpuNotCompiled();
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlOneTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, int, DiagMatr1) )



/*
 * TWO-TARGET DIAG MATRIX
 */


template <int NumCtrls> 
void gpu_statevec_anyCtrlTwoTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // diagonal matrices are always embarrassingly parallel, but alas
    // cuQuantum's API cannot handle when the targeted qubits are prefix,
    // since they appear larger than the local statevector (unknown by it
    // to be distributed). cuQuantum would require we gratuitously swap
    // such qubits into the suffix qubits before invocation, which is unnece-
    // ssary. As such, we will only invoke cuQuantum when the targeted qubits
    // are both within the suffix substate, otherwise we fall back to using 
    // our custom kernels which never require comm.

#if COMPILE_CUQUANTUM

    if (util_areAllQubitsInSuffix({targ1,targ2}, qureg)) {

        // we never conjugate DiagMatr2 at this level; the caller will have already conjugated
        bool conj = false;

        // we can pass 1D CPU array directly to cuQuantum, and it will recognise host pointers
        cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, ctrls, ctrlStates, {targ1, targ2}, toCuQcomps(matr.elems), conj);

        // explicitly return to avoid re-simulation below
        return;
    }

#endif 

// note preprocessors are not exclusive
#if COMPILE_CUDA

    /// @todo
    /// when NumCtrls==0, a Thrust functor would be undoubtedly more
    /// efficient (because of improved parallelisation granularity) 

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devints deviceCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    auto elems = unpackMatrixToCuQcomps(matr);

    kernel_statevec_anyCtrlTwoTargDiagMatr_sub <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qureg.rank, qureg.logNumAmpsPerNode,
        getPtr(deviceCtrls), ctrls.size(), ctrlStateMask, targ1, targ2,
        elems[0], elems[1], elems[2], elems[3]
    );

    // explicitly return to avoid runtime error below
    return;

#endif

    // only reachable when nothing above simulated
    error_gpuSimButGpuNotCompiled();
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlTwoTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, int, int, DiagMatr2) )



/*
 * ANY-TARGET DIAG MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj, bool HasPower>
void gpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);
    assert_exponentMatchesTemplateParam(exponent, HasPower);

    // diagonal matrices are always embarrassingly parallel, but alas
    // cuQuantum's API cannot handle when the targeted qubits are prefix,
    // since they appear larger than the local statevector (unknown by it
    // to be distributed). cuQuantum would require we gratuitously swap
    // such qubits into the suffix qubits before invocation, which is unnece-
    // ssary. As such, we will only invoke cuQuantum when the targeted qubits
    // are both within the suffix substate, otherwise we fall back to using 
    // our custom kernels which never require comm. Furthermore, cuQuantum
    // cannot handle when exponent != 1, for which we also fallback to custom.

#if COMPILE_CUQUANTUM

    // cuQuantum cannot handle HasPower, in which case we fall back to custom kernel
    if (!HasPower && util_areAllQubitsInSuffix(targs, qureg)) {
        cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, ctrls, ctrlStates, targs, toCuQcomps(util_getGpuMemPtr(matr)), ApplyConj);

        // must return to avoid re-simulation below
        return;
    }

#endif

// note preprocessors are not exclusive
#if COMPILE_CUDA

    /// @todo
    /// when NumCtrls==0, a Thrust functor would be undoubtedly more
    /// efficient (because of improved parallelisation granularity) 

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devints deviceTargs = targs;
    devints deviceCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    kernel_statevec_anyCtrlAnyTargDiagMatr_sub <NumCtrls, NumTargs, ApplyConj, HasPower> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qureg.rank, qureg.logNumAmpsPerNode,
        getPtr(deviceCtrls), ctrls.size(), ctrlStateMask, getPtr(deviceTargs), targs.size(), 
        toCuQcomps(util_getGpuMemPtr(matr)), toCuQcomp(exponent)
    );

    // must return to avoid runtime error below
    return;

#endif

    // only reachable when nothing above simulated
    error_gpuSimButGpuNotCompiled();
}


INSTANTIATE_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, gpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr, qcomp) )



/*
 * ALL-TARGS DIAGONAL MATRIX
 */


template <bool HasPower>
void gpu_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    assert_exponentMatchesTemplateParam(exponent, HasPower);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    // we always use Thrust because we are doubtful that cuQuantum's
    // diagonal-matrix facilities are optimised for the all-qubit case

    thrust_statevec_allTargDiagMatr_sub<HasPower>(qureg, matr, toCuQcomp(exponent));

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <bool HasPower, bool MultiplyOnly>
void gpu_densmatr_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    assert_exponentMatchesTemplateParam(exponent, HasPower);

    // in theory, we could use cuQuantum when HasPower=MultiplyOnly=true,
    // treating FullStateDiagMatr like an N/2-qubit DiagMatr upon a SV,
    // but this scenario is not worth the code complication

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    kernel_densmatr_allTargDiagMatr_sub <HasPower, MultiplyOnly> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qureg.rank, qureg.logNumAmpsPerNode,
        toCuQcomps(util_getGpuMemPtr(matr)), matr.numElems, toCuQcomp(exponent)
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template void gpu_statevec_allTargDiagMatr_sub<true >(Qureg, FullStateDiagMatr, qcomp);
template void gpu_statevec_allTargDiagMatr_sub<false>(Qureg, FullStateDiagMatr, qcomp);

template void gpu_densmatr_allTargDiagMatr_sub<true, true>  (Qureg, FullStateDiagMatr, qcomp);
template void gpu_densmatr_allTargDiagMatr_sub<true, false> (Qureg, FullStateDiagMatr, qcomp);
template void gpu_densmatr_allTargDiagMatr_sub<false, true> (Qureg, FullStateDiagMatr, qcomp);
template void gpu_densmatr_allTargDiagMatr_sub<false, false>(Qureg, FullStateDiagMatr, qcomp);



/*
 * PAULI TENSOR AND GADGET
 */


template <int NumCtrls, int NumTargs> 
void gpu_statevector_anyCtrlPauliTensorOrGadget_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> x, vector<int> y, vector<int> z, qcomp ampFac, qcomp pairAmpFac) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(x.size() + y.size(), NumTargs);

    // we do not make use of cuQuantum's custatevecApplyGeneralizedPermutationMatrix() to effect
    // a pauli tensor because we wish to avoid creating the (2^#paulis) large permutation matrix.
    // We also cannot make use of cuQuantum's custatevecApplyPauliRotation() because it cannot
    // handle Pauli operators upon the prefix substate as our singly-communicating method does.
    // This is true even if we passed down the gadget phase to this function; cuStateVec would
    // exact amp -> a amp + b other_amp for the wrong b, which we cannot thereafter remedy.

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qcomp powI   = util_getPowerOfI(y.size());
    auto targsXY = util_getConcatenated(x, y);
    auto maskXY  = util_getBitMask(targsXY);
    auto maskYZ  = util_getBitMask(util_getConcatenated(y, z));

    devints deviceTargs   = targsXY;
    devints deviceQubits  = util_getSorted(ctrls, targsXY);
    qindex qubitStateMask = util_getBitMask(ctrls, ctrlStates, targsXY, vector<int>(targsXY.size(),0));

    // unlike the analogous cpu routine, this function has only a single parallelisation
    // granularity; where every pair-of-amps is modified by an independent thread, despite
    // that many threads share a common i0 value (appearing in the kernel). This turns out
    // faster than when giving threads many pair-amps to modify, due to memory movements

    qindex numThreads = (qureg.numAmpsPerNode / powerOf2(ctrls.size())) / 2; // divides evenly
    qindex numBlocks = getNumBlocks(numThreads);
    kernel_statevector_anyCtrlPauliTensorOrGadget_subA <NumCtrls, NumTargs> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        getPtr(deviceQubits), ctrls.size(), qubitStateMask, 
        getPtr(deviceTargs), deviceTargs.size(),
        maskXY, maskYZ, toCuQcomp(powI), toCuQcomp(ampFac), toCuQcomp(pairAmpFac)
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <int NumCtrls> 
void gpu_statevector_anyCtrlPauliTensorOrGadget_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> x, vector<int> y, vector<int> z, qcomp ampFac, qcomp pairAmpFac, qindex bufferMaskXY) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    qcomp powI = util_getPowerOfI(y.size());
    auto maskXY = util_getBitMask(util_getConcatenated(x, y));
    auto maskYZ = util_getBitMask(util_getConcatenated(y, z));

    devints sortedCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    kernel_statevector_anyCtrlPauliTensorOrGadget_subB <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        getPtr(sortedCtrls), ctrls.size(), ctrlStateMask,
        maskXY, maskYZ, bufferMaskXY,
        toCuQcomp(powI), toCuQcomp(ampFac), toCuQcomp(pairAmpFac)
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, gpu_statevector_anyCtrlPauliTensorOrGadget_subA, (Qureg, vector<int>, vector<int>, vector<int>, vector<int>, vector<int>, qcomp, qcomp) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevector_anyCtrlPauliTensorOrGadget_subB, (Qureg, vector<int>, vector<int>, vector<int>, vector<int>, vector<int>, qcomp, qcomp, qindex) )



/*
 * PHASE TENSOR AND GADGET
 */


template <int NumCtrls> 
void gpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qcomp fac0, qcomp fac1) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devints sortedCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    qindex targMask = util_getBitMask(targs);

    kernel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub <NumCtrls> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        getPtr(sortedCtrls), ctrls.size(), ctrlStateMask, targMask,
        toCuQcomp(fac0), toCuQcomp(fac1)
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub, (Qureg, vector<int>, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * QUREG COMBINATION
 */


void gpu_statevec_setQuregToSuperposition_sub(qcomp facOut, Qureg outQureg, qcomp fac1, Qureg inQureg1, qcomp fac2, Qureg inQureg2) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_statevec_setQuregToSuperposition_sub(toCuQcomp(facOut), outQureg, toCuQcomp(fac1), inQureg1, toCuQcomp(fac2), inQureg2);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inQureg) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_densmatr_mixQureg_subA(outProb, outQureg, inProb, inQureg);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_mixQureg_subB(qreal outProb, Qureg outQureg, qreal inProb, Qureg inQureg) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = outQureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    kernel_densmatr_mixQureg_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        outProb, toCuQcomps(outQureg.gpuAmps), inProb, toCuQcomps(inQureg.gpuAmps),
        numThreads, inQureg.numAmps
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_mixQureg_subC(qreal outProb, Qureg outQureg, qreal inProb) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = outQureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    kernel_densmatr_mixQureg_subC <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        outProb, toCuQcomps(outQureg.gpuAmps), inProb, toCuQcomps(outQureg.gpuCommBuffer),
        numThreads, outQureg.rank, powerOf2(outQureg.numQubits), outQureg.logNumAmpsPerNode        
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * ONE-QUBIT DEPHASING
 */


void gpu_densmatr_oneQubitDephasing_subA(Qureg qureg, int ketQubit, qreal prob) {

#if COMPILE_CUQUANTUM

    // gauranteed that corresponding braQubit is in suffix, so always safe to call cuQuantum
    cuquantum_densmatr_oneQubitDephasing_subA(qureg, ketQubit, prob);

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / 4;
    qindex numBlocks = getNumBlocks(numThreads);

    auto fac = util_getOneQubitDephasingFactor(prob);
    int braQubit = util_getBraQubit(ketQubit, qureg);

    kernel_densmatr_oneQubitDephasing_subA <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, ketQubit, braQubit, fac
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_oneQubitDephasing_subB(Qureg qureg, int ketQubit, qreal prob) {

#if COMPILE_CUQUANTUM 

    // gauranteed that corresponding braQubit is in prefix; however, cuQuantum effects
    // the gate as a phase*Id gate on any qubit, so just picks one in suffix
    cuquantum_densmatr_oneQubitDephasing_subB(qureg, ketQubit, prob);

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numThreads);

    auto fac = util_getOneQubitDephasingFactor(prob);
    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);

    kernel_densmatr_oneQubitDephasing_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, ketQubit, braBit, fac
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * TWO-QUBIT DEPHASING
 */


void gpu_densmatr_twoQubitDephasing_subA(Qureg qureg, int ketQubitA, int ketQubitB, qreal prob) {

#if COMPILE_CUQUANTUM

    // gauranteed that both corresponding braQubits are in prefix, so safe to invoke cuQuantum
    cuquantum_densmatr_twoQubitDephasing_subA(qureg, ketQubitA, ketQubitB, prob);

#elif COMPILE_CUDA

    // the rank-agnostic version is identical to the subB algorithm below, because the
    // queried bits of the global index i below will always be in the suffix substate.
    gpu_densmatr_twoQubitDephasing_subB(qureg, ketQubitA, ketQubitB, prob);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_twoQubitDephasing_subB(Qureg qureg, int ketQubitA, int ketQubitB, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM 

    qindex numThreads = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    auto term = util_getTwoQubitDephasingTerm(prob);
    int braQubitA = util_getBraQubit(ketQubitA, qureg);
    int braQubitB = util_getBraQubit(ketQubitB, qureg);

    kernel_densmatr_twoQubitDephasing_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qureg.rank, qureg.logNumAmpsPerNode, // numAmps, not numCols
        ketQubitA, ketQubitB, braQubitA, braQubitB, term
    );
    
#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * ONE-QUBIT DEPOLARISING
 */


void gpu_densmatr_oneQubitDepolarising_subA(Qureg qureg, int ketQubit, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 4;
    qindex numBlocks = getNumBlocks(numThreads);

    int braQubit = util_getBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitDepolarisingFactors(prob);

    kernel_densmatr_oneQubitDepolarising_subA <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, ketQubit, braQubit, factors.c1, factors.c2, factors.c3
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_oneQubitDepolarising_subB(Qureg qureg, int ketQubit, qreal prob) {
    
#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitDepolarisingFactors(prob);

    kernel_densmatr_oneQubitDepolarising_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        ketQubit, braBit, factors.c1, factors.c2, factors.c3
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * TWO-QUBIT DEPOLARISING
 */


void gpu_densmatr_twoQubitDepolarising_subA(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braQb2 = util_getBraQubit(ketQb2, qureg);
    auto c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    kernel_densmatr_twoQubitDepolarising_subA <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        ketQb1, ketQb2, braQb1, braQb2, c3
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_twoQubitDepolarising_subB(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 16;
    qindex numBlocks = getNumBlocks(numThreads);

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braQb2 = util_getBraQubit(ketQb2, qureg);
    auto factors = util_getTwoQubitDepolarisingFactors(prob);

    // each kernel invocation sums all 4 amps together, so adjusts c1
    qreal altc1 = factors.c1 - factors.c2;

    kernel_densmatr_twoQubitDepolarising_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        ketQb1, ketQb2, braQb1, braQb2, altc1, factors.c2
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_twoQubitDepolarising_subC(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);
    auto c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    kernel_densmatr_twoQubitDepolarising_subC <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        ketQb1, ketQb2, braQb1, braBit2, c3
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_twoQubitDepolarising_subD(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 8;
    qindex numBlocks = getNumBlocks(numThreads);
    qindex offset = getBufferRecvInd();

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);
    auto factors = util_getTwoQubitDepolarisingFactors(prob);

    kernel_densmatr_twoQubitDepolarising_subD <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[offset], numThreads,
        ketQb1, ketQb2, braQb1, braBit2, factors.c1, factors.c2
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_twoQubitDepolarising_subE(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    qreal fac0 = 1 + factors.c3;
    qreal fac1 = factors.c1 - fac0;

    kernel_densmatr_twoQubitDepolarising_subE <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        ketQb1, ketQb2, braBit1, braBit2, fac0, fac1
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_twoQubitDepolarising_subF(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 4;
    qindex numBlocks = getNumBlocks(numThreads);
    qindex offset = getBufferRecvInd();

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);
    auto c2 = util_getTwoQubitDepolarisingFactors(prob).c2;

    kernel_densmatr_twoQubitDepolarising_subF <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[offset], numThreads,
        ketQb1, ketQb2, braBit1, braBit2, c2
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * PAULI CHANNEL
 */


void gpu_densmatr_oneQubitPauliChannel_subA(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 4;
    qindex numBlocks = getNumBlocks(numThreads);

    int braQubit = util_getBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitPauliChannelFactors(pI, pX, pY, pZ);

    kernel_densmatr_oneQubitPauliChannel_subA <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, ketQubit, braQubit, 
        factors.c1, factors.c2, factors.c3, factors.c4
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_oneQubitPauliChannel_subB(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {
    
#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitPauliChannelFactors(pI, pX, pY, pZ);

    kernel_densmatr_oneQubitPauliChannel_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        ketQubit, braBit, factors.c1, factors.c2, factors.c3, factors.c4
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * AMPLITUDE DAMPING
 */


void gpu_densmatr_oneQubitDamping_subA(Qureg qureg, int ketQubit, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 4;
    qindex numBlocks = getNumBlocks(numThreads);

    int braQubit = util_getBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitDampingFactors(prob);

    kernel_densmatr_oneQubitDamping_subA <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads,
        ketQubit, braQubit, prob, factors.c1, factors.c2
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_oneQubitDamping_subB(Qureg qureg, int qubit, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numThreads);

    auto c2 = util_getOneQubitDampingFactors(prob).c2;

    kernel_densmatr_oneQubitDamping_subB <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qubit, c2
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_oneQubitDamping_subC(Qureg qureg, int ketQubit, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numThreads);

    auto braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto c1 = util_getOneQubitDampingFactors(prob).c1;

    kernel_densmatr_oneQubitDamping_subC <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, ketQubit, braBit, c1
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_densmatr_oneQubitDamping_subD(Qureg qureg, int qubit, qreal prob) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numThreads);
    qindex recvInd = getBufferRecvInd();

    kernel_densmatr_oneQubitDamping_subD <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[recvInd], numThreads, 
        qubit, prob
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 * PARTIAL TRACE
 */


template <int NumTargs> 
void gpu_densmatr_partialTrace_sub(Qureg inQureg, Qureg outQureg, vector<int> targs, vector<int> pairTargs) {

    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = outQureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    devints devTargs = targs;
    devints devPairTargs = pairTargs;
    devints devAllTargs = util_getSorted(targs, pairTargs);

    kernel_densmatr_partialTrace_sub <NumTargs> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(inQureg.gpuAmps), toCuQcomps(outQureg.gpuAmps), numThreads,
        getPtr(devTargs), getPtr(devPairTargs), getPtr(devAllTargs), targs.size()
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, gpu_densmatr_partialTrace_sub, (Qureg, Qureg, vector<int>, vector<int>) )



/*
 * PROBABILITIES
 */


qreal gpu_statevec_calcTotalProb_sub(Qureg qureg) {

#if COMPILE_CUQUANTUM
    return cuquantum_statevec_calcTotalProb_sub(qureg);

#elif COMPILE_CUDA
    return thrust_statevec_calcTotalProb_sub(qureg);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


qreal gpu_densmatr_calcTotalProb_sub(Qureg qureg) {

#if COMPILE_CUQUANTUM || COMPILE_CUDA
    return thrust_densmatr_calcTotalProb_sub(qureg);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template <int NumQubits> 
qreal gpu_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

#if COMPILE_CUQUANTUM

    // cuQuantum disregards NumQubits compile-time param
    return cuquantum_statevec_calcProbOfMultiQubitOutcome_sub(qureg, qubits, outcomes);

#elif COMPILE_CUDA 

    return thrust_statevec_calcProbOfMultiQubitOutcome_sub<NumQubits>(qureg, qubits, outcomes);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template <int NumQubits> 
qreal gpu_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

#if COMPILE_CUQUANTUM || COMPILE_CUDA 

    return thrust_densmatr_calcProbOfMultiQubitOutcome_sub<NumQubits>(qureg, qubits, outcomes);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template <int NumQubits> 
void gpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

#if COMPILE_CUQUANTUM

    /// @todo
    /// cuQuantum assumes all qubits are local (since it does not consult rank) 
    /// and so cannot be used when any qubits are in prefix, despite that the 
    /// algorithm is always embarrassingly parallel. In theory we can workaround
    /// this by only passing suffix-qubits to the backend, and shuffling the
    /// resulting outProbs array within localiser.cpp. Experiment with this!

    if (util_areAllQubitsInSuffix(qubits, qureg)) {

        // cuQuantum discards NumQubits template param and creates
        // its own temporary GPU memory (to store probs) if necessary
        cuquantum_statevec_calcProbsOfAllMultiQubitOutcomes_sub(outProbs, qureg, qubits);

        // explicitly return to avoid re-simulation below
        return;
    }

#endif

// note preprocessors are not exclusive
#if COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numThreads);

    // allocate exponentially-big temporary memory (error if failed)
    devints devQubits = qubits;
    devreals devProbs = getDeviceRealsVec(powerOf2(qubits.size())); // throws

    kernel_statevec_calcProbsOfAllMultiQubitOutcomes_sub<NumQubits> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        getPtr(devProbs), toCuQcomps(qureg.gpuAmps), numThreads, 
        qureg.rank, qureg.logNumAmpsPerNode, getPtr(devQubits), devQubits.size()
    );

    // overwrite outProbs with GPU memory
    copyFromDeviceVec(devProbs, outProbs);

    // explicitly return to avoid error-msg below
    return;

#endif

    // should be unreadchable
    error_gpuSimButGpuNotCompiled();
}


template <int NumQubits> 
void gpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    // we decouple numColsPerNode and numThreads for clarity
    // (and in case parallelisation granularity ever changes);
    qindex numThreads = powerOf2(qureg.logNumColsPerNode);
    qindex numBlocks = getNumBlocks(numThreads);
    
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);

    // allocate exponentially-big temporary memory (error if failed)
    devints devQubits = qubits;
    devreals devProbs = getDeviceRealsVec(powerOf2(qubits.size())); // throws

    kernel_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<NumQubits> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        getPtr(devProbs), toCuQcomps(qureg.gpuAmps), 
        numThreads, firstDiagInd, numAmpsPerCol,
        qureg.rank, qureg.logNumAmpsPerNode, 
        getPtr(devQubits), devQubits.size()
    );

    // overwrite outProbs with GPU memory
    copyFromDeviceVec(devProbs, outProbs);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qreal, gpu_statevec_calcProbOfMultiQubitOutcome_sub, (Qureg, vector<int>, vector<int>) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qreal, gpu_densmatr_calcProbOfMultiQubitOutcome_sub, (Qureg, vector<int>, vector<int>) )

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, gpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub, (qreal* outProbs, Qureg, vector<int>) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, gpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub, (qreal* outProbs, Qureg, vector<int>) )



/*
 * INNER PRODUCTS
 */


qcomp gpu_statevec_calcInnerProduct_sub(Qureg quregA, Qureg quregB) {

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp prod = thrust_statevec_calcInnerProduct_sub(quregA, quregB);
    return toQcomp(prod);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


qreal gpu_densmatr_calcHilbertSchmidtDistance_sub(Qureg quregA, Qureg quregB) {

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    return thrust_densmatr_calcHilbertSchmidtDistance_sub(quregA, quregB);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template <bool Conj>
qcomp gpu_densmatr_calcFidelityWithPureState_sub(Qureg rho, Qureg psi) {

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp fid = thrust_densmatr_calcFidelityWithPureState_sub<Conj>(rho, psi);
    return toQcomp(fid);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template qcomp gpu_densmatr_calcFidelityWithPureState_sub<true >(Qureg, Qureg);
template qcomp gpu_densmatr_calcFidelityWithPureState_sub<false>(Qureg, Qureg);



/*
 * PAULI EXPECTATION VALUES
 */


qreal gpu_statevec_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

#if COMPILE_CUQUANTUM

    return cuquantum_statevec_calcExpecAnyTargZ_sub(qureg, targs);

#elif COMPILE_CUDA

    return thrust_statevec_calcExpecAnyTargZ_sub(qureg, targs);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


qcomp gpu_densmatr_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp value = thrust_densmatr_calcExpecAnyTargZ_sub(qureg, targs);
    return toQcomp(value);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


qcomp gpu_statevec_calcExpecPauliStr_subA(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

#if COMPILE_CUQUANTUM

    return cuquantum_statevec_calcExpecPauliStr_subA(qureg, x, y, z);

#elif COMPILE_CUDA

    cu_qcomp value = thrust_statevec_calcExpecPauliStr_subA(qureg, x, y, z);
    return toQcomp(value);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


qcomp gpu_statevec_calcExpecPauliStr_subB(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp value = thrust_statevec_calcExpecPauliStr_subB(qureg, x, y, z);
    return toQcomp(value);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


qcomp gpu_densmatr_calcExpecPauliStr_sub(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {
    
#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp value = thrust_densmatr_calcExpecPauliStr_sub(qureg, x, y, z);
    return toQcomp(value);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}



/*
 * DIAGONAL MATRIX EXPECTATION VALUES
 */


template <bool HasPower, bool UseRealPow> 
qcomp gpu_statevec_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {
    assert_exponentMatchesTemplateParam(exponent, HasPower, UseRealPow);

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp expo = toCuQcomp(exponent);
    cu_qcomp value = thrust_statevec_calcExpecFullStateDiagMatr_sub<HasPower,UseRealPow>(qureg, matr, expo);
    return toQcomp(value);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template <bool HasPower, bool UseRealPow>
qcomp gpu_densmatr_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {
    assert_exponentMatchesTemplateParam(exponent, HasPower, UseRealPow);

#if COMPILE_CUQUANTUM || COMPILE_CUDA

    cu_qcomp expo = toCuQcomp(exponent);
    cu_qcomp value = thrust_densmatr_calcExpecFullStateDiagMatr_sub<HasPower,UseRealPow>(qureg, matr, expo);
    return toQcomp(value);

#else
    error_gpuSimButGpuNotCompiled();
    return -1;
#endif
}


template qcomp gpu_statevec_calcExpecFullStateDiagMatr_sub<true, true >(Qureg, FullStateDiagMatr, qcomp);
template qcomp gpu_statevec_calcExpecFullStateDiagMatr_sub<true, false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp gpu_statevec_calcExpecFullStateDiagMatr_sub<false,false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp gpu_statevec_calcExpecFullStateDiagMatr_sub<false,true >(Qureg, FullStateDiagMatr, qcomp); // uncallable

template qcomp gpu_densmatr_calcExpecFullStateDiagMatr_sub<true, true >(Qureg, FullStateDiagMatr, qcomp);
template qcomp gpu_densmatr_calcExpecFullStateDiagMatr_sub<true, false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp gpu_densmatr_calcExpecFullStateDiagMatr_sub<false,false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp gpu_densmatr_calcExpecFullStateDiagMatr_sub<false,true >(Qureg, FullStateDiagMatr, qcomp); // uncallable


/*
 * PROJECTORS
 */


template <int NumQubits> 
void gpu_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {

    // all qubits are in suffix
    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

#if COMPILE_CUQUANTUM

    // cuQuantum disregards NumQubits template param
    cuquantum_statevec_multiQubitProjector_sub(qureg, qubits, outcomes, prob);

#elif COMPILE_CUDA

    qreal renorm = 1 / std::sqrt(prob);
    thrust_statevec_multiQubitProjector_sub<NumQubits>(qureg, qubits, outcomes, renorm);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <int NumQubits> 
void gpu_densmatr_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {

    // qubits are unconstrained, and can include prefix qubits
    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qreal renorm = 1 / prob;
    thrust_densmatr_multiQubitProjector_sub<NumQubits>(qureg, qubits, outcomes, renorm);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, gpu_statevec_multiQubitProjector_sub, (Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, gpu_densmatr_multiQubitProjector_sub, (Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) )



/*
 * STATE INITIALISATION
 */


void gpu_statevec_initUniformState_sub(Qureg qureg, qcomp amp) {
#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_statevec_initUniformState(qureg, toCuQcomp(amp));

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_statevec_initDebugState_sub(Qureg qureg) {
#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_statevec_initDebugState_sub(qureg);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


void gpu_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg qureg) {

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    thrust_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(qureg);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}
