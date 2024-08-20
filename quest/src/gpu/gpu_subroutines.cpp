/** @file
 * CUDA GPU-accelerated definitions of the subroutines called by
 * accelerator.cpp. This file is always compiled, even when GPU
 * acceleration is disabled and when parsed by a non-CUDA compiler,
 * and so uses precompiler guards to disable CUDA-only code. Note
 * that COMPILE_CUDA=1 whenever COMPILE_CUQUANTUM=1, but we may
 * still use superfluous (COMPILE_CUDA || COMPILE_CUQUANTUM) guards
 * to communicate when there is no bespoke cuQuantum routine.
 * 
 * This file contains host definitions and associated memory and
 * thread management, and invokes custom kernels defined in
 * kernels.hpp which is never parsed by non-CUDA compilers. This
 * file also invokes Thrust and cuQuantum routines, defined in
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
 * When compiling for AMD GPUs, the CUDA symbols invoked herein are
 * mapped to HIP symbols by cuda_to_hip.h 
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/comm/comm_indices.hpp"

#if COMPILE_CUDA
    #include "quest/src/gpu/gpu_types.hpp"
    #include "quest/src/gpu/gpu_kernels.hpp"
    #include "quest/src/gpu/gpu_thrust.hpp"
#endif

#if COMPILE_CUQUANTUM
    #include "quest/src/gpu/gpu_cuquantum.hpp"
#endif

#include <vector>

using std::vector;



/*
 * COMMUNICATION BUFFER PACKING
 */


template <int NumQubits>
void gpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates) {

    assert_numQubitsMatchesQubitStatesAndTemplateParam(qubits.size(), qubitStates.size(), NumQubits);

#if COMPILE_CUDA || COMPILE_CUQUANTUM

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(qubits.size());
    qindex numBlocks = getNumBlocks(numThreads);
    qindex sendInd = getSubBufferSendInd(qureg);

    devicevec sortedQubits = util_getSorted(qubits);
    qindex qubitStateMask  = util_getBitMask(qubits, qubitStates);

    kernel_statevec_packAmpsIntoBuffer <NumQubits> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), &toCuQcomps(qureg.gpuCommBuffer)[sendInd], numThreads, 
        getPtr(sortedQubits), qubits.size(), qubitStateMask
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )



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

    devicevec sortedQubits = util_getSorted(ctrls, {targ2, targ1});
    qindex qubitStateMask  = util_getBitMask(ctrls, ctrlStates, {targ2, targ1}, {0, 1});

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

    devicevec sortedCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask  = util_getBitMask(ctrls, ctrlStates);

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

    devicevec sortedQubits = util_getSorted(ctrls, {targ});
    qindex qubitStateMask  = util_getBitMask(ctrls, ctrlStates, {targ}, {targState});

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

    devicevec sortedQubits = util_getSorted(ctrls, {targ});
    qindex qubitStateMask  = util_getBitMask(ctrls, ctrlStates, {targ}, {0});

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

    devicevec sortedCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask  = util_getBitMask(ctrls, ctrlStates);

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
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs>
void gpu_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

#if COMPILE_CUQUANTUM

    cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(qureg, ctrls, ctrlStates, targs, toCuQcomps(matr.gpuElems));

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size() + targs.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devicevec deviceTargs  = targs;
    devicevec deviceQubits = util_getSorted(ctrls, targs);
    qindex qubitStateMask  = util_getBitMask(ctrls, ctrlStates, targs, vector<int>(targs.size(),0));

    // TODO: allocating per-thread private memory here is ridiculously slow and pointless, and makes redundant
    // the persistent GPU memory in the CompMatr. Revise this to use persistent env-attached memory, as per 
    // discussions with Richard Meister.

    cu_qcomp *cache;
    qindex numTargAmps = powerOf2(targs.size());
    qindex gridSize = numBlocks * NUM_THREADS_PER_BLOCK;
    qindex cacheSize = numTargAmps * gridSize;
    cudaMalloc(&cache, cacheSize * sizeof *cache);

    kernel_statevec_anyCtrlAnyTargDenseMatr_sub <NumCtrls, NumTargs> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), cache, numThreads,
        getPtr(deviceQubits), ctrls.size(), qubitStateMask, getPtr(deviceTargs), targs.size(),
        toCuQcomps(matr.gpuElems));

    cudaFree(cache);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, gpu_statevec_anyCtrlAnyTargDenseMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, CompMatr) )



/*
 * DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs>
void gpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

#if COMPILE_CUQUANTUM

    cuquantum_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, ctrls, ctrlStates, targs, toCuQcomps(matr.gpuElems));

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThreads);

    devicevec deviceTargs = targs;
    devicevec deviceCtrls = util_getSorted(ctrls);
    qindex ctrlStateMask  = util_getBitMask(ctrls, ctrlStates);

    kernel_statevec_anyCtrlAnyTargDiagMatr_sub <NumCtrls, NumTargs> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThreads, qureg.rank, qureg.logNumAmpsPerNode,
        getPtr(deviceCtrls), ctrls.size(), ctrlStateMask, getPtr(deviceTargs), targs.size(), 
        toCuQcomps(matr.gpuElems)
    );

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, gpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr) )
