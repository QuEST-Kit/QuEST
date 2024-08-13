/** @file
 * CUDA GPU-accelerated definitions of the subroutines called by
 * accelerator.cpp. This file contains the host definitions and
 * associated memory and thread management, and invocations of
 * Thrust and cuQuantum subroutines (if the latter is enabled).
 * All custom kernels are defined in kernels.hpp, which is never
 * parsed by non-CUDA compilers and ignored when not compiling GPU.
 * When compiling for AMD GPUs, the CUDA symbols invoked herein are
 * mapped to HIP symbols by cuda_to_hip.h
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/indexer.hpp"

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



// TODO: hmm

#include "quest/src/core/accelerator.hpp"






/*
 * BUFFER PACKING
 *
 * which are templated and require explicit instantiation below
 */


template <int NumCtrls>
void gpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {
    
    // there is no cuQuantum facility for packing

#if COMPILE_CUDA


    // TODO:
    //      I think we have to device-malloc space for ctrls.
    //      or like, copy it into a single compile-time array



    auto [qubits, stateMask] = getSortedQubitsAndMask(ctrls, ctrlStates, {}, {});
    qindex numThreads = qureg.numAmpsPerNode / powerOf2(NumCtrls);

    // cast qcomp to cu_qcomp
    qindex numBlocks = getNumBlocks(numThreads);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);
    cu_qcomp* buff = toCuQcomps(qureg.gpuCommBuffer);

    kernel_statevec_packAmpsIntoBuffer 
        <NumCtrls> 
        <<<numBlocks, NUM_THREADS_PER_BLOCK>>> 
        (amps, buff, qubits.data(), stateMask, numThreads);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(void, gpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>));






/*
 * TODO SECTION
 */


template <int NumCtrls>
void gpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
#if COMPILE_CUQUANTUM

    // ignore ctrlFlag
    cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(qureg, ctrls, ctrlStates, {targ}, unpackMatrixToCuQcomps(matr).data());

#elif COMPILE_CUDA

    // prepare parameters needed for optimal indexing logic with the given ctrls

    qindex numBlocks = getNumBlocks(params.numInds);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);

    // there are sufficiently few matrix elements (<256 bytes) to pass directly to kernel
    cu_qcomp m00, m01, m10, m11;
    unpackMatrixToCuQcomps(matr, m00,m01,m10,m11);

    // use ctrlFlag to dispatch to optimised kernel
    kernel_statevec_anyCtrlOneTargDenseMatr_subA <ctrlFlag> <<<numBlocks,NUM_THREADS_PER_BLOCK>>> (amps, params, targ, m00,m01,m10,m11);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(void, gpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1));



template <int NumCtrls>
void gpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {
    
    // there is no cuQuantum facility for subB

#if COMPILE_CUDA

    // prepare parameters needed for optimal indexing logic with the given ctrls
    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(qureg.numAmpsPerNode, ctrls, ctrlStates);

    // cast qcomp to cu_qcomp
    qindex numBlocks = getNumBlocks(params.numInds);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);
    cu_qcomp* buff = toCuQcomps(qureg.gpuCommBuffer);
    cu_qcomp f0 = toCuQcomp(fac0);
    cu_qcomp f1 = toCuQcomp(fac1);

    // use ctrlFlag to dispatch to optimised kernel
    kernel_statevec_anyCtrlOneTargDenseMatr_subB <ctrlFlag> <<<numBlocks,NUM_THREADS_PER_BLOCK>>> (amps, buff, params, f0, f1);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(void, gpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp));







// template <int NumCtrls>
// void gpu_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {
// #if COMPILE_CUQUANTUM

//     // ignore ctrlFlag
//     cuquantum_statevec_anyCtrlAnyTargDiagMatr_subA(qureg, ctrls, ctrlStates, {targ}, toCuQcomps(matr.elems));

// #elif COMPILE_CUDA

//     qindex numThreads = qureg.numAmpsPerNode / powerOf2(NumCtrls);
//     qindex numBlocks = getNumBlocks(numThreads);
//     cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);
//     cu_qcomp d0 = toCuQcomp(matr.elems[0]);
//     cu_qcomp d1 = toCuQcomp(matr.elems[1]);

//     // TODO: fix kernel args etc

//     // use ctrlFlag to dispatch to optimised kernel
//     kernel_statevec_anyCtrlOneTargDiagMatr_sub <ctrlFlag> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (amps, params, targ, d0, d1);

// #else
//     error_gpuSimButGpuNotCompiled();
// #endif
// }


// template <int NumCtrls, int NumTargs>
// void gpu_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {

//     // TODO
// }

