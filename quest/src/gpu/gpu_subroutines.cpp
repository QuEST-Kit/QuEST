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
using namespace index_flags;



/*
 * ANY-CTRL ONE-TARG MATRIX TEMPLATES
 */


template <CtrlFlag ctrlFlag>
void gpu_statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
#if COMPILE_CUQUANTUM

    // TODO: call one of the cuQuantum funcs
    // e.g.
    // cuquantum_statevec_oneTargetGate_subA(qureg, targ, matrix);

#elif COMPILE_CUDA

    // each thread modifies 2 amps, and each ctrl qubit halves needed threads 
    qindex numThrds = qureg.numAmpsPerNode / powerOf2(1 + ctrls.size());
    qindex numBlocks = getNumBlocks(numThrds);

    // there are sufficiently few matrix elements (<256 bytes) to pass directly to kernel
    cu_qcomp m00, m01, m10, m11;
    unpackMatrixToCuQcomps(matr, m00,m01,m10,m11);

    // prepare the parameters necessary for efficiently computing amp indices inside kernels
    CtrlTargIndParams params = getParamsInformingIndsWhereCtrlsAreActiveAndTargIsOne(ctrls, ctrlStates, targ);

    kernel_statevector_anyCtrlOneTargCompMatr_subA <ctrlFlag> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThrds, params, targ, m00,m01,m10,m11);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


template <CtrlFlag ctrlFlag>
void gpu_statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {
#if COMPILE_CUQUANTUM

    // TODO: call one of the cuQuantum funcs

#elif COMPILE_CUDA

    // each thread modifies 1 amp, and each ctrl qubit halves needed threads 
    qindex numThrds = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    qindex numBlocks = getNumBlocks(numThrds);

    // there are sufficiently few matrix elements (<256 bytes) to pass directly to kernel
    cu_qcomp d0, d1;
    unpackMatrixToCuQcomps(matr, d0,d1);

    // prepare the parameters necessary for efficiently computing amp indices inside kernels
    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(ctrls, ctrlStates);

    kernel_statevector_anyCtrlOneTargDiagMatr_subA <ctrlFlag> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (
        toCuQcomps(qureg.gpuAmps), numThrds, params, targ, d0,d1);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}



/*
 *  INSTANTIATING TEMPLATES
 */

INSTANTIATE_TEMPLATED_VOID_FUNC_WITH_ALL_CTRL_FLAGS(
    gpu_statevector_anyCtrlOneTargMatrix_subA, 
    (Qureg, vector<int>, vector<int>, int, CompMatr1))

INSTANTIATE_TEMPLATED_VOID_FUNC_WITH_ALL_CTRL_FLAGS(
    gpu_statevector_anyCtrlOneTargMatrix_subA, 
    (Qureg, vector<int>, vector<int>, int, DiagMatr1))
