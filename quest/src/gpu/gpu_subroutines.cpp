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
#include "quest/include/structures.h"

#include "quest/src/core/errors.hpp"

#if ENABLE_GPU_ACCELERATION
    #include "quest/src/gpu/gpu_types.hpp"
#endif



void gpu_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix) {
#if ENABLE_GPU_ACCELERATION
    
    qindex numIts = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numIts);

    // unpack matrix from qcomp into cu_qcomp
    cu_qcomp m00,m01,m10,m11;
    m00 = toCuQcomp(matrix.elems[0][0]);
    m01 = toCuQcomp(matrix.elems[0][1]);
    m10 = toCuQcomp(matrix.elems[1][0]);
    m11 = toCuQcomp(matrix.elems[1][1]);

    kernel_statevec_oneTargetGate_subA<<<numBlocks, NUM_THREADS_PER_BLOCK>>>(
        toCuQcomps(qureg.gpuAmps), numIts, target, m00,m01,m10,m11);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}

void gpu_statevec_oneTargetGate_subB(Qureg qureg, qcomp fac0, qcomp fac1) {
#if ENABLE_GPU_ACCELERATION
    
    qindex numIts = qureg.numAmpsPerNode;
    qindex numBlocks = getNumBlocks(numIts);
    
    kernel_statevec_oneTargetGate_subB<<<numBlocks, NUM_THREADS_PER_BLOCK>>>(
        toCuQcomps(qureg.gpuAmps), 
        toCuQcomps(qureg.gpuCommBuffer), 
        numIts, toCuQcomp(fac0), toCuQcomp(fac1));

#else
    error_gpuSimButGpuNotCompiled();
#endif
}