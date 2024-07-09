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

#include "modes.h"
#include "types.h"
#include "qureg.h"
#include "structures.h"

#include "../core/errors.hpp"

#if ENABLE_GPU_ACCELERATION
    #include "../gpu/gpu_types.hpp"
    #include "../gpu/gpu_kernels.hpp"
    #include "../gpu/gpu_thrust.hpp"
#endif

#if ENABLE_CUQUANTUM
    #include "../gpu/gpu_cuquantum.hpp"
#endif



void gpu_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix) {
#if ENABLE_CUQUANTUM

    cuquantum_statevec_oneTargetGate_subA(qureg, target, matrix);

#elif ENABLE_GPU_ACCELERATION

    qindex numIts = qureg.numAmpsPerNode / 2;
    qindex numBlocks = getNumBlocks(numIts);

    // sufficiently few matrix elements (<256 bytes) to pass directly to kernel
    cu_qcomp m00, m01, m10, m11;
    unpackMatrixToCuQcomps(matrix, m00,m01,m10,m11);

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
