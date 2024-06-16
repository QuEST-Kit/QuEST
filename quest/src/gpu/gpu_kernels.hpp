/** @file
 * Custom CUDA kernels invoked by gpu.cpp, usually only necessary when
 * there is no equivalent utility in Thrust (or cuQuantum, when it is
 * targeted). This file is only ever included when ENABLE_GPU_ACCELERATION=1 
 * so it can safely invoke CUDA signatures without guards.
 */

#ifndef GPU_KERNELS_HPP
#define GPU_KERNELS_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/gpu/gpu_types.hpp"

#if ! ENABLE_GPU_ACCELERATION
    #error "A file being compiled somehow included gpu_kernels.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif



/*
 * THREAD MANAGEMENT
 */

const int NUM_THREADS_PER_BLOCK = 128;

__forceinline__ __device__ qindex getThreadInd() {
    return blockIdx.x*blockDim.x + threadIdx.x;
}

__host__ qindex getNumBlocks(qindex numIts) {
    return ceil(numIts / (qreal) NUM_THREADS_PER_BLOCK);
}



/*
 * OPERATORS 
 */


__global__ void kernel_statevec_oneTargetGate_subA(cu_qcomp* amps, qindex numIts, int target, cu_qcomp m00, cu_qcomp m01, cu_qcomp m10, cu_qcomp m11) {

    qindex j = getThreadInd();
    if (j >= numIts) 
        return;

    qindex i0 = insertBit(j, target, 0);
    qindex i1 = flipBit(i0, target);

    cu_qcomp amp0 = amps[i0];
    cu_qcomp amp1 = amps[i1];

    amps[i0] = m00*amp0 + m01*amp1;
    amps[i1] = m10*amp0 + m11*amp1;
}

__global__ void kernel_statevec_oneTargetGate_subB(cu_qcomp* amps, cu_qcomp* buffer, qindex numIts, cu_qcomp fac0, cu_qcomp fac1) {

    qindex j = getThreadInd();
    if (j >= numIts) 
        return;

    cu_qcomp amp0 = amps[j];
    cu_qcomp amp1 = buffer[j];

    amps[j] = fac0*amp0 + fac1*amp1;
}



#endif // GPU_KERNELS_HPP