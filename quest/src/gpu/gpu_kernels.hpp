/** @file
 * Custom CUDA kernels invoked by gpu_subroutines.cpp, usually only necessary 
 * when there is no equivalent utility in Thrust (or cuQuantum, when it is
 * targeted). This file is only ever included when COMPILE_CUDA=1 
 * so it can safely invoke CUDA signatures without guards.
 * Some kernels are templated to compile-time optimise their bitwise
 * and indexing logic depending on the number of control qubits.
 */

#ifndef GPU_KERNELS_HPP
#define GPU_KERNELS_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/indexer.hpp"
#include "quest/src/gpu/gpu_types.hpp"

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_kernels.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

using namespace index_flags;



/*
 * THREAD MANAGEMENT
 */


const int NUM_THREADS_PER_BLOCK = 128;


__forceinline__ __device__ qindex getThreadInd() {
    return blockIdx.x*blockDim.x + threadIdx.x;
}


__host__ qindex getNumBlocks(qindex numIts) {

    // TODO:
    // improve this with cudaOccupancyMaxPotentialBlockSize(),
    // making it function specific

    return ceil(numIts / (qreal) NUM_THREADS_PER_BLOCK);
}



/* 
 * COMMUNICATION BUFFER PACKING
 */


template <CtrlFlag ctrlFlag>
qindex kernel_statevec_packAmpsIntoBuffer(cu_qcomp* amps, cu_qcomp* buffer, CtrlIndParams params) {

    qindex n = getThreadInd();
    if (n >= params.numInds) 
        return;

    qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);
    buffer[n] = amps[i];
}



/*
 * ANY-CTRL ONE-TARG MATRIX TEMPLATES
 */


template <CtrlFlag ctrlFlag>
__global__ void kernel_statevec_anyCtrlOneTargDenseMatr_subA(
    cu_qcomp* amps, CtrlTargIndParams params, int targ, 
    cu_qcomp m00, cu_qcomp m01, cu_qcomp m10, cu_qcomp m11
) {
    qindex n = getThreadInd();
    if (n >= params.numInds) 
        return;

    // each thread modifies two amps
    qindex i1 = getNthIndWhereCtrlsAreActiveAndTargIsOne<ctrlFlag>(n, params);
    qindex i0 = flipBit(i1, targ);

    // note they are likely strided and not adjacent
    cu_qcomp amp0 = amps[i0];
    cu_qcomp amp1 = amps[i1];

    amps[i0] = m00*amp0 + m01*amp1;
    amps[i1] = m10*amp0 + m11*amp1;
}


template <CtrlFlag ctrlFlag>
__global__ void kernel_statevec_anyCtrlOneTargDiagMatr_subA(
    cu_qcomp* amps, CtrlIndParams params, 
    int targ, cu_qcomp d0, cu_qcomp d1
) {
    qindex n = getThreadInd();
    if (n >= params.numInds) 
        return;

    // each thread modifies one amp, multiplying by d0 or d1
    qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);
    amps[i] *= d0 + (d1-d0)*getBit(i, targ);
}


template <CtrlFlag ctrlFlag>
__global__ void kernel_statevec_anyCtrlOneTargDenseMatr_subB(
    cu_qcomp* amps, cu_qcomp* buffer, CtrlIndParams params, 
    cu_qcomp fac0, cu_qcomp fac1
) {
    qindex n = getThreadInd();
    if (n >= params.numInds) 
        return;

    // each thread modifies one amp, using one buffer amp
    qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);
    qindex j = n + params.numInds;

    amps[i] = fac0*amps[i] + fac1*buffer[j];
}



#endif // GPU_KERNELS_HPP