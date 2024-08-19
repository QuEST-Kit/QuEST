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

#include "quest/src/gpu/gpu_types.hpp"

#if ! COMPILE_CUDA
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

    // TODO:
    // improve this with cudaOccupancyMaxPotentialBlockSize(),
    // making it function specific

    return ceil(numIts / (qreal) NUM_THREADS_PER_BLOCK);
}



/* 
 * COMMUNICATION BUFFER PACKING
 */


template <int NumCtrls>
__global__ void kernel_statevec_packAmpsIntoBuffer(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* qubits, int numQubits, qindex mask
) {
    qindex n = getThreadInd();
    if (n >= numThreads) 
        return;

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumCtrls, numQubits);

    // i = nth local index where qubits are active
    qindex i = insertBitsWithMaskedValues(n, qubits, numBits, mask);

    // caller offsets buffer by sub-buffer send-index
    buffer[n] = amps[i];
}



/*
 * SWAPS
 */

template <int NumCtrls> 
__global__ void kernel_statevec_anyCtrlSwap_subA(
    cu_qcomp* amps, qindex numThreads, 
    int* ctrlsAndTargs, int numCtrls, qindex ctrlsAndTargsMask, int targ1, int targ2
) {
    qindex n = getThreadInd();
    if (n >= numThreads) 
        return;

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    int numQubitBits = 2 + numCtrlBits;

    // i01 = nth local index where ctrls are active, targ2=0 and targ1=1
    qindex i01 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numQubitBits, ctrlsAndTargsMask);
    qindex i10 = flipTwoBits(i01, targ2, targ1);

    // swap amps
    cu_qcmop amp01 = amps[i01];
    amps[i01] = amps[i10];
    amps[i10] = amp01;
}


template <int NumCtrls> 
__global__ void kernel_statevec_anyCtrlSwap_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* ctrls, int numCtrls, qindex mask
) {
    qindex n = getThreadInd();
    if (n >= numThreads)
        return;

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrls are active
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, mask);

    // caller offsets buffer if necessary
    amps[i] = buffer[n];
}


template <int NumCtrls> 
__global__ void kernel_statevec_anyCtrlSwap_subC(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* ctrlsAndTarg, int numCtrls, qindex ctrlsAndTargMask
) {
    qindex n = getThreadInd();
    if (n >= numThreads)
        return;

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    int numQubitBits = numCtrlBits + 1;

    // i = nth local index where ctrls and targ are in specified states
    qindex i = insertBitsWithMaskedValues(n, ctrlsAndTarg, numQubitBits, ctrlsAndTargMask);

    // caller offsets buffer if necessary
    amps[i] = buffer[n];
}



/*
 * ONE-TARGET DENSE MATRIX
 */


template <int NumCtrls>
__global__ void kernel_statevec_anyCtrlOneTargDenseMatr_subA(
    cu_qcomp* amps, qindex numThreads, 
    int* ctrlsAndTarg, int numCtrls, qindex mask, int targ, 
    cu_qcomp m00, cu_qcomp m01, cu_qcomp m10, cu_qcomp m11
) {
    qindex n = getThreadInd();
    if (n >= numThreads) 
        return;

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i0 = nth local index where ctrls are active and targ is 0
    qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTarg, numCtrlBits + 1, mask);
    qindex i1 = flipBit(i0, targ);

    // note amps are strided by 2^targ
    cu_qcomp amp0 = amps[i0];
    cu_qcomp amp1 = amps[i1];

    amps[i0] = m00*amp0 + m01*amp1;
    amps[i1] = m10*amp0 + m11*amp1;
}


template <int NumCtrls>
__global__ void kernel_statevec_anyCtrlOneTargDenseMatr_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* ctrls, int numCtrls, qindex mask,
    cu_qcomp fac0, cu_qcomp fac1
) {
    qindex n = getThreadInd();
    if (n >= numThreads) 
        return;

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrl bits are active
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, mask);

    // caller offsets buffer by receive-index
    amps[i] = fac0*amps[i] + fac1*buffer[n];
}



/*
 * ANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs>
__global__ void kernel_statevec_anyCtrlAnyTargDenseMatr_subA(
    cu_qcomp* amps, qindex numThreads,
    int* qubits, int numQubits, qindex mask, int* targs, int numTargs,
    cu_qcomp* elems
) {
    // TODO
}



/*
 * DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs>
__global__ void kernel_statevec_anyCtrlAnyTargDiagMatr_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex mask, int* targs, int numTargs,
    cu_qcomp* elems
) {
    qindex n = getThreadInd();
    if (n >= numThreads) 
        return;

    // use template params to compile-time unroll loops in insertBits() and getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numTargs);

    // j = nth local index where ctrls are active (in the specified states)
    qindex j = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, mask);

    // i = global index corresponding to j
    qindex i = concatenateBits(rank, j, logNumAmpsPerNode);

    // t = value of targeted bits, which may be in the prefix substate
    qindex t = getValueOfBits(i, targs, numTargBits);

    amps[i] *= elems[t];
}



#endif // GPU_KERNELS_HPP