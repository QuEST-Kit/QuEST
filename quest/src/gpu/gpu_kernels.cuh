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
#include "quest/src/gpu/gpu_types.cuh"

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


#define GET_THREAD_IND(var, numThreads) \
    qindex var = getThreadInd(); \
    if (var >= numThreads) \
        return;



/*
 * KERNEL PRIMITIVES
 */


__forceinline__ __device__  int cudaGetBitMaskParity(qindex mask) {

    // we cannot use bitwise's getBitMaskParity()'s host-only GCC call
    return __popcll(mask) & 1;
}



/* 
 * COMMUNICATION BUFFER PACKING
 */


template <int NumCtrls>
__global__ void kernel_statevec_packAmpsIntoBuffer(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* qubits, int numQubits, qindex qubitStateMask
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumCtrls, numQubits);

    // i = nth local index where qubits are active
    qindex i = insertBitsWithMaskedValues(n, qubits, numBits, qubitStateMask);

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
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    int numQubitBits = 2 + numCtrlBits;

    // i01 = nth local index where ctrls are active, targ2=0 and targ1=1
    qindex i01 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numQubitBits, ctrlsAndTargsMask);
    qindex i10 = flipTwoBits(i01, targ2, targ1);

    // swap amps
    cu_qcomp amp01 = amps[i01];
    amps[i01] = amps[i10];
    amps[i10] = amp01;
}


template <int NumCtrls> 
__global__ void kernel_statevec_anyCtrlSwap_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* ctrls, int numCtrls, qindex ctrlStateMask
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrls are active
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // caller offsets buffer if necessary
    amps[i] = buffer[n];
}


template <int NumCtrls> 
__global__ void kernel_statevec_anyCtrlSwap_subC(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int* ctrlsAndTarg, int numCtrls, qindex ctrlsAndTargMask
) {
    GET_THREAD_IND(n, numThreads);

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
    int* ctrlsAndTarg, int numCtrls, qindex ctrlStateMask, int targ, 
    cu_qcomp m00, cu_qcomp m01, cu_qcomp m10, cu_qcomp m11
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i0 = nth local index where ctrls are active and targ is 0
    qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTarg, numCtrlBits + 1, ctrlStateMask);
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
    int* ctrls, int numCtrls, qindex ctrlStateMask,
    cu_qcomp fac0, cu_qcomp fac1
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrl bits are active
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // caller offsets buffer by receive-index
    amps[i] = fac0*amps[i] + fac1*buffer[n];
}



/*
 * ANY-TARGET DENSE MATRIX
 */


__forceinline__ __device__ qindex getStrideOfGlobalThreadArr() {
    return gridDim.x * blockDim.x;
}


__forceinline__ __device__ qindex getThreadsNthGlobalArrInd(qindex n, qindex threadInd, qindex stride) {
    return (n * stride) + threadInd;
}


__forceinline__ __device__ qindex getFlattenedMatrInd(qindex row, qindex col, qindex dim) {
    return (row * col) + dim;
}


template <int NumCtrls, int NumTargs>
__global__ void kernel_statevec_anyCtrlAnyTargDenseMatr_sub(
    cu_qcomp* amps, cu_qcomp* cache, qindex numThreads,
    int* ctrlsAndTargs, int numCtrls, qindex ctrlsAndTargsMask, int* targs, int numTargs,
    cu_qcomp* flatMatrElems
) {
    GET_THREAD_IND(n, numThreads);

    // use template params to compile-time unroll loops in insertBits() and setBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numTargs);
    int numQubitBits = numCtrlBits + numTargBits;
    qindex numTargAmps = powerOf2(numTargBits);

    // i0 = nth local index where ctrls are active and targs are all zero
    qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numQubitBits, ctrlsAndTargsMask);
    qindex stride = getStrideOfGlobalThreadArr();

    // collect and cache all to-be-modified amps (loop might be unrolled)
    for (qindex k=0; k<numTargAmps; k++) {

        // i = nth local index where ctrls are active and targs form value k
        qindex i = setBits(i0, targs, numTargBits, k); // loop may be unrolled
        qindex j = getThreadsNthGlobalArrInd(k, n, stride);
        cache[j] = amps[i];
    }

    // modify each amplitude (loop might be unrolled)
    for (qindex k=0; k<numTargAmps; k++) {

        // i = nth local index where ctrls are active and targs form value k
        qindex i = setBits(i0, targs, numTargBits, k); // loop may be unrolled
        amps[i] = {0, 0}; // zero comp literal
    
        // loop may be unrolled
        for (qindex l=0; l<numTargAmps; l++) {

            qindex j = getThreadsNthGlobalArrInd(l, n, stride);
            qindex h = getFlattenedMatrInd(k, l, numTargAmps);
            amps[i] += flatMatrElems[h] * cache[j];
        }
    }
}



/*
 * DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs>
__global__ void kernel_statevec_anyCtrlAnyTargDiagMatr_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex ctrlStateMask, int* targs, int numTargs,
    cu_qcomp* elems
) {
    GET_THREAD_IND(n, numThreads);

    // use template params to compile-time unroll loops in insertBits() and getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numTargs);

    // j = nth local index where ctrls are active (in the specified states)
    qindex j = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // i = global index corresponding to j
    qindex i = concatenateBits(rank, j, logNumAmpsPerNode);

    // t = value of targeted bits, which may be in the prefix substate
    qindex t = getValueOfBits(i, targs, numTargBits);

    amps[i] *= elems[t];
}



/*
 * DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs> 
__global__ void kernel_statevector_anyCtrlPauliTensorOrGadget_subA(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrlsAndTargs, int numCtrls, qindex ctrlsAndTargsStateMask, 
    int* suffixTargsXY, int numSuffixTargsXY,
    qindex suffixMaskXY, qindex allMaskYZ, 
    cu_qcomp powI, cu_qcomp thisAmpFac, cu_qcomp otherAmpFac
) {
    GET_THREAD_IND(n, numThreads);

    // use template params to compile-time unroll loops in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numSuffixTargsXY);

    // compiler will infer these at compile-time if possible
    int numQubitBits = numCtrlBits + numTargBits;
    qindex numTargAmps = powerOf2(numTargBits);

    // each inner iteration modifies 2 amplitudes (may be compile-time sized) 
    qindex numInnerIts = numTargAmps / 2;

    // i0 = nth local index where ctrls are active and targs are all zero
    qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numQubitBits, ctrlsAndTargsStateMask);

    // loop may be unrolled
    for (qindex m=0; m<numInnerIts; m++) {

        // iA = nth local index where targs have value m, iB = (last - nth) such index
        qindex iA = setBits(i0, suffixTargsXY, numSuffixTargsXY, m);
        qindex iB = flipBits(iA, suffixMaskXY);

        // jA = global index corresponding to iA
        qindex jA = concatenateBits(rank, iA, logNumAmpsPerNode);
        qindex jB = concatenateBits(rank, iB, logNumAmpsPerNode);

        // determine whether to multiply amps by +-1 or +-i
        cu_qcomp pmPowA = powI * (1. - 2. * cudaGetBitMaskParity(jA & allMaskYZ));
        cu_qcomp pmPowB = powI * (1. - 2. * cudaGetBitMaskParity(jB & allMaskYZ));

        cu_qcomp ampA = amps[iA];
        cu_qcomp ampB = amps[iB];

        // mix or swap scaled amp pair
        amps[iA] = (thisAmpFac * ampA) + (otherAmpFac * pmPowB * ampB);
        amps[iB] = (thisAmpFac * ampB) + (otherAmpFac * pmPowA * ampA);
    }
}


template <int NumCtrls>
__global__ void kernel_statevector_anyCtrlPauliTensorOrGadget_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex ctrlStateMask,
    qindex suffixMaskXY, qindex bufferMaskXY, qindex allMaskYZ, 
    cu_qcomp powI, cu_qcomp thisAmpFac, cu_qcomp otherAmpFac
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrl bits are in specified states
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // j = buffer index of amp to be mixed with i
    qindex j = flipBits(n, bufferMaskXY);

    // k = global index of amp at buffer index j
    qindex k = concatenateBits(rank, flipBits(i, suffixMaskXY), logNumAmpsPerNode);

    // determine whether to multiply buffer amp by +-1 or +-i
    int negParity = cudaGetBitMaskParity(k & allMaskYZ);
    cu_qcomp pmPowI = powI * (1. - 2. * negParity);

    amps[i] *= thisAmpFac;
    amps[i] += otherAmpFac * pmPowI * buffer[j];
}



#endif // GPU_KERNELS_HPP