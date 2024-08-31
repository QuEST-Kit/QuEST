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


__global__ void kernel_statevec_packPairSummedAmpsIntoBuffer(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int qubit1, int qubit2, int qubit3, int bit2
) {
    GET_THREAD_IND(n, numThreads);

    // i000 = nth local index where all qubits are 0
    qindex i000 = insertThreeZeroBits(n, qubit3, qubit2, qubit1);
    qindex i0b0 = setBit(i000, qubit2, bit2);
    qindex i1b1 = flipTwoBits(i0b0, qubit3, qubit1);

    buffer[n] = amps[i0b0] + amps[i1b1];
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
        amps[i] = {0, 0}; // zero cu_comp literal
    
        // loop may be unrolled
        for (qindex l=0; l<numTargAmps; l++) {

            qindex j = getThreadsNthGlobalArrInd(l, n, stride);
            qindex h = getFlattenedMatrInd(k, l, numTargAmps);
            amps[i] += flatMatrElems[h] * cache[j];
        }
    }
}



/*
 * ANY-TARG DIAGONAL MATRIX
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
 * PAULI/PHASE TENSORS/GADGETS
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


template <int NumCtrls>
__global__ void kernel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex ctrlStateMask, qindex targMask,
    cu_qcomp fac0, cu_qcomp fac1
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrl bits are in specified states
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // j = global index corresponding to i
    qindex j = concatenateBits(rank, i, logNumAmpsPerNode);

    // apply phase to amp depending on parity of targets in global index 
    int p = cudaGetBitMaskParity(j & targMask);

    cu_qcomp facs[] = {fac0, fac1};
    amps[j] *= facs[p];
}



/*
 * DEPHASING
 */


__global__ void kernel_densmatr_oneQubitDephasing_subA(
    cu_qcomp* amps, qindex numThreads, 
    int ketQubit, int braQubit, qreal fac
) {
    GET_THREAD_IND(n, numThreads);

    // TODO:
    // each kernel modifies two amps strided by 2^qureg.numQubits, which is terrible!
    // we can easy template this kernel to modify only 1 thread-local amp, and invoke
    // two kernels at launch. Benchmark this and update

    // i01 = nth local index of |*0*><*1*|
    qindex i01 = insertTwoBits(n, braQubit, 0, ketQubit, 1);
    qindex i10 = insertTwoBits(n, braQubit, 1, ketQubit, 0);

    amps[i01] *= fac;
    amps[i10] *= fac;
}


__global__ void kernel_densmatr_oneQubitDephasing_subB(
    cu_qcomp* amps, qindex numThreads, 
    int ketQubit, int braBit, qreal fac
) {
    GET_THREAD_IND(n, numThreads);

    // TODO:
    // we're just modifying every 2*2^(ketQubit)-th element; turn
    // this into a trivial thrust call (to reduce boilerplate)

    // i = nth local index where bra-qubit differs from ket-qubit
    qindex i = insertBit(n, ketQubit, ! braBit);
    amps[i] *= fac;
}


// there is no bespoke kernel_densmatr_twoQubitDephasing_subA(), since _subB() is directly callable


__global__ void kernel_densmatr_twoQubitDephasing_subB(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode, // numAmps, not numCols
    int ketQubit1, int ketQubit2, int braQubit1, int braQubit2, qreal term
) {
    GET_THREAD_IND(n, numThreads);

    // i = global index of nth local amp
    qindex i = concatenateBits(rank, n, logNumAmpsPerNode);

    int bitA = getBit(i, ketQubit1) ^ getBit(i, braQubit1);
    int bitB = getBit(i, ketQubit2) ^ getBit(i, braQubit2);

    // determine whether or not to modify this amplitude...
    int flag = bitA | bitB;

    // by multiplying by 1 or (1 + term)
    amps[n] *= (term * flag) + 1;
}



/*
 * ONE-QUBIT DEPOLARISING
 */


__global__ void kernel_densmatr_oneQubitDepolarising_subA(
    cu_qcomp* amps, qindex numThreads, 
    int ketQubit, int braQubit, qreal facAA, qreal facBB, qreal facAB
) {
    GET_THREAD_IND(n, numThreads);

    // i00 = nth local index where both qubits are 0
    qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
    qindex i01 = flipBit(i00, ketQubit);
    qindex i10 = flipBit(i00, braQubit);
    qindex i11 = flipBit(i01, braQubit);

    // modify 4 amps, mixing a pair, and scaling the other
    cu_qcomp amp00 = amps[i00];
    amps[i00] = (facAA * amp00) + (facBB * amps[i11]);
    amps[i01] *= facAB;
    amps[i10] *= facAB;
    amps[i11] = (facAA * amps[i11]) + (facBB * amp00);
}


__global__ void kernel_densmatr_oneQubitDepolarising_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int braBit, int ketQubit, qreal facAA, qreal facBB, qreal facAB
) {
    GET_THREAD_IND(n, numThreads);

    // iAA = nth local index where ket qubit agrees with bra qubit
    qindex iAA = insertBit(n, ketQubit, braBit);
    amps[iAA] *= facAA;
    amps[iAA] += facBB * buffer[n];

    // iAB = nth local index where ket qubit disagrees with bra qubit
    qindex iAB = insertBit(n, ketQubit, ! braBit);
    amps[iAB] *= facAB;
}



/*
 * TWO-QUBIT DEPOLARISING
 */


__global__ void kernel_densmatr_twoQubitDepolarising_subA(
    cu_qcomp* amps, qindex numThreads, 
    int ketQb1, int ketQb2, int braQb1, int braQb2, qreal c3
) {
    GET_THREAD_IND(n, numThreads);

    // determine whether to modify amp
    int flag1 = !(getBit(n, ketQb1) ^ getBit(n, braQb1));
    int flag2 = !(getBit(n, ketQb2) ^ getBit(n, braQb2));
    int mod   = !(flag1 & flag2);

    // multiply amp by 1 or (1 + c3)
    amps[n] *= 1 + c3 * mod;
}


__global__ void kernel_densmatr_twoQubitDepolarising_subB(
    cu_qcomp* amps, qindex numThreads, 
    int ketQb1, int ketQb2, int braQb1, int braQb2, qreal c1, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i0000 = nth local index where all bra = ket = 00
    qindex i0000 = insertFourZeroBits(n, braQb2, braQb1, ketQb2, ketQb1);
    qindex i0101 = flipTwoBits(i0000, braQb1, ketQb1);
    qindex i1010 = flipTwoBits(i0000, braQb2, ketQb2);
    qindex i1111 = flipTwoBits(i0101, braQb2, ketQb2);
    
    // mix 1/16 of all amps in groups of 4
    cu_qcomp term = amps[i0000] + amps[i0101] + amps[i1010] + amps[i1111];

    amps[i0000] = c1*amps[i0000] + c2*term;
    amps[i0101] = c1*amps[i0101] + c2*term;
    amps[i1010] = c1*amps[i1010] + c2*term;
    amps[i1111] = c1*amps[i1111] + c2*term;
}


__global__ void kernel_densmatr_twoQubitDepolarising_subC(
    cu_qcomp* amps, qindex numThreads, 
    int ketQb1, int ketQb2, int braQb1, int braBit2, qreal c3
) {
    GET_THREAD_IND(n, numThreads);

    // TODO:
    // this kernel modifies every amplitude, but I think only
    // 25% are actually being changed; fix this by dispatching
    // 25% fewer kernels which go straight to the modified amps

    // decide whether or not to modify nth local
    bool flag1 = getBit(n, ketQb1) == getBit(n, braQb1); 
    bool flag2 = getBit(n, ketQb2) == braBit2;
    bool mod   = !(flag1 & flag2);

    // scale amp by 1 or (1 + c3)
    amps[n] *= 1 + c3 * mod;
}


__global__ void kernel_densmatr_twoQubitDepolarising_subD(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int ketQb1, int ketQb2, int braQb1, int braBit2, qreal c1, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i000 = nth local index where all suffix bits are 0
    qindex i000 = insertThreeZeroBits(n, braQb1, ketQb2, ketQb1);
    qindex i0b0 = setBit(i000, ketQb2, braBit2);
    qindex i1b1 = flipTwoBits(i0b0, braQb1, ketQb1);

    // mix pair of amps using buffer
    cu_qcomp amp0b0 = amps[i0b0];
    cu_qcomp amp1b1 = amps[i1b1];

    amps[i0b0] = c1*amp0b0 + c2*(amp1b1 + buffer[n]);
    amps[i1b1] = c1*amp1b1 + c2*(amp0b0 + buffer[n]);
}


__global__ void kernel_densmatr_twoQubitDepolarising_subE(
    cu_qcomp* amps, qindex numThreads, 
    int ketQb1, int ketQb2, int braBit1, int braBit2, qreal c3
) {
    GET_THREAD_IND(n, numThreads);

    // TODO:
    // this kernel modifies every amplitude, but I think only
    // 25% are actually being changed; fix this by dispatching
    // 25% fewer kernels which go straight to the modified amps

    // choose whether to modify amp
    bool flag1 = getBit(n, ketQb1) == braBit1; 
    bool flag2 = getBit(n, ketQb2) == braBit2;
    bool mod   = !(flag1 & flag2);
    
    // multiply amp by 1 or (1 + c3)
    amps[n] *=  1 + c3 * mod;
}


__global__ void kernel_densmatr_twoQubitDepolarising_subF(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int ketQb1, int ketQb2, int braBit1, int braBit2, qreal c1, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i = nth local index where suffix ket qubits equal prefix bra qubits
    qindex i = insertTwoBits(n, ketQb2, braBit2, ketQb1, braBit1);

    // mix local amp with received buffer amp
    amps[i] = c1*amps[i] + c2*buffer[n];
}


__global__ void kernel_densmatr_twoQubitDepolarising_subG(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int ketQb1, int ketQb2, int braBit1, int braBit2, qreal c1, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i = nth local index where suffix ket qubits equal prefix bra qubits
    qindex i = insertTwoBits(n, ketQb2, braBit2, ketQb1, braBit1);

    // overwrite local amp with buffer amp
    amps[i] = (c2 / c1) * buffer[n];
}



/*
 * PAULI CHANNEL
 */


__global__ void kernel_densmatr_oneQubitPauliChannel_subA(
    cu_qcomp* amps, qindex numThreads, int ketQubit, int braQubit, 
    qreal facAA, qreal facBB, qreal facAB, qreal facBA
) {
    GET_THREAD_IND(n, numThreads);

    // i00 = nth local index where both qubits are 0
    qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
    qindex i01 = flipBit(i00, ketQubit);
    qindex i10 = flipBit(i00, braQubit);
    qindex i11 = flipBit(i01, braQubit);

    // modify 4 amps in 2 separable pairs
    cu_qcomp amp00 = amps[i00];
    cu_qcomp amp01 = amps[i01];
    cu_qcomp amp10 = amps[i10];
    cu_qcomp amp11 = amps[i11];

    amps[i00] = (facAA * amp00) + (facBB * amp11);
    amps[i01] = (facAB * amp01) + (facBA * amp10);
    amps[i10] = (facAB * amp10) + (facBA * amp01);
    amps[i11] = (facAA * amp11) + (facBB * amp00);
}


__global__ void kernel_densmatr_oneQubitPauliChannel_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int ketQubit, int braBit, qreal facAA, qreal facBB, qreal facAB, qreal facBA
) {
    GET_THREAD_IND(n, numThreads);

    // iAA = nth local index where ket qubit is the same as bra, i.e. |.A.><.A.|
    qindex iAA = insertBit(n, ketQubit, braBit);

    // iAB = nth local index where ket qubit is different from bra, i.e. |.A.><.B.|
    qindex iAB = flipBit(iAA, ketQubit);

    // jBB = buffer index of amp to be mixed with iAA's amp, i.e. |.B.><.B.|
    qindex jBB = iAB;
    qindex jBA = iAA;

    // mix each local amp with a received buffer amp, but not each other
    amps[iAA] = (facAA * amps[iAA]) + (facBB * buffer[jBB]);
    amps[iAB] = (facAB * amps[iAB]) + (facBA * buffer[jBA]);
}



/*
 * AMPLITUDE DAMPING
 */


__global__ void kernel_densmatr_oneQubitDamping_subA(
    cu_qcomp* amps, qindex numThreads,
    int ketQubit, int braQubit, qreal prob, qreal c1, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i00 = nth local index where bra and ket qubits are 0
    qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
    qindex i01 = flipBit(i00, ketQubit);
    qindex i10 = flipBit(i00, braQubit);
    qindex i11 = flipBit(i01, braQubit);
    
    // mix both-zero amp with both-one amp (but not vice versa)
    amps[i00] += prob * amps[i11];

    // scale other amps
    amps[i01] *= c1;
    amps[i10] *= c1;
    amps[i11] *= c2;
}


__global__ void kernel_densmatr_oneQubitDamping_subB(
    cu_qcomp* amps, qindex numThreads,
    int qubit, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // TODO:
    // this extremely simple kernel can be definitely
    // be replaced with a Thrust invocation, to reduce
    // boilerplate

    // i = nth local index where qubit=1
    qindex i = insertBit(n, qubit, 1);
    amps[i] *= c2;
}


__global__ void kernel_densmatr_oneQubitDamping_subC(
    cu_qcomp* amps, qindex numThreads,
    int ketQubit, int braBit, qreal c1
) {
    GET_THREAD_IND(n, numThreads);

    // TODO:
    // this extremely simple kernel can be definitely
    // be replaced with a Thrust invocation, to reduce
    // boilerplate

    // i = nth local index where ket differs from bra
    qindex i = insertBit(n, ketQubit, ! braBit);
    amps[i] *= c1;
}


__global__ void kernel_densmatr_oneQubitDamping_subD(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads,
    int qubit, qreal prob
) {
    GET_THREAD_IND(n, numThreads);

    // i = nth local index where ket is 0
    qindex i = insertBit(n, qubit, 0);
    amps[i] += prob * buffer[n];
}



/*
 * PARTIAL TRACE
 */


template <int NumTargs>
__global__ void kernel_densmatr_partialTrace_sub(
    cu_qcomp* ampsIn, cu_qcomp* ampsOut, qindex numThreads,
    int* ketTargs, int* pairTargs, int* allTargs, int numKetTargs
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll below loops
    SET_VAR_AT_COMPILE_TIME(int, numTargPairs, NumTargs, numKetTargs);

    // may be inferred at compile-time
    int numAllTargs = 2*numTargPairs;
    qindex numIts = powerOf2(numTargPairs);

    // TODO:
    // this implementation assumes that the number of amps in outQureg equals or exceeds the 
    // number of CUDA cores, which may not be true when tracing out almost all qubits. We 
    // should change the parallelisation axis in this scenario, or preclude it with validation!

    // k = nth local index of inQureg where all targs and pairs are zero
    qindex k = insertBits(n, allTargs, numAllTargs, 0); // loop may be unrolled

    // each outQureg amp results from summing 2^targs inQureg amps
    cu_qcomp outAmp = {0, 0}; // zero cu_comp literal

    // loop may be unrolled
    for (qindex j=0; j<numIts; j++) {

        // i = nth local index of inQureg where targs=j and pairTargs=j
        qindex i = k;
        i = setBits(i, ketTargs, numTargPairs, j); // loops may be unrolled
        i = setBits(i, pairTargs, numTargPairs, j);

        outAmp += ampsIn[i];
    }

    ampsOut[n] = outAmp;
}



#endif // GPU_KERNELS_HPP