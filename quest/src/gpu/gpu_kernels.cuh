/** @file
 * Custom CUDA kernels invoked by gpu_subroutines.cpp, usually only necessary 
 * when there is no equivalent utility in Thrust (or cuQuantum, when it is
 * targeted). 
 * 
 * This file is only ever included when COMPILE_CUDA=1 so it can safely invoke 
 * CUDA signatures without guards. Some kernels are templated to compile-time 
 * optimise their bitwise and indexing logic depending on the number of qubits.
 * This file is a header since only ever included by gpu_subroutines.cpp.
 * 
 * When compiling for AMD GPUs, the CUDA symbols invoked herein are
 * mapped to HIP symbols by cuda_to_hip.h 
 * 
 * @author Tyson Jones
 * @author Ania (Anna) Brown (developed QuEST v1 logic)
 */

#ifndef GPU_KERNELS_HPP
#define GPU_KERNELS_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/gpu/gpu_types.cuh"

// kernels/thrust must use cu_qcomp, never qcomp
#define USE_CU_QCOMP
#include "quest/src/core/fastmath.hpp"
#undef USE_CU_QCOMP

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_kernels.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

// cuda keyword 'register' is misinterpreted by HIP
#if defined(__NVCC__)
    #define REGISTER register
#elif defined(__HIP__)
    #define REGISTER
#endif



/*
 * THREAD MANAGEMENT
 */


const int NUM_THREADS_PER_BLOCK = 128;


__forceinline__ __device__ qindex getThreadInd() {
    return blockIdx.x*blockDim.x + threadIdx.x;
}


__host__ qindex getNumBlocks(qindex numThreads) {

    /// @todo
    /// improve this with cudaOccupancyMaxPotentialBlockSize(),
    /// making it function specific

    // CUDA ceil
    return ceil(numThreads / static_cast<qreal>(NUM_THREADS_PER_BLOCK));
}


#define GET_THREAD_IND(var, numThreads) \
    qindex var = getThreadInd(); \
    if (var >= numThreads) \
        return;



/*
 * KERNEL PRIMITIVES
 */


__forceinline__ __device__ int cudaGetBitMaskParity(qindex mask) {

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
 * TWO-TARGET DENSE MATRIX
 */


template <int NumCtrls>
__global__ void kernel_statevec_anyCtrlTwoTargDenseMatr_sub(
    cu_qcomp* amps, qindex numThreads, 
    int* ctrlsAndTarg, int numCtrls, qindex ctrlStateMask, int targ1, int targ2,
    cu_qcomp m00, cu_qcomp m01, cu_qcomp m02, cu_qcomp m03,
    cu_qcomp m10, cu_qcomp m11, cu_qcomp m12, cu_qcomp m13,
    cu_qcomp m20, cu_qcomp m21, cu_qcomp m22, cu_qcomp m23,
    cu_qcomp m30, cu_qcomp m31, cu_qcomp m32, cu_qcomp m33
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i00 = nth local index where ctrls are active and both targs are 0
    qindex i00 = insertBitsWithMaskedValues(n, ctrlsAndTarg, numCtrlBits + 2, ctrlStateMask);
    qindex i01 = flipBit(i00, targ1);
    qindex i10 = flipBit(i00, targ2);
    qindex i11 = flipBit(i01, targ2);

    // note amps00 and amps01 are strided by 2^targ1, and amps00 and amps10 are strided by 2^targ2
    cu_qcomp amp00 = amps[i00];
    cu_qcomp amp01 = amps[i01];
    cu_qcomp amp10 = amps[i10];
    cu_qcomp amp11 = amps[i11];

    // amps[i_n] = sum_j elems[n][j] amp[i_n]
    amps[i00] = m00*amp00 + m01*amp01 + m02*amp10 + m03*amp11;
    amps[i01] = m10*amp00 + m11*amp01 + m12*amp10 + m13*amp11;
    amps[i10] = m20*amp00 + m21*amp01 + m22*amp10 + m23*amp11;
    amps[i11] = m30*amp00 + m31*amp01 + m32*amp10 + m33*amp11;
}



/*
 * ANY-TARGET DENSE MATRIX
 */


__forceinline__ __device__ qindex getThreadsNthGlobalArrInd(qindex n, qindex threadInd, qindex numThreads) {

    // threads store their i-th element contiguously to one another,
    // so that neighbouring threads in a warp access neighbouring 
    // elements of global memory (i.e. coalesce), for caching efficiency
    return (n * numThreads) + threadInd;
}


template <int NumCtrls, int NumTargs, bool ApplyConj>
__global__ void kernel_statevec_anyCtrlFewTargDenseMatr(
    cu_qcomp* amps, qindex numThreads,
    int* ctrlsAndTargs, int numCtrls, qindex ctrlsAndTargsMask, int* targs,
    cu_qcomp* flatMatrElems
) {
    GET_THREAD_IND(n, numThreads);

    // it is gauranteed that NumTargs <= 5, such that the thread-private array is
    // <= 2^5 = 32 amps <= 512 bytes, which aggregated between all threads in the
    // block (assumed ~128) is 64 KiB, and which should be small enough to fit into
    // SM registers without spillage into slow local memory. This is despite it
    // exceeding the maximum per-block shared memory of 48 KiB. Access to this cache
    // must be strictly through compile-time-known indices, otherwise it will auto-
    // spill to local memory). Hence, this _subA() function is not a subroutine 
    // despite some logic being common to non-compile-time _subB(), and hence
    // why the loops below are explicitly compile-time unrolled
    REGISTER cu_qcomp privateCache[1 << NumTargs];

    // we know NumTargs <= 5, though NumCtrls is permitted anything (including -1)
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    constexpr qindex numTargAmps = (1 << NumTargs); // explicit, in lieu of powerOf2

    // i0 = nth local index where ctrls are active and targs are all zero
    qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numCtrlBits + NumTargs, ctrlsAndTargsMask); // loop may be unrolled

    // populate cache (force unroll to ensure compile-time cache indices)
    #pragma unroll  
    for (qindex k=0; k<numTargAmps; k++) {

        // i = nth local index where ctrls are active and targs form value k
        qindex i = setBits(i0, targs, NumTargs, k); // loop will be unrolled

        // write to thread-private cache at compile-time known index
        privateCache[k] = amps[i];
    }

    // modify each amplitude (let compiler decide whether to unroll, to avoid 2^10 = 1024 expansion)
    for (qindex k=0; k<numTargAmps; k++) {

        // i = nth local index where ctrls are active and targs form value k
        qindex i = setBits(i0, targs, NumTargs, k); // loop will be unrolled
        amps[i] = getCuQcomp(0, 0);
    
        // force unroll to ensure compile-time cache indices
        #pragma unroll
        for (qindex l=0; l<numTargAmps; l++) {

            // h = flat index of matrix's (k,l)-th element
            qindex h = fast_getMatrixFlatIndex(k, l, numTargAmps);

            // optionally conjugate matrix elem
            cu_qcomp elem = flatMatrElems[h];
            if constexpr (ApplyConj)
                elem.y *= -1;

            // thread-private cache is accessed with compile-time known index
            amps[i] = amps[i] + (elem * privateCache[l]);
        }
    }
}


template <int NumCtrls, bool ApplyConj>
__global__ void kernel_statevec_anyCtrlManyTargDenseMatr(
    cu_qcomp* globalCache,
    cu_qcomp* amps, qindex numThreads, qindex numBatchesPerThread,
    int* ctrlsAndTargs, int numCtrls, qindex ctrlsAndTargsMask, 
    int* targs, int numTargBits, qindex numTargAmps,
    cu_qcomp* flatMatrElems
) {
    GET_THREAD_IND(t, numThreads);

    // NumCtrls might be compile-time known, but numTargBits>5 is always unknown/runtime
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // unlike all other kernels, each thread modifies multiple batches of amplitudes
    for (qindex b=0; b<numBatchesPerThread; b++) {

        // n = this thread's b-th batch index
        qindex n = t + b * numThreads;

        // i0 = nth local index where ctrls are active and targs are all zero
        qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numCtrlBits + numTargBits, ctrlsAndTargsMask);

        // collect and cache all to-be-modified amps (loop might be unrolled)        
        for (qindex k=0; k<numTargAmps; k++) {

            // i = nth local index where ctrls are active and targs form value k
            qindex i = setBits(i0, targs, numTargBits, k); // loop may be unrolled

            // j = index of k-th element of thread's private cache partition
            qindex j = getThreadsNthGlobalArrInd(k, t, numThreads);
            globalCache[j] = amps[i];
        }

        // modify each amplitude (loop might be unrolled)
        for (qindex k=0; k<numTargAmps; k++) {

            // i = nth local index where ctrls are active and targs form value k
            qindex i = setBits(i0, targs, numTargBits, k); // loop may be unrolled
            amps[i] = getCuQcomp(0, 0);
        
            for (qindex l=0; l<numTargAmps; l++) {
                qindex j = getThreadsNthGlobalArrInd(l, t, numThreads);
                qindex h = fast_getMatrixFlatIndex(k, l, numTargAmps);

                // optionally conjugate matrix elem
                cu_qcomp elem = flatMatrElems[h];
                if constexpr (ApplyConj)
                    elem.y *= -1;

                amps[i] = amps[i] + (elem * globalCache[j]);

                /// @todo
                /// qureg.cpuAmps[i] is being serially updated by only this thread,
                /// so is a candidate for Kahan summation for improved numerical
                /// stability. Explore whether this is time-free and worthwhile!
            }
        }
    }
}



/*
 * ONE-TARG DIAGONAL MATRIX
 */


template <int NumCtrls>
__global__ void kernel_statevec_anyCtrlOneTargDiagMatr_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex ctrlStateMask, int targ, 
    cu_qcomp m1, cu_qcomp m2
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// we have implemented a custom kernel, rather than a thrust
    /// functor, for efficient treatment of control qubits (even
    /// when not exploiting the compile-time parameter NumCtrls).
    /// Our kernel enumerates only amps which satisfy the control
    /// condition, whereas a natural Thrust functor would involve
    /// enumerating all amplitudes and skipping some via a condition.
    /// This might still be beneficial in memory-bandwidth-bound
    /// regimes, but is expected inferior for many control qubits.
    /// We should verify this!

    // use template params to compile-time unroll loops in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // j = nth local index where ctrls are active (in the specified states)
    qindex j = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // i = global index corresponding to j
    qindex i = concatenateBits(rank, j, logNumAmpsPerNode);

    int b = getBit(i, targ);
    amps[j] = amps[j] * (m1 + b * (m2 - m1));
}



/*
 * TWO-TARG DIAGONAL MATRIX
 */


template <int NumCtrls>
__global__ void kernel_statevec_anyCtrlTwoTargDiagMatr_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex ctrlStateMask, int targ1, int targ2,
    cu_qcomp m1, cu_qcomp m2, cu_qcomp m3, cu_qcomp m4
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// we have implemented a custom kernel, rather than a thrust
    /// functor, for efficient treatment of control qubits (even
    /// when not exploiting the compile-time parameter NumCtrls).
    /// Our kernel enumerates only amps which satisfy the control
    /// condition, whereas a natural Thrust functor would involve
    /// enumerating all amplitudes and skipping some via a condition.
    /// This might still be beneficial in memory-bandwidth-bound
    /// regimes, but is expected inferior for many control qubits.
    /// We should verify this!

    // use template params to compile-time unroll loops in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // j = nth local index where ctrls are active (in the specified states)
    qindex j = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // i = global index corresponding to j
    qindex i = concatenateBits(rank, j, logNumAmpsPerNode);

    // k = local elem index
    int k = getTwoBits(i, targ2, targ1);
    cu_qcomp elems[] = {m1, m2, m3, m4};
    amps[j] = amps[j] * elems[k];
}



/*
 * ANY-TARG DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj, bool HasPower>
__global__ void kernel_statevec_anyCtrlAnyTargDiagMatr_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    int* ctrls, int numCtrls, qindex ctrlStateMask, int* targs, int numTargs,
    cu_qcomp* elems, cu_qcomp exponent
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// we have implemented a custom kernel, rather than a thrust
    /// functor, for efficient treatment of control qubits (even
    /// when not exploiting the compile-time parameter NumCtrls).
    /// Our kernel enumerates only amps which satisfy the control
    /// condition, whereas a natural Thrust functor would involve
    /// enumerating all amplitudes and skipping some via a condition.
    /// This might still be beneficial in memory-bandwidth-bound
    /// regimes, but is expected inferior for many control qubits.
    /// We should verify this!

    // use template params to compile-time unroll loops in insertBits() and getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numTargs);

    // j = nth local index where ctrls are active (in the specified states)
    qindex j = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // i = global index corresponding to j
    qindex i = concatenateBits(rank, j, logNumAmpsPerNode);

    // t = value of targeted bits, which may be in the prefix substate
    qindex t = getValueOfBits(i, targs, numTargBits);

    cu_qcomp elem = elems[t];

    if constexpr (HasPower)
        elem = getCompPower(elem, exponent);

    if constexpr (ApplyConj)
        elem.y *= -1;

    amps[j] = amps[j] * elem;
}



/*
 * ALL-TARG DIAGONAL MATRIX
 */


template <bool HasPower, bool MultiplyOnly>
__global__ void kernel_densmatr_allTargDiagMatr_sub(
    cu_qcomp* amps, qindex numThreads, int rank, qindex logNumAmpsPerNode,
    cu_qcomp* elems, qindex numElems, cu_qcomp exponent
) {
    GET_THREAD_IND(n, numThreads);

    // i = global row of nth local index
    qindex i = n % numElems;
    cu_qcomp fac = elems[i];

    if constexpr (HasPower)
        fac = getCompPower(fac, exponent);

    if constexpr (!MultiplyOnly) {

        // m = global index corresponding to n
        qindex m = concatenateBits(rank, n, logNumAmpsPerNode);

        // j = global column corresponding to n
        qindex j = m / numElems;
        cu_qcomp term = elems[j];

        if constexpr(HasPower)
            term = getCompPower(term, exponent);

        // conj after pow
        term.y *= -1;
        fac = fac * term;
    }

    amps[n] = amps[n] * fac;
}



/*
 * PAULI TENSORS/GADGET
 */


template <int NumCtrls, int NumTargs> 
__global__ void kernel_statevector_anyCtrlPauliTensorOrGadget_subA(
    cu_qcomp* amps, qindex numThreads,
    int* ctrlsAndTargs, int numCtrls, qindex ctrlsAndTargsStateMask, 
    int* targsXY, int numXY, qindex maskXY, qindex maskYZ, 
    cu_qcomp powI, cu_qcomp ampFac, cu_qcomp pairAmpFac
) {
    GET_THREAD_IND(t, numThreads);

    // use template params to compile-time unroll loops in insertBits() and setBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numXY);

    // n = local index of amp sub-batch with common i0, v = value of target bits
    qindex numInnerIts = powerOf2(numTargBits) / 2;
    qindex n = t / numInnerIts;
    qindex v = t % numInnerIts;

    // i0 = nth local index where ctrls are active and targs are all zero (loop therein may be unrolled)
    qindex i0 = insertBitsWithMaskedValues(n, ctrlsAndTargs, numCtrlBits + numTargBits, ctrlsAndTargsStateMask);

    // iA = nth local index where targs have value v, iB = (last - nth) such index
    qindex iA = setBits(i0, targsXY, numTargBits, v); // may be unrolled
    qindex iB = flipBits(iA, maskXY);

    // determine whether to multiply amps by +-1 or +-i
    int parA = cudaGetBitMaskParity(iA & maskYZ);
    int parB = cudaGetBitMaskParity(iB & maskYZ);
    cu_qcomp coeffA = powI * fast_getPlusOrMinusOne(parA);
    cu_qcomp coeffB = powI * fast_getPlusOrMinusOne(parB);

    cu_qcomp ampA = amps[iA];
    cu_qcomp ampB = amps[iB];

    // mix or swap scaled amp pair
    amps[iA] = (ampFac * ampA) + (pairAmpFac * coeffB * ampB);
    amps[iB] = (ampFac * ampB) + (pairAmpFac * coeffA * ampA);
}


template <int NumCtrls>
__global__ void kernel_statevector_anyCtrlPauliTensorOrGadget_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads,
    int* ctrls, int numCtrls, qindex ctrlStateMask,
    qindex maskXY, qindex maskYZ, qindex bufferMaskXY,
    cu_qcomp powI, cu_qcomp thisAmpFac, cu_qcomp otherAmpFac
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrl bits are in specified states
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // j = buffer index of amp to be mixed with i
    qindex j = flipBits(n, bufferMaskXY);

    // k = local index of j-th buffer amplitude in its original node
    qindex k = flipBits(i, maskXY);

    // determine whether to multiply buffer amp by +-1 or +-i
    int par = cudaGetBitMaskParity(k & maskYZ);
    cu_qcomp coeff = powI * fast_getPlusOrMinusOne(par);

    amps[i] = (thisAmpFac * amps[i]) + (otherAmpFac * coeff * buffer[j]);
}



/*
 * PHASE TENSORS/GADGET
 */


template <int NumCtrls>
__global__ void kernel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    cu_qcomp* amps, qindex numThreads,
    int* ctrls, int numCtrls, qindex ctrlStateMask, qindex targMask,
    cu_qcomp fac0, cu_qcomp fac1
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, numCtrls);

    // i = nth local index where ctrl bits are in specified states
    qindex i = insertBitsWithMaskedValues(n, ctrls, numCtrlBits, ctrlStateMask);

    // apply phase to amp depending on parity of targets in global index 
    int p = cudaGetBitMaskParity(i & targMask);

    cu_qcomp facs[] = {fac0, fac1};
    amps[i] = amps[i] * facs[p];
}



/*
 * QUREG COMBINATION 
 */


// kernel_densmatr_mixQureg_subA() is avoided; we instead use
// Thrust for this common circumstances (mixing density matrices),
// which should be significantly more optimisex


__global__ void kernel_densmatr_mixQureg_subB(
    qreal outProb, cu_qcomp* outAmps, qreal inProb, cu_qcomp* inAmps,
    qindex numThreads, qindex numInAmps
) {
    GET_THREAD_IND(n, numThreads);

    // (i,j) = row & column of outAmps corresponding to n
    qindex i = n % numInAmps;
    qindex j = n / numInAmps;

    cu_qcomp iAmp = inAmps[i];
    cu_qcomp jAmp = inAmps[j]; jAmp.y *= -1; // conj
    
    outAmps[n] = (outProb * outAmps[n]) + (inProb * iAmp * jAmp);
}


__global__ void kernel_densmatr_mixQureg_subC(
    qreal outProb, cu_qcomp* outAmps, qreal inProb, cu_qcomp* inAmps,
    qindex numThreads, int rank, qindex numInAmps, qindex logNumOutAmpsPerNode
) {
    GET_THREAD_IND(n, numThreads);

    // m = global index of local 'out' index n
    qindex m = concatenateBits(rank, n, logNumOutAmpsPerNode);

    // (i,j) = row & column of outAmps corresponding to n
    qindex i = m % numInAmps;
    qindex j = m / numInAmps;

    cu_qcomp iAmp = inAmps[i];
    cu_qcomp jAmp = inAmps[j]; jAmp.y *= -1; // conj
    
    outAmps[n] = (outProb * outAmps[n]) + (inProb * iAmp * jAmp);
}



/*
 * DEPHASING
 */


__global__ void kernel_densmatr_oneQubitDephasing_subA(
    cu_qcomp* amps, qindex numThreads, 
    int ketQubit, int braQubit, qreal fac
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// each kernel modifies two amps strided by 2^qureg.numQubits, which is terrible!
    /// we can easy template this kernel to modify only 1 thread-local amp, and invoke
    /// two kernels at launch. Benchmark this and update

    // i01 = nth local index of |*0*><*1*|
    qindex i01 = insertTwoBits(n, braQubit, 0, ketQubit, 1);
    qindex i10 = insertTwoBits(n, braQubit, 1, ketQubit, 0);

    amps[i01] = amps[i01] * fac;
    amps[i10] = amps[i10] * fac;
}


__global__ void kernel_densmatr_oneQubitDephasing_subB(
    cu_qcomp* amps, qindex numThreads, 
    int ketQubit, int braBit, qreal fac
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// this extremely simple kernel can be definitely
    /// be replaced with a Thrust invocation, to reduce
    /// boilerplate

    // i = nth local index where bra-qubit differs from ket-qubit
    qindex i = insertBit(n, ketQubit, ! braBit);
    amps[i] = amps[i] * fac;
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
    amps[n] = amps[n] * ((term * flag) + 1);
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
    amps[i01] = amps[i01] * facAB;
    amps[i10] = amps[i10] * facAB;
    amps[i11] = (facAA * amps[i11]) + (facBB * amp00);
}


__global__ void kernel_densmatr_oneQubitDepolarising_subB(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int ketQubit, int braBit, qreal facAA, qreal facBB, qreal facAB
) {
    GET_THREAD_IND(n, numThreads);

    // iAA = nth local index where ket qubit agrees with bra qubit
    qindex iAA = insertBit(n, ketQubit, braBit);
    amps[iAA] = (facAA * amps[iAA]) + (facBB * buffer[n]);

    // iAB = nth local index where ket qubit disagrees with bra qubit
    qindex iAB = insertBit(n, ketQubit, ! braBit);
    amps[iAB] = facAB * amps[iAB];
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
    amps[n] = amps[n] * (1 + c3 * mod);
}


__global__ void kernel_densmatr_twoQubitDepolarising_subB(
    cu_qcomp* amps, qindex numThreads, 
    int ketQb1, int ketQb2, int braQb1, int braQb2, qreal c1alt, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i0000 = nth local index where all bra = ket = 00
    qindex i0000 = insertFourZeroBits(n, braQb2, braQb1, ketQb2, ketQb1);
    qindex i0101 = flipTwoBits(i0000, braQb1, ketQb1);
    qindex i1010 = flipTwoBits(i0000, braQb2, ketQb2);
    qindex i1111 = flipTwoBits(i0101, braQb2, ketQb2);
    
    // mix 1/16 of all amps in groups of 4
    cu_qcomp term = amps[i0000] + amps[i0101] + amps[i1010] + amps[i1111];

    amps[i0000] = c1alt*amps[i0000] + c2*term;
    amps[i0101] = c1alt*amps[i0101] + c2*term;
    amps[i1010] = c1alt*amps[i1010] + c2*term;
    amps[i1111] = c1alt*amps[i1111] + c2*term;
}


__global__ void kernel_densmatr_twoQubitDepolarising_subC(
    cu_qcomp* amps, qindex numThreads, 
    int ketQb1, int ketQb2, int braQb1, int braBit2, qreal c3
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// this kernel modifies every amplitude, but I think only
    /// 25% are actually being changed; fix this by dispatching
    /// 25% fewer kernels which go straight to the modified amps

    // decide whether or not to modify nth local
    bool flag1 = getBit(n, ketQb1) == getBit(n, braQb1); 
    bool flag2 = getBit(n, ketQb2) == braBit2;
    bool mod   = !(flag1 & flag2);

    // scale amp by 1 or (1 + c3)
    amps[n] = amps[n] * (1 + c3 * mod);
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
    int ketQb1, int ketQb2, int braBit1, int braBit2, qreal fac0, qreal fac1
) {
    GET_THREAD_IND(n, numThreads);

    // choose factor by which to scale amp
    bool same1 = getBit(n, ketQb1) == braBit1; 
    bool same2 = getBit(n, ketQb2) == braBit2;
    bool flag = (same1 & same2);

    // scale amp by c1 or (1+c3)
    amps[n] = amps[n] * (fac1 * flag + fac0);
}


__global__ void kernel_densmatr_twoQubitDepolarising_subF(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads, 
    int ketQb1, int ketQb2, int braBit1, int braBit2, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    // i = nth local index where suffix ket qubits equal prefix bra qubits
    qindex i = insertTwoBits(n, ketQb2, braBit2, ketQb1, braBit1);

    // mix local amp with received buffer amp
    amps[i] = amps[i] + (c2 * buffer[n]);
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
    amps[i00] = amps[i00] + (prob * amps[i11]);

    // scale other amps
    amps[i01] = amps[i01] * c1;
    amps[i10] = amps[i10] * c1;
    amps[i11] = amps[i11] * c2;
}


__global__ void kernel_densmatr_oneQubitDamping_subB(
    cu_qcomp* amps, qindex numThreads,
    int qubit, qreal c2
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// this extremely simple kernel can be definitely
    /// be replaced with a Thrust invocation, to reduce
    /// boilerplate

    // i = nth local index where qubit=1
    qindex i = insertBit(n, qubit, 1);
    amps[i] = amps[i] * c2;
}


__global__ void kernel_densmatr_oneQubitDamping_subC(
    cu_qcomp* amps, qindex numThreads,
    int ketQubit, int braBit, qreal c1
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// this extremely simple kernel can be definitely
    /// be replaced with a Thrust invocation, to reduce
    /// boilerplate

    // i = nth local index where ket differs from bra
    qindex i = insertBit(n, ketQubit, ! braBit);
    amps[i] = amps[i] * c1;
}


__global__ void kernel_densmatr_oneQubitDamping_subD(
    cu_qcomp* amps, cu_qcomp* buffer, qindex numThreads,
    int qubit, qreal prob
) {
    GET_THREAD_IND(n, numThreads);

    // i = nth local index where ket is 0
    qindex i = insertBit(n, qubit, 0);
    amps[i] = amps[i] + (prob * buffer[n]);
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

    /// @todo
    /// this implementation assumes that the number of amps in outQureg equals or exceeds the 
    /// number of CUDA cores, which may not be true when tracing out almost all qubits. We 
    /// should change the parallelisation axis in this scenario, or preclude it with validation!

    // k = nth local index of inQureg where all targs and pairs are zero
    qindex k = insertBits(n, allTargs, numAllTargs, 0); // loop may be unrolled

    // each outQureg amp results from summing 2^targs inQureg amps
    cu_qcomp outAmp = getCuQcomp(0, 0);

    // loop may be unrolled
    for (qindex j=0; j<numIts; j++) {

        // i = nth local index of inQureg where targs=j and pairTargs=j
        qindex i = k;
        i = setBits(i, ketTargs,  numTargPairs, j); // loops may be unrolled
        i = setBits(i, pairTargs, numTargPairs, j);

        outAmp = outAmp + ampsIn[i];
    }

    ampsOut[n] = outAmp;
}



/*
 * PROBABILITIES
 */


template <int NumQubits>
__global__ void kernel_statevec_calcProbsOfAllMultiQubitOutcomes_sub(
    qreal* outProbs, cu_qcomp* amps, qindex numThreads, 
    int rank, qindex logNumAmpsPerNode,
    int* qubits, int numQubits
) {
    GET_THREAD_IND(n, numThreads);

    /// @todo
    /// it might be possible to replace this custom kernel 
    /// with an invocation of Thrust's reduce_by_key(),
    /// where the key is j as computed below. Look into
    /// whether this is worthwhile and faster!

    // use template param to compile-time unroll below loops
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, numQubits);

    qreal prob = getCompNorm(amps[n]);

    // i = global index corresponding to n
    qindex i = concatenateBits(rank, n, logNumAmpsPerNode);

    // j = outcome index corresponding to prob
    qindex j = getValueOfBits(i, qubits, numBits); // loop therein may be unrolled

    atomicAdd(&outProbs[j], prob);
}


template <int NumQubits>
__global__ void kernel_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(
    qreal* outProbs, cu_qcomp* amps, qindex numThreads, 
    qindex firstDiagInd, qindex numAmpsPerCol,
    int rank, qindex logNumAmpsPerNode,
    int* qubits, int numQubits
) {
    GET_THREAD_IND(n, numThreads);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, numQubits);

    // i = index of nth local diagonal elem
    qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);
    qreal prob = getCompReal(amps[i]);

    // j = global index of i
    qindex j = concatenateBits(rank, i, logNumAmpsPerNode);

    // k = outcome index corresponding to 
    qindex k = getValueOfBits(j, qubits, numBits); // loop therein may be unrolled

    atomicAdd(&outProbs[k], prob);
}


#endif // GPU_KERNELS_HPP
