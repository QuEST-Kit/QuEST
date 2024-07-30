/** @file
 * Inlined bitwise operations used by all deployment modes for fast,
 * low-level processing of basis state indices (qindex). 
 */

#ifndef BITWISE_HPP
#define BITWISE_HPP

#include "types.h"



/*
 * This header forcefully inlines all functions so that it
 * can be included by multiple independently-compiled source
 * files without symbol duplication; especially important for
 * their invocation by CUDA kernels. It also provides a 
 * performance benefit when these functions are called in tight loops.
 */

#if COMPILE_CUDA && (defined(__NVCC__) || defined(__HIPCC__))
    #define INLINE __forceinline__ __device__ __host__
#else
    #define INLINE inline __attribute__((always_inline))
#endif



/* 
 * Performance-critical functions (called in tight loops)
 */


#define QINDEX_ONE 1ULL


INLINE qindex powerOf2(int exponent) {
    
    return QINDEX_ONE << exponent;
}


INLINE bool isPowerOf2(qindex number) {

    return (number > 0) && ((number & (number - QINDEX_ONE)) == 0);
}


INLINE int getBit(qindex number, int bitIndex) {
    
    return (number >> bitIndex) & QINDEX_ONE;
}


INLINE qindex flipBit(qindex number, int bitIndex) {
    
    return number ^ (QINDEX_ONE << bitIndex);
}


INLINE qindex insertBit(qindex number, int bitIndex, int bitValue) {
    
    qindex left = (number >> bitIndex) << (bitIndex + 1);
    qindex middle = bitValue << bitIndex;
    qindex right = number & ((QINDEX_ONE << bitIndex) - 1);
    return left | middle | right;
}


INLINE qindex insertBits(qindex number, int* bitIndices, int numIndices, int bitValue) {
    
    // bitIndices must be strictly increasing
    for (int i=0; i<numIndices; i++)
        number = insertBit(number, bitIndices[i], bitValue);
        
    return number;
}


INLINE qindex setBit(qindex number, int bitIndex, int bitValue) {
    
    qindex mask = bitValue << bitIndex;
    return (number & (~mask)) | mask;
}


INLINE qindex setBits(qindex number, int* bitIndices, int numIndices, qindex bitsValue) {
    
    for (int i=0; i<numIndices; i++) {
        int bit = getBit(bitsValue, i);
        number = setBit(number, bitIndices[i], bit);
    }
    
    return number;
}


INLINE int getBitMaskParity(qindex mask) {
    
    // this compiler extension may not be defined on all platforms
    return __builtin_parity(mask);
}



/*
 * Convenience wrappers around performance-critical functions
 */
 

INLINE qindex insertTwoBits(qindex number, int highInd, int highBit, int lowInd, int lowBit) {
    
    number = insertBit(number, lowInd, lowBit);
    number = insertBit(number, highInd, highBit);
    return number;
}


INLINE qindex insertThreeZeroBits(qindex number, int i3, int i2, int i1) {
    
    number = insertTwoBits(number, i2, 0, i1, 0);
    number = insertBit(number, i3, 0);
    return number;
}


INLINE qindex insertFourZeroBits(qindex number, int i4, int i3, int i2, int i1) {
    
    number = insertTwoBits(number, i2, 0, i1, 0);
    number = insertTwoBits(number, i4, 0, i3, 0);
    return number;
}

INLINE qindex flipTwoBits(qindex number, int i1, int i0) {
    
    number = flipBit(number, i1);
    number = flipBit(number, i0);
    return number;
}



/* 
 * Non-performance critical convenience functions, which should
 * not be used in exponentially-big tight-loops. We inline anyway
 * to avoid symbol duplication issues
 */


INLINE int getNextLeftmostZeroBit(qindex mask, int bitInd) {

    bitInd--;
    while (getBit(mask, bitInd))
        bitInd--;
    
    return bitInd;
}


INLINE bool allBitsAreOne(qindex number, int* bitIndices, int numIndices) {
    
    for (int i=0; i<numIndices; i++)
        if (!getBit(number, bitIndices[i]))
            return false;
            
    return true;
}


INLINE qindex getBitMask(int* bitIndices, int numIndices) {
    
    qindex mask = 0;
    for (int i=0; i<numIndices; i++)
        mask = flipBit(mask, bitIndices[i]);
        
    return mask;
}


INLINE int logBase2(qindex powerOf2) {
    
    int expo = 0;
    while (getBit(powerOf2, 0) != 1) {
        expo++;
        powerOf2 >>= 1;
    }

    return expo;
}



#endif // BITWISE_HPP