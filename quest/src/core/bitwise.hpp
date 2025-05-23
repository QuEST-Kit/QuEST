/** @file
 * Bitwise operations used by all deployment modes for fast,
 * low-level processing of basis state indices (qindex). 
 * 
 * @author Tyson Jones
 * @author Erich Essmann (improved OS agnosticism)
 * @author James Richings (patched setBit)
 */

#ifndef BITWISE_HPP
#define BITWISE_HPP

#ifdef _MSC_VER
  #include <intrin.h>
#endif

#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"



/* 
 * PERFORMANCE-CRITICAL FUNCTIONS
 *
 * which are called in hot loops loops (like by OpenMP threads and
 * CUDA kernels) so are aggressively inlined.
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


INLINE qindex getBitsLeftOfIndex(qindex number, int bitIndex) {

    return number >> (bitIndex + 1);
}


INLINE qindex getBitsRightOfIndex(qindex number, int bitIndex) {

    qindex mask = (QINDEX_ONE << bitIndex) - 1;
    return number & mask;
}


INLINE int getTwoAdjacentBits(qindex number, qindex lowerBitInd) {

    return (number >> lowerBitInd) & 3;
}


INLINE qindex flipBits(qindex number, qindex mask) {

    return number ^ mask;
}


INLINE qindex flipBit(qindex number, int bitIndex) {
    
    qindex mask = QINDEX_ONE << bitIndex;
    return flipBits(number, mask);
}


INLINE qindex concatenateBits(qindex prefix, qindex suffix, int numBitsInSuffix) {

    return (prefix << numBitsInSuffix) | suffix;
}


INLINE qindex concatenateBits(qindex pref, qindex mid, int numMidBits, qindex suf, int numSufBits) {

    int numRight = numMidBits + numSufBits;
    qindex right = concatenateBits(mid, suf, numSufBits);
    qindex all = concatenateBits(pref, right, numRight);
    return all;
}


INLINE qindex insertBit(qindex number, int bitIndex, int bitValue) {
    
    qindex left  = getBitsLeftOfIndex (number, bitIndex-1); // include bit at bitIndex
    qindex right = getBitsRightOfIndex(number, bitIndex);
    qindex all = concatenateBits(left, bitValue, 1, right, bitIndex);
    return all;
}


INLINE qindex setBit(qindex number, int bitIndex, int bitValue) {
    
    // beware that shifting the raw int would overflow (#623)
    qindex bitInPlace = ((qindex) bitValue) << bitIndex;
    qindex oneInPlace = QINDEX_ONE << bitIndex;
    return (number & ~oneInPlace) | bitInPlace;
}


INLINE int getBitMaskParity(qindex mask) {

    // Try a builtin if on GCC/Clang and it is available
#if defined(__has_builtin)
#if __has_builtin(__builtin_parityll)
    return __builtin_parityll(mask);
#endif
#elif defined(__GNUC__) || defined(__clang__)
    // Older GCC/Clang typically have __builtin_parityll by default
    return __builtin_parityll(mask);
#endif
    // Use MSVC-specific popcount intrinsic if available
#ifdef _MSC_VER
    return __popcnt64(mask) & 1;
#endif

    // Fallback: a portable nibble-based trick for parity
    //
    // Explanation:
    //   - XOR the upper half into the lower half until you’re down to 4 bits.
    //   - Then use a 16-bit constant (0x6996) to map 4-bit values to their parity.
    static_assert(sizeof(qindex) >= 4, "qindex must be at least 32 bits");
    if constexpr (sizeof(qindex) >= 8) {
        // If 64-bit, fold mask[63..32] into mask[31..0]
        mask ^= (mask >> 32);
    }
    // Then fold mask[31..16] into mask[15..0], etc.
    mask ^= (mask >> 16);
    mask ^= (mask >> 8);
    mask ^= (mask >> 4);
    mask &= 0xF;
    return (0x6996 >> mask) & 1;
}



/* 
 * LOOPED PERFORMANCE-CRITICAL FUNCTIONS
 *
 * wherein the runtime loops can damage performance when they
 * are embedded in exponentially-large hot loops. As such, these
 * functions should be called with compile-time loop sizes
 * (e.g. through function template parameters, or constexpr)
 * to trigger automatic loop unrolling.
 */


INLINE qindex insertBits(qindex number, int* bitIndices, int numIndices, int bitValue) {
    
    // bitIndices must be strictly increasing
    for (int i=0; i<numIndices; i++)
        number = insertBit(number, bitIndices[i], bitValue);
        
    return number;
}


INLINE qindex setBits(qindex number, int* bitIndices, int numIndices, qindex bitsValue) {
    
    // bitIndices are arbitrarily ordered, which does not affect number
    for (int i=0; i<numIndices; i++) {
        int bit = getBit(bitsValue, i);
        number = setBit(number, bitIndices[i], bit);
    }
    
    return number;
}


INLINE qindex getValueOfBits(qindex number, int* bitIndices, int numIndices) {

    // bits are arbitrarily ordered, which affects value
    qindex value = 0;

    for (int i=0; i<numIndices; i++)
        value |= getBit(number, bitIndices[i]) << i;

    return value;
}



/*
 * PERFORMANCE-CRITICAL CONVENIENCE FUNCTIONS
 *
 * which merely improve caller's code readability
 */


INLINE qindex insertBitsWithMaskedValues(qindex number, int* bitInds, int numBits, qindex mask) {

    // bitInds must be sorted (increasing), and mask must be zero everywhere except bitInds
    return mask | insertBits(number, bitInds, numBits, 0);
}


INLINE int getTwoBits(qindex number, int highInd, int lowInd) {

    int b1 = getBit(number, lowInd);
    int b2 = getBit(number, highInd);
    int v = concatenateBits(b2, b1, 1);
    return v;
}


INLINE qindex insertTwoBits(qindex number, int highInd, int highBit, int lowInd, int lowBit) {
    
    // assumes highInd > lowInd
    number = insertBit(number, lowInd, lowBit);
    number = insertBit(number, highInd, highBit);
    return number;
}


INLINE qindex insertThreeZeroBits(qindex number, int i3, int i2, int i1) {
    
    // assumes i3 > i2 > i1
    number = insertTwoBits(number, i2, 0, i1, 0);
    number = insertBit(number, i3, 0);
    return number;
}


INLINE qindex insertFourZeroBits(qindex number, int i4, int i3, int i2, int i1) {
    
    // assumes i4 > i3 > i2 > i1
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
 * SLOW FUNCTIONS
 *
 * which should never be called in hot loops, but which are
 * inlined anyway to avoid symbol duplication. Some may
 * actually be fast (i.e. use a fixed number of integer 
 * operations) but are used exclusively outside of hot loops
 * in the source code, so their speed is inconsequential.
 */


INLINE qindex flipBits(qindex number, int* bitIndices, int numIndices) {

    for (int i=0; i<numIndices; i++)
        number = flipBit(number, bitIndices[i]);
    
    return number;
}


INLINE int getIndOfNextLeftmostZeroBit(qindex mask, int bitInd) {

    bitInd--;
    while (getBit(mask, bitInd))
        bitInd--;
    
    return bitInd;
}


INLINE int getIndOfNextRightmostZeroBit(qindex mask, int bitInd) {

    bitInd++;
    while (getBit(mask, bitInd))
        bitInd++;

    return bitInd;
}


INLINE bool allBitsAreOne(qindex number, int* bitIndices, int numIndices) {
    
    for (int i=0; i<numIndices; i++)
        if (!getBit(number, bitIndices[i]))
            return false;
            
    return true;
}


INLINE qindex getBitMask(int* bitIndices, int* bitValues, int numIndices) {

    qindex mask = 0;
    for (int i=0; i<numIndices; i++)
        mask = setBit(mask, bitIndices[i], bitValues[i]); 

    return mask;
}


INLINE qindex getBitMask(int* bitIndices, int numIndices) {
    
    qindex mask = 0;
    for (int i=0; i<numIndices; i++)
        mask = flipBit(mask, bitIndices[i]);
        
    return mask;
}


INLINE qindex removeBits(qindex number, int* bitInds, int numInds) {

    // assumes bitIndices are strictly increasing without duplicates
    int numRemoved = 0;

    // remove each bit in-turn
    for (int i=0; i<numInds; i++) {

        // removal of bits invalidates bitInds
        int shiftedInd = bitInds[i] - (numRemoved++);

        qindex lowerBits = getBitsRightOfIndex(number, shiftedInd);
        qindex upperBits = getBitsLeftOfIndex(number, shiftedInd);
        number = concatenateBits(upperBits, lowerBits, shiftedInd);
    }

    return number;
}


INLINE int logBase2(qindex powerOf2) {
    
    int expo = 0;
    while (getBit(powerOf2, 0) != 1) {
        expo++;
        powerOf2 >>= 1;
    }

    return expo;
}


INLINE qindex getIntegerFromBits(int* bits, int numBits) {

    // first bit is treated as least significant
    qindex value = 0;

    for (int i=0; i<numBits; i++)
        value |= bits[i] << i;

    return value;
}


INLINE void getBitsFromInteger(int* bits, qindex number, int numBits) {

    for (int i=0; i<numBits; i++)
        bits[i] = getBit(number, i);
}



#endif // BITWISE_HPP