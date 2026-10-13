/** @file
 * Bitwise operations used by all deployment modes for fast,
 * low-level processing of basis state indices (qindex). 
 * 
 * @author Tyson Jones
 * @author Erich Essmann (improved OS agnosticism)
 * @author James Richings (patched setBit)
 * @author PoJen Wang (added BMI2 intrinsics)
 */

#ifndef BITWISE_HPP
#define BITWISE_HPP

#include "quest/include/config.h"
#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"

#if QUEST_COMPILE_BMI2
    #include <immintrin.h>
#endif

#ifdef _MSC_VER
  #include <intrin.h>
#endif



/* 
 * PERFORMANCE-CRITICAL FUNCTIONS
 *
 * which are called in hot loops loops (like by OpenMP threads and
 * CUDA kernels) so are aggressively inlined.
 */


// alternatives to type-unsafe literals 
constexpr qindex QINDEX_ZERO = 0; // used by gpu_thrust.cuh
constexpr qindex QINDEX_ONE  = 1; // used only here


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


INLINE qindex insertBits(qindex number, const int* bitIndices, int numIndices, int bitValue) {
    
    // bitIndices must be strictly increasing
    for (int i=0; i<numIndices; i++)
        number = insertBit(number, bitIndices[i], bitValue);
        
    return number;
}


INLINE qindex setBits(qindex number, const int* bitIndices, int numIndices, qindex bitsValue) {
    
    // bitIndices are arbitrarily ordered, which does not affect number
    for (int i=0; i<numIndices; i++) {
        int bit = getBit(bitsValue, i);
        number = setBit(number, bitIndices[i], bit);
    }
    
    return number;
}


INLINE qindex getValueOfBits(qindex number, const int* bitIndices, int numIndices) {

    // indices are arbitrarily ordered, which affects value; if the indices are
    // known to be sorted, callers should instead use getValueOfPossiblySortedBits()
    // which may (if available) use an optimised intrinsic, eliminating the below
    // loop (though which will anyway be unrolled when numIndices is compile-time)
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


INLINE qindex insertBitsWithMaskedValues(qindex number, const int* bitInds, int numBits, qindex mask) {

    // there exists an overload of insertBitsWithMaskedValues() below which 
    // additionally accepts a (seemingly) superfluous mask encoding bitInds, 
    // and which will use a CPU intrinsic when available
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
 * INTRINSIC-BASED PERFORMANCE-CRITICAL FUNCTIONS
 *
 * which are alternatives to the above functions, and which use 
 * intrinsics for acceleration with specific compilers and on
 * specific CPUs. When the intrinsic is not available, these
 * fallback to the above looped functions.
 */


INLINE qindex getValueOfPossiblySortedBits(qindex number, const bool isSorted, qindex sortedIndsMask, const int* unsortedInds, int numInds) {

    // must not expose BMI2 to GPU backend
#if QUEST_COMPILE_BMI2 && !defined(__NVCC__) && !defined(__HIP__)

    // The BMI2 intrinsic is only usable when inds are sorted (such that the
    // mask is usable). When isSorted is compile-time known, the below branch
    // is eliminated. Otherwise, isSorted is fixed across the caller's hot 
    // loops, and a smart compiler will duplicate the loop and move the branch
    // outside of it. Otherwise, a sensible CPU's branch prediction will
    // eliminate the branch during big hot loops. Otherwise, a very stoopid 
    // compiler and CPU combo will slow small-Qureg simulation via this branch!
    
    return (isSorted)?
        _pext_u64(number, sortedIndsMask):
        getValueOfBits(number, unsortedInds, numInds);

#else

    // suppress unused-var warning
    (void) isSorted;
    (void) sortedIndsMask;

    return getValueOfBits(number, unsortedInds, numInds);

#endif 
}


INLINE qindex insertBitsWithMaskedValues(qindex number, const int* bitInds, int numBits, qindex bitIndsMask, qindex bitValuesMask) {

    // This is an overload of insertBitsWithMaskedValues() above, which accepts the seemingly
    // gratuitous bitIndsMask (which just compactly encodes bitInds), so that a BMI2 intrinsic
    // can be used when available, falling back to the existing looped version. Note bitInds 
    // is always assumed/required to be sorted, regardless of bitIndsMask/instrinsics usage

    // must not expose BMI2 to GPU backend
#if QUEST_COMPILE_BMI2 && !defined(__NVCC__) && !defined(__HIP__)

    // the BMI2 intrinsic only consults bitIndsMask (suppress unused-var warning) 
    (void) bitInds;
    (void) numBits;

    // _pdep_u64 scatters number's bits into set-positions of bitIndsMask, hence "~"
    return bitValuesMask | _pdep_u64(number, ~bitIndsMask);

#else

    // the platform-agnostic version loops through bitInds (and will unroll when numBits is compile-time known)
    (void) bitIndsMask;

    return insertBitsWithMaskedValues(number, bitInds, numBits, bitValuesMask);

#endif
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


INLINE qindex flipBits(qindex number, const int* bitIndices, int numIndices) {

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


INLINE bool allBitsAreOne(qindex number, const int* bitIndices, int numIndices) {
    
    for (int i=0; i<numIndices; i++)
        if (!getBit(number, bitIndices[i]))
            return false;
            
    return true;
}


INLINE qindex getBitMask(const int* bitIndices, const int* bitValues, int numIndices) {

    qindex mask = 0;
    for (int i=0; i<numIndices; i++)
        mask = setBit(mask, bitIndices[i], bitValues[i]); 

    return mask;
}


INLINE qindex getBitMask(const int* bitIndices, int numIndices) {
    
    qindex mask = 0;
    for (int i=0; i<numIndices; i++)
        mask = flipBit(mask, bitIndices[i]);
        
    return mask;
}


INLINE qindex removeBits(qindex number, const int* bitInds, int numInds) {

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


INLINE qindex getIntegerFromBits(const int* bits, int numBits) {

    // first bit is treated as least significant
    qindex value = 0;

    for (int i=0; i<numBits; i++)
        value |= bits[i] << i;

    return value;
}


INLINE void setToBitsOfInteger(int* bits, qindex number, int numBits) {

    for (int i=0; i<numBits; i++)
        bits[i] = getBit(number, i);
}



#endif // BITWISE_HPP
