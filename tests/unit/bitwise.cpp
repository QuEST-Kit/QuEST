/** @file
 * Unit tests for the BMI2 bit gather/scatter helpers added for issue #717
 * (quest/src/core/bitwise.hpp). The optimisation changes only how basis-state
 * indices are computed, so the new mask-accepting helpers must be bit-for-bit
 * identical to the original scalar routines they accelerate. These tests assert
 * exactly that, over exhaustive-small and randomised inputs.
 *
 * When this translation unit is compiled with BMI2 enabled (see
 * tests/unit/CMakeLists.txt) the helpers exercise the _pext_u64 / _pdep_u64
 * intrinsics; otherwise they exercise the scalar fallback. Both must agree with
 * the originals, so the same assertions hold either way.
 *
 * @author (issue #717 contribution)
 *
 * @defgroup unitbitwise Bitwise
 * @ingroup unittests
 */

#include "quest/src/core/bitwise.hpp"

#include <catch2/catch_test_macros.hpp>

#include <vector>
#include <random>
#include <algorithm>

namespace {

    // k distinct indices in [0, maxBit), returned strictly increasing
    std::vector<int> randomIncreasingInds(std::mt19937_64& rng, int k, int maxBit) {
        std::vector<int> pool(maxBit);
        for (int i=0; i<maxBit; i++)
            pool[i] = i;
        std::shuffle(pool.begin(), pool.end(), rng);
        std::vector<int> inds(pool.begin(), pool.begin() + k);
        std::sort(inds.begin(), inds.end());
        return inds;
    }
}

TEST_CASE( "issue #717 helpers compiled path", "[bitwise]" ) {

    // surfaced in the CI log so it is clear which path these tests exercised
#ifdef QUEST_BITWISE_USE_BMI2
    WARN( "bitwise helpers compiled with BMI2 PEXT/PDEP enabled" );
#else
    WARN( "bitwise helpers compiled with the scalar fallback (BMI2 not targeted)" );
#endif
    SUCCEED();
}

TEST_CASE( "getValueOfBitsFromSortedPosMask matches getValueOfBits", "[bitwise]" ) {

    std::mt19937_64 rng(0x717ULL);

    for (int k=0; k<=12; k++) {
        for (int trial=0; trial<200; trial++) {

            std::vector<int> inds = (k==0)
                ? std::vector<int>{}
                : randomIncreasingInds(rng, k, 50);

            qindex posMask = getBitMask(inds.data(), k);

            for (int s=0; s<8; s++) {
                qindex number = (qindex) (rng() & ((1ULL<<50) - 1));   // bits live in [0,50)
                REQUIRE(
                    getValueOfBitsFromSortedPosMask(number, posMask, inds.data(), k) ==
                    getValueOfBits(number, inds.data(), k) );
            }
        }
    }
}

TEST_CASE( "insertBitsWithMaskedValuesAndPosMask matches insertBitsWithMaskedValues", "[bitwise]" ) {

    std::mt19937_64 rng(0x718ULL);

    for (int k=0; k<=12; k++) {
        for (int trial=0; trial<200; trial++) {

            std::vector<int> inds = (k==0)
                ? std::vector<int>{}
                : randomIncreasingInds(rng, k, 50);

            qindex posMask = getBitMask(inds.data(), k);

            // per the original contract, the value mask is zero except at the inserted positions
            qindex valueMask = ((qindex) rng()) & posMask;

            for (int s=0; s<8; s++) {
                qindex number = (qindex) (rng() & ((1ULL<<40) - 1));   // avoid shifting bits past bit 63
                REQUIRE(
                    insertBitsWithMaskedValuesAndPosMask(number, valueMask, posMask, inds.data(), k) ==
                    insertBitsWithMaskedValues(number, inds.data(), k, valueMask) );
            }
        }
    }
}

TEST_CASE( "helpers match at boundary bit positions", "[bitwise]" ) {

    // Deterministic coverage of the awkward positions the randomised tests above never reach:
    // the 32-bit word boundary (31/32) and the high bits 61/62/63 — bit 63 being the sign bit of the
    // signed qindex, where the scalar (arithmetic-shift) and BMI2 (unsigned PEXT/PDEP) paths are most
    // likely to disagree if anything is wrong.
    const std::vector<std::vector<int>> indexSets = {
        {31}, {32}, {63}, {31, 32}, {62, 63}, {0, 63},
        {0, 31, 32, 63}, {30, 31, 32, 33}, {59, 60, 61, 62, 63},
    };
    const std::vector<unsigned long long> numbers = {
        0ULL,
        ~0ULL,                          // all bits set
        1ULL << 63,                     // only the sign bit
        (1ULL << 63) | 1ULL,            // sign bit + bit 0
        0x00000000FFFFFFFFULL,          // low 32
        0xFFFFFFFF00000000ULL,          // high 32
        (1ULL << 31) | (1ULL << 32),    // straddle the word boundary
        0xAAAAAAAAAAAAAAAAULL,          // alternating
        0x5555555555555555ULL,
    };

    for (const auto& inds : indexSets) {
        int k = (int) inds.size();
        qindex posMask = getBitMask(inds.data(), k);

        for (unsigned long long raw : numbers) {

            // gather: any 64-bit input is valid (reads bits, incl. bit 63 of a negative qindex)
            qindex g = (qindex) raw;
            REQUIRE(
                getValueOfBitsFromSortedPosMask(g, posMask, inds.data(), k) ==
                getValueOfBits(g, inds.data(), k) );

            // insert: keep the input within its low (64-k) significant bits so the scalar reference
            // is well-defined (no shift past bit 63); still lets a high input bit land on position 63.
            unsigned long long fitMask = (k == 0) ? ~0ULL : ((1ULL << (64 - k)) - 1);
            qindex n = (qindex) (raw & fitMask);
            for (qindex valueMask : { (qindex) 0, (qindex) (g & posMask) }) {
                REQUIRE(
                    insertBitsWithMaskedValuesAndPosMask(n, valueMask, posMask, inds.data(), k) ==
                    insertBitsWithMaskedValues(n, inds.data(), k, valueMask) );
            }
        }
    }
}

TEST_CASE( "isStrictlyIncreasing detects order", "[bitwise]" ) {

    int sorted[] = {0, 2, 5, 9};
    int equalAdj[] = {0, 2, 2, 9};
    int decreasing[] = {9, 5, 2, 0};

    REQUIRE( isStrictlyIncreasing(sorted, 4) );
    REQUIRE_FALSE( isStrictlyIncreasing(equalAdj, 4) );
    REQUIRE_FALSE( isStrictlyIncreasing(decreasing, 4) );

    // trivially ordered for 0 or 1 elements
    REQUIRE( isStrictlyIncreasing(sorted, 1) );
    REQUIRE( isStrictlyIncreasing(sorted, 0) );
}
