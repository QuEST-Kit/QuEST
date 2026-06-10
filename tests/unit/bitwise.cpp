/** @file
 * Unit tests of internal bitwise subroutines.
 *
 * @defgroup unitbitwise Bitwise
 * @ingroup unittests
 */

#include <catch2/catch_test_macros.hpp>

#include "quest/src/core/bitwise.hpp"
#include "tests/utils/macros.hpp"

#include <vector>

using std::vector;


/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[bitwise]"


static qindex getRefInsertedBits(qindex number, const vector<int>& bitIndices, int bitValue) {

    qindex out = 0;
    int srcInd = 0;
    int nextIns = 0;

    for (int dstInd=0; dstInd<63; dstInd++) {
        if (nextIns < static_cast<int>(bitIndices.size()) && dstInd == bitIndices[nextIns]) {
            if (bitValue)
                out |= QINDEX_ONE << dstInd;
            nextIns++;
        } else {
            if ((number >> srcInd) & QINDEX_ONE)
                out |= QINDEX_ONE << dstInd;
            srcInd++;
        }
    }

    return out;
}


static qindex getRefValueOfBits(qindex number, const vector<int>& bitIndices) {

    qindex out = 0;

    for (int i=0; i<static_cast<int>(bitIndices.size()); i++)
        if ((number >> bitIndices[i]) & QINDEX_ONE)
            out |= QINDEX_ONE << i;

    return out;
}


/**
 * TESTS
 *
 * @ingroup unitbitwise
 * @{
 */


TEST_CASE( "insertBits", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {
        vector<qindex> numbers = {0, 1, 2, 5, 21, 0x12345, 0x6DB6DB};
        vector<vector<int>> bitIndices = {
            {},
            {0},
            {1},
            {0, 1},
            {1, 3, 5},
            {0, 2, 6, 9},
            {4, 8, 12, 20}
        };

        for (auto number: numbers)
            for (auto& inds: bitIndices)
                for (int bitValue: {0, 1})
                    REQUIRE( insertBits(number, inds.data(), static_cast<int>(inds.size()), bitValue) == getRefInsertedBits(number, inds, bitValue) );
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


TEST_CASE( "getValueOfBits", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {
        vector<qindex> numbers = {0, 1, 2, 5, 0x12345, 0xAAAAAAAA, 0x55555555};
        vector<vector<int>> bitIndices = {
            {},
            {0},
            {1},
            {0, 1, 2, 3},
            {3, 1, 4, 0},
            {20, 0, 16, 8, 4}
        };

        for (auto number: numbers)
            for (auto& inds: bitIndices)
                REQUIRE( getValueOfBits(number, inds.data(), static_cast<int>(inds.size())) == getRefValueOfBits(number, inds) );
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


/** @} (end defgroup) */
