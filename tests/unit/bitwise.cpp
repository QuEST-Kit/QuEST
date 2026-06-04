/** @file
 * Unit tests of internal bitwise helpers.
 *
 * @defgroup unitbitwise Bitwise
 * @ingroup unittests
 */

#include "quest/src/core/bitwise.hpp"

#include <catch2/catch_test_macros.hpp>

#include "tests/utils/macros.hpp"



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[bitwise]"


static qindex getReferenceInsertBits(qindex number, const int* bitIndices, int numIndices, int bitValue) {

    for (int i=0; i<numIndices; i++)
        number = insertBit(number, bitIndices[i], bitValue);

    return number;
}


static qindex getReferenceValueOfBits(qindex number, const int* bitIndices, int numIndices) {

    qindex value = 0;

    for (int i=0; i<numIndices; i++)
        value |= getBit(number, bitIndices[i]) << i;

    return value;
}



/**
 * TESTS
 *
 * @ingroup unitbitwise
 * @{
 */


TEST_CASE( "insertBits", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int bitInds[] = {1, 3, 6, 9, 12, 20};
        int numInds = 6;
        qindex number = 0b101101001011;

        REQUIRE( insertBits(number, bitInds, numInds, 0) == getReferenceInsertBits(number, bitInds, numInds, 0) );
        REQUIRE( insertBits(number, bitInds, numInds, 1) == getReferenceInsertBits(number, bitInds, numInds, 1) );
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


TEST_CASE( "getValueOfBits", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qindex number = 0b101101101001011;

        SECTION( "increasing indices" ) {

            int bitInds[] = {0, 2, 5, 8, 11, 14};
            int numInds = 6;
            REQUIRE( getValueOfBits(number, bitInds, numInds) == getReferenceValueOfBits(number, bitInds, numInds) );
        }

        SECTION( "arbitrarily ordered indices" ) {

            int bitInds[] = {14, 0, 8, 2, 11, 5};
            int numInds = 6;
            REQUIRE( getValueOfBits(number, bitInds, numInds) == getReferenceValueOfBits(number, bitInds, numInds) );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


TEST_CASE( "insertBitsWithMaskedValues", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        int bitInds[] = {2, 4, 7, 10, 13, 21};
        int numInds = 6;
        qindex number = 0b1101001011;
        qindex mask = 0;

        for (int i=0; i<numInds; i++)
            mask |= QINDEX_ONE << bitInds[i];

        qindex expected = mask | getReferenceInsertBits(number, bitInds, numInds, 0);
        REQUIRE( insertBitsWithMaskedValues(number, bitInds, numInds, mask) == expected );
    }

    SECTION( LABEL_VALIDATION ) {

        // no validation!
        SUCCEED( );
    }
}


/** @} (end defgroup) */
