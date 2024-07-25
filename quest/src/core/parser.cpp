/** @file
 * Internal functions which parse strings given by users.
 * Performance here is not critical, so we opt to use
 * sub-optimal but clear, defensively functions.
 */

#include "quest/src/core/parser.hpp"
#include "quest/src/core/errors.hpp"
/*
 * PAULI STRING PARSING
 */


int parser_getPauliIntFromChar(char ch) {

    // must cover every char in parser_RECOGNISED_PAULI_CHARS
    switch (ch) {
        case '0': case 'i': case 'I': return 0;
        case '1': case 'x': case 'X': return 1;
        case '2': case 'y': case 'Y': return 2;
        case '3': case 'z': case 'Z': return 3;
    }

    // should be unreachable after validation
    error_attemptedToParseUnrecognisedPauliChar();
    return -1;
}



