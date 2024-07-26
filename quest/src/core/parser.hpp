/** @file
 * Internal signatures for parsing user-given strings.
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include "quest/include/paulis.h"

#include <string>



/*
 * PARSING INDIVIDUAL PAULIS
 */

const std::string parser_RECOGNISED_PAULI_CHARS = "0123ixyzIXYZ";

int parser_getPauliIntFromChar(char ch);



/*
 * PARSING PAULI STRING SUMS
 */

PauliStrSum parser_validateAndParsePauliStrSum(std::string lines, bool rightIsLeastSignificant, const char* caller);



/*
 * FILE IO
 */

bool parser_canReadFile(std::string fn);

std::string parser_loadFile(std::string fn);



#endif // PARSER_HPP