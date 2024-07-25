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



#endif // PARSER_HPP