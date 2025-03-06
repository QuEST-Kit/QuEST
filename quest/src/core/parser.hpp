/** @file
 * Internal signatures for parsing user-given strings.
 * 
 * @author Tyson Jones
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include "quest/include/paulis.h"

#include <string>

using std::string;



/*
 * PARSING INDIVIDUAL PAULIS
 */

const string parser_RECOGNISED_PAULI_CHARS = "0123ixyzIXYZ";

int parser_getPauliIntFromChar(char ch);



/*
 * PARSING PAULI STRING SUMS
 */

PauliStrSum parser_validateAndParsePauliStrSum(string lines, bool rightIsLeastSignificant, const char* caller);



/*
 * FILE IO
 */

bool parser_canReadFile(string fn);

string parser_loadFile(string fn);



#endif // PARSER_HPP