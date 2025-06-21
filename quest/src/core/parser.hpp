/** @file
 * Internal signatures for parsing user-given strings.
 * 
 * @author Tyson Jones
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include <string>

using std::string;



/*
 * PARSING NUMBERS
 */

bool parser_isAnySizedReal(string str);
bool parser_isAnySizedComplex(string str);

bool parser_isValidReal(string str);
bool parser_isValidComplex(string str);

qreal parser_parseReal(string str);
qcomp parser_parseComplex(string str);



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