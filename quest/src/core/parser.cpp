/** @file
 * Internal functions which parse strings given by users.
 * 
 * Performance here is not critical, so we opt to use
 * sub-optimal but clear, defensively-designed functions.
 * We also use quite stringent internal error checking
 * due to the many risks and uncertainties caused by 
 * parsing user-given strings with regex.
 * 
 * @author Tyson Jones
 */

#include "quest/include/config.h"
#include "quest/include/precision.h"
#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/parser.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/validation.hpp"

#include <regex>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

using std::regex;
using std::vector;
using std::string;
using std::smatch;
using std::ifstream;
using std::out_of_range;
using std::stringstream;
using std::sregex_iterator;
using std::invalid_argument;



/*
 * REGULAR EXPRESSIONS
 *
 * which capture our supported PauliStrSum string syntax
 */


namespace patterns {

    // utilities
    string group(string in) { return "(?:" + in + ")"; } // groups sub-patterns to control eval order
    string capt (string in) { return "("   + in + ")"; } // captures sub-patterns for later extraction
    string opt  (string in) { return group(in) + "?"; }  // groups and makes optional a sub-pattern

    // sub-numbers
    string mantissa = "[0-9]+" + opt("[.][0-9]*") + "|[.][0-9]+"; // e.g.  0  0.  .0  0.0
    string exponent = "[eE][+-]?[0-9]+";                          // e.g.  e5  e-5  E5  E-5  e-0
    string imagsymb = "[ijIJ]";

    // unsigned
    string ureal = group(mantissa) + opt(exponent); // e.g.  4  5E0  .1E-10  

    // constants we may wish to generalise
    string space = "[ \\t]";  // single horizontal whitespace char
    string sign  = "[+-]";

    // optional
    string optSpace = space + "*";
    string optSign  = optSpace + opt(sign) + optSpace; // optional +- with any spacing
    
    // component delimiter
    string delim = optSpace + sign + optSpace; // mandatory +- with any spacing

    // full complex; only one component given, which is captured with sign (excluding imag symb)
    string real = capt(optSign + ureal)            + optSpace;
    string imag = capt(optSign + ureal) + imagsymb + optSpace;

    // full complex; both components given, real always before imag, individually captured with signs (but without imag symb)
    string comp = capt(optSign + ureal) + capt(delim + ureal) + imagsymb + optSpace;

    // full complex; any format, importantly in order of decreasing specificity. do not consult for captured groups
    string num = group(comp) + "|" + group(imag) + "|" + group(real);

    // no capturing because 'num' pollutes captured groups, and pauli syntax overlaps real integers
    string pauli = "[" + parser_RECOGNISED_PAULI_CHARS + "]";
    string paulis = group(optSpace + pauli + optSpace) + "+";
    string weightedPaulis = "^" + group(num) + space + optSpace + paulis + "$";
}


namespace regexes {

    // instantiation here negates risk of later runtime error due to invalid regex
    regex real(patterns::real);
    regex imag(patterns::imag);
    regex comp(patterns::comp);
    regex num(patterns::num);
    regex paulis(patterns::paulis);
    regex weightedPaulis(patterns::weightedPaulis);
}



/*
 * HANDLING WHITESPACE
 *
 * in a verbose, centralised way to ensure we never tolerate space characters
 * in our regex which we do not later permit in our parsing. As such,
 * isWhiteSpace() must return 'true' for every character in patterns:space
 */


bool isWhiteSpace(char ch) {

    // matches patterns:space and more (additionally; newlines)
    return isspace(ch); 
}

bool isNotWhiteSpace(char ch) {

    return !isWhiteSpace(ch);
}

bool isOnlyWhiteSpace(string str) {

    return all_of(str.begin(), str.end(), isWhiteSpace);
}


void removeWhiteSpace(string &line) {

    // modifies line var, including removing newlines. we don't need to do
    // this for pattern matching (our regex patterns include all permitted,
    // frivolous whitespace), but we do need to perform it before passing
    // matched strings to functions like std::stold()
    line.erase(remove_if(line.begin(), line.end(), isWhiteSpace), line.end());
}



/*
 * STRING PARTITIONING
 */


void separateStringIntoCoeffAndPaulis(string line, string &coeff, string &paulis) {

    // absolutely gauranteed to match due to prior matching of regexes::line
    smatch match;

    // locate the coefficient
    regex_search(line, match, regexes::num);
    coeff = match.str(0);

    // the remainder of the line must be the paulis, but we explicitly match to
    // regex just in case we later permit additional substrings like comments. 
    // we must match only the substring AFTER the coeff, since valid coeffs (e.g.
    // 1) can be mistaken for a Pauli code.
    line = line.substr(match.position(0) + match.length(0));
    regex_search(line, match, regexes::paulis);
    paulis = match.str(0);
}


int getNumPaulisInLine(string line) {

    // simply count the non-whitespace chars in the paulis substring
    string coeff, paulis;
    separateStringIntoCoeffAndPaulis(line, coeff, paulis);
    return count_if(paulis.begin(), paulis.end(), isNotWhiteSpace);
}



/*
 * REAL NUMBER PARSING
 */


qreal precisionAgnosticStringToFloat(string str) {

    // remove whitespace which stold() et al cannot handle after the sign.
    // beware this means that e.g. "1 0" (invalid number) would become "10"
    // (valid) so this function cannot be used for duck-typing, though that
    // is anyway the case since stold() et al permit "10abc"
    removeWhiteSpace(str);

    // below throws exception when the (prefix) of str cannot be/fit into a qreal
    if (FLOAT_PRECISION == 1) return static_cast<qreal>(std::stof (str));
    if (FLOAT_PRECISION == 2) return static_cast<qreal>(std::stod (str));
    if (FLOAT_PRECISION == 4) return static_cast<qreal>(std::stold(str));

    // unreachable
    return -1;
}


bool parser_isAnySizedReal(string str) {

    // we assume that all strings which match the regex can be parsed by
    // precisionAgnosticStringToFloat() above (once whitespace is removed)
    // EXCEPT strings which contain a number too large to store in the qreal 
    // type (as is separately checked below). Note it is insufficient to merely 
    // duck-type using stold() et al because such functions permit non-numerical 
    // characters to follow the contained number which are silently removed (grr!)
    smatch match;
    return regex_match(str, match, regexes::real);
}


bool parser_isValidReal(string str) {

    // reject str if it doesn't match regex
    if (!parser_isAnySizedReal(str))
        return false;

    // check number is in-range of qreal via duck-typing
    try {
        precisionAgnosticStringToFloat(str);
    } catch (const out_of_range&) {
        return false;

    // error if our regex permitted an unparsable string
    } catch (const invalid_argument&) {
        error_attemptedToParseRealFromInvalidString();
    }

    return true;
}


qreal parser_parseReal(string str) {

    try {
        return precisionAgnosticStringToFloat(str);
    } catch (const invalid_argument&) {
        error_attemptedToParseRealFromInvalidString();
    } catch (const out_of_range&) {
        error_attemptedToParseOutOfRangeReal();
    }

    // unreachable
    return -1;
}



/*
 * COMPLEX NUMBER PARSING
 */


bool parser_isAnySizedComplex(string str) {

    // we assume that all strings which match the regex can be parsed to
    // a qcomp (once whitespace is removed) EXCEPT strings which contain a 
    // number too large to store in the qcomp type (as is separately checked 
    // below). Note it is insufficient to merely duck-type each component using
    // using stold() et al because such functions permit non-numerical chars to 
    // follow the contained number (grr!)
    smatch match;

    // must match real, imaginary or complex number regex
    if (regex_match(str, match, regexes::real)) return true;
    if (regex_match(str, match, regexes::imag)) return true;
    if (regex_match(str, match, regexes::comp)) return true;

    return false;
}


bool parser_isValidComplex(string str) {

    // reject str if it doesn't match complex regex
    if (!parser_isAnySizedComplex(str))
        return false;

    // we've so far gauranteed str has a valid form, but we must now check 
    // each included complex component (which we enumerate) is in range of a qreal
    sregex_iterator it(str.begin(), str.end(), regexes::real);
    sregex_iterator end;

    // valid coeffs contain 1 or 2 reals, never 0, which regex should have caught
    if (it == end)
        error_attemptedToParseComplexFromInvalidString();

    // for each of the 1 or 2 components...
    for (; it != end; it++) {

        // check component is in-range of qreal via duck-typing
        try {
            precisionAgnosticStringToFloat(it->str(0));
        } catch (const out_of_range&) {
            return false;

        // error if our regex permitted an unparsable component
        } catch (const invalid_argument&) {
            error_attemptedToParseComplexFromInvalidString();
        }
    }

    // report that each/all detected components of str can form a valid qcomp
    return true;
}


qcomp parser_parseComplex(string str) {

    if (!parser_isValidComplex(str))
        error_attemptedToParseComplexFromInvalidString();

    // we are gauranteed to fully match real, imag or comp after prior validation
    smatch match;

    // extract and parse components and their signs (excluding imaginary symbol)
    if (regex_match(str, match, regexes::real))
        return qcomp(parser_parseReal(match.str(1)), 0);

    if (regex_match(str, match, regexes::imag))
        return qcomp(0, parser_parseReal(match.str(1)));

    if (regex_match(str, match, regexes::comp))
        return qcomp(
            parser_parseReal(match.str(1)),
            parser_parseReal(match.str(2)));
    
    // should be unreachable
    error_attemptedToParseComplexFromInvalidString();
    return qcomp(0,0); 
}



/*
 * VALIDATION
 *
 * which checks user-given strings are correctly formatted, which are
 * defined here (in lieu of inside validation.cpp) because they make
 * extensive use of regex and parsing.
 */


bool isInterpretablePauliStrSumLine(string line) {

    // checks whether line has the expected format; a real, imaginary or
    // complex number (expressed as an integer, decimal, or in scientific 
    // notation) followed by 1 or more space characters, then one or
    // more pauli codes/chars. It does NOT determine whether the coeff
    // can actually be instantiated as a qcomp
    return regex_match(line, regexes::weightedPaulis);
}


bool isPauliStrSumCoeffWithinQcompRange(string line) {

    // it is gauranteed that line is interpretable and contains a regex-matching
    // coefficient, but we must additionally verify it is within range of qreal.
    // So we duck type each of the 1 or 2 matches with the real regex (i.e. one or 
    // both of the real and imaginary components of a complex coeff).

    // process only the coeff, since pauli codes may resemble coeffs (e.g. 1)
    string coeff, _; 
    separateStringIntoCoeffAndPaulis(line, coeff, _); // discard pauli substr

    // beautiful iterator boilerplate
    sregex_iterator it(coeff.begin(), coeff.end(), regexes::real);
    sregex_iterator end;

    // valid coeffs contain 1 or 2 reals, never 0
    if (it == end)
        return false;

    // enumerate all matches of 'real' regex in line
    for (; it != end; it++) {

        // remove whitespace (stold cannot handle space between sign and number)
        string match = it->str(0);
        removeWhiteSpace(match);

        // return false if number cannot become a qreal
        try {
            precisionAgnosticStringToFloat(match);
        } catch (const out_of_range&) {
            return false;
        } catch (const invalid_argument&) { // should be impossible (indicates bad regex)
            return false;
        }
    }

    // report that each double in coeff substr was parsable
    return true;
}


void assertStringIsValidPauliStrSum(string lines, const char* caller) {

    int numPaulis = 0;
    bool nonEmptyLine = false;
    qindex lineIndex = 0;

    // parse each line in-turn
    stringstream stream(lines);

    for (string line; getline(stream, line); lineIndex++) {

        // permit and skip empty lines
        if (isOnlyWhiteSpace(line))
            continue;
        else
            nonEmptyLine = true;
        
        // assert the line is interpretable
        bool validLine = isInterpretablePauliStrSumLine(line);
        validate_parsedPauliStrSumLineIsInterpretable(validLine, line, lineIndex, caller);

        // assert the coeff is parsable (e.g. doesn't exceed valid number range)
        bool validCoeff = isPauliStrSumCoeffWithinQcompRange(line);
        validate_parsedPauliStrSumCoeffWithinQcompRange(validCoeff, line, lineIndex, caller);

        // assert the line has a consistent number of Paulis as previous
        int numLinePaulis = getNumPaulisInLine(line);
        if (!numPaulis)
            numPaulis = numLinePaulis;
        validate_parsedPauliStrSumLineHasConsistentNumPaulis(numPaulis, numLinePaulis, line, lineIndex, caller);
    }

    // ensure we parsed at least 1 Pauli string
    validate_parsedStringIsNotEmpty(nonEmptyLine, caller);
}



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



/*
 * PAULI STRING SUM PARSING
 */


PauliStr parsePaulis(string paulis, bool rightIsLeastSignificant) {

    // remove whitespace to make string compatible with getPauliStr()
    removeWhiteSpace(paulis);

    // default creator treats rightmost pauli as least significant
    if (rightIsLeastSignificant)
        return getPauliStr(paulis);

    // otherwise pass qubits {0, 1, ... }
    vector<int> qubits(paulis.size());
    for (size_t i=0; i<paulis.size(); i++)
        qubits[i] = i;

    return getPauliStr(paulis, qubits);
}


void parseWeightedPaulis(string line, qcomp &coeff, PauliStr &pauli, bool rightIsLeastSignificant) {

    // separate line into substrings
    string coeffStr, pauliStr;
    separateStringIntoCoeffAndPaulis(line, coeffStr, pauliStr);

    // parse each, overwriting calller primitives
    coeff = parser_parseComplex(coeffStr);
    pauli = parsePaulis(pauliStr, rightIsLeastSignificant);
}


qindex getNumLines(string lines) {

    /// @todo is this platform agnostic?
    char newline = '\n';

    return 1 + count(lines.begin(), lines.end(), newline);
}


PauliStrSum parser_validateAndParsePauliStrSum(string lines, bool rightIsLeastSignificant, const char* caller) {
    assertStringIsValidPauliStrSum(lines, caller);

    // allocate space for as many strings as there are lines, though we might collect fewer (we skip empty lines)
    qindex numLines = getNumLines(lines);
    vector<qcomp>    coeffs;   coeffs.reserve(numLines);
    vector<PauliStr> strings; strings.reserve(numLines);
    
    // parse each line in-turn
    stringstream stream(lines);
    for (string line; getline(stream, line); ) {

        // skip empty lines
        if (isOnlyWhiteSpace(line))
            continue;

        qcomp coeff;
        PauliStr string;
        parseWeightedPaulis(line, coeff, string, rightIsLeastSignificant); // validates

        coeffs.push_back(coeff);
        strings.push_back(string);
    }

    // invoke other API function, which will run additional validation (e.g. checking
    // memory allocations succeeded) and may fail, reporting that 'createPauliStrSum()'
    // failed, rather than caller. That's a minor, acceptable evil. Furthermore, this
    // call creates new memory alongside our existing vectors, meaning total memory
    // is temporarily double than that which is strictly necessary; also acceptable.
    PauliStrSum out = createPauliStrSum(strings, coeffs);
    return out;
}



/*
 * FILE IO
 */


bool parser_canReadFile(string fn) {

    ifstream file(fn);
    return file.good();
}


string parser_loadFile(string fn) {

    // ensure file is still readable since validation
    ifstream file(fn);
    if (!file.good())
        error_couldNotReadFile();

    // load entire file into string
    stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}
