/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef FORMATTER_HPP
#define FORMATTER_HPP

#include "quest/include/types.h"
#include "quest/include/matrices.h"

#include <vector>
#include <tuple>
#include <complex>
#include <string>
#include <sstream>
#include <type_traits>



/*
 * USER-CONFIGURABLE GLOBALS
 */

void form_setMaxNumPrintedMatrixElems(qindex num);



/*
 * TYPE NAME STRINGS
 */

std::string form_getQrealType();

std::string form_getQcompType();

std::string form_getQindexType();

std::string form_getFloatPrecisionFlag();



/*
 * STRING CASTS
 */


// type T can be int, qindex, float, double, long double.
// this more specific form_str() definition is called when passing complex
template <typename T>
std::string form_str(std::complex<T> num);

// type T can be any non-complex type recognised by ostringstream!
template<typename T> 
std::string form_str(T expr) {

    // write to buffer (rather than use to_string()) so that floating-point numbers
    // are automatically converted to scientific notation when necessary
    std::ostringstream buffer;
    buffer << expr;
    return buffer.str();
}



/*
 * SUBSTRINGS USED BY REPORTERS
 */

namespace form_substrings { 
    extern std::string eq;
    extern std::string mu;
    extern std::string by;
    extern std::string pn;
    extern std::string pg;
    extern std::string pm;
    extern std::string bt;
    extern std::string na;
    extern std::string un;
}



/*
 * DEFAULT INDENTS OF PRINTERS
 *
 * which doesn't need to be publicly accessible but we want them as
 * default arguments to the public signatures, so whatever. Declared
 * static to avoid symbol duplication by repeated header inclusion.
 */

static std::string defaultMatrIndent  = "    ";
static std::string defaultTableIndent = "  ";



/*
 * MATRIX PRINTING
 */

void form_printMatrixInfo(std::string nameStr, int numQubits, qindex dim, size_t elemMem, size_t otherMem, int numNodes=1);

void form_printMatrix(CompMatr1 matr, std::string indent=defaultMatrIndent);
void form_printMatrix(CompMatr2 matr, std::string indent=defaultMatrIndent);
void form_printMatrix(CompMatr  matr, std::string indent=defaultMatrIndent);
void form_printMatrix(DiagMatr1 matr, std::string indent=defaultMatrIndent);
void form_printMatrix(DiagMatr2 matr, std::string indent=defaultMatrIndent);
void form_printMatrix(DiagMatr  matr, std::string indent=defaultMatrIndent);
void form_printMatrix(FullStateDiagMatr matr, std::string indent=defaultMatrIndent);


/*
 * TABLE PRINTING
 */

void form_printTable(std::string title, std::vector<std::tuple<std::string, std::string>>   rows, std::string indent=defaultTableIndent);
void form_printTable(std::string title, std::vector<std::tuple<std::string, long long int>> rows, std::string indent=defaultTableIndent);



#endif // FORMATTER_HPP