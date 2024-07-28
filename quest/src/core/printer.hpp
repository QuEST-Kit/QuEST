/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef PRINTER_HPP
#define PRINTER_HPP

#include "quest/include/types.h"
#include "quest/include/matrices.h"
#include "quest/include/paulis.h"

#include <vector>
#include <tuple>
#include <complex>
#include <string>
#include <sstream>
#include <type_traits>



/*
 * USER-CONFIGURABLE GLOBALS
 */

void printer_setMaxNumPrintedItems(qindex num);



/*
 * TYPE NAME STRINGS
 */

std::string printer_getQrealType();

std::string printer_getQcompType();

std::string printer_getQindexType();

std::string printer_getFloatPrecisionFlag();



/*
 * STRING CASTS
 */


// type T can be int, qindex, float, double, long double.
// this more specific printer_toStr() definition is called when passing complex
template <typename T>
std::string printer_toStr(std::complex<T> num);

// type T can be any non-complex type recognised by ostringstream!
template<typename T> 
std::string printer_toStr(T expr) {

    // write to buffer (rather than use to_string()) so that floating-point numbers
    // are automatically converted to scientific notation when necessary
    std::ostringstream buffer;
    buffer << expr;
    return buffer.str();
}



/*
 * SUBSTRINGS USED BY REPORTERS
 */

namespace printer_substrings { 
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
 * PRIMITIVE PRINTING
 */

void print(const char* str);
void print(std::string str);
void print(qcomp num);



/*
 * MATRIX PRINTING
 */

void print_matrixInfo(std::string nameStr, int numQubits, qindex dim, size_t elemMem, size_t otherMem, int numNodes=1);

void print_matrix(CompMatr1 matr, std::string indent=defaultMatrIndent);
void print_matrix(CompMatr2 matr, std::string indent=defaultMatrIndent);
void print_matrix(CompMatr  matr, std::string indent=defaultMatrIndent);
void print_matrix(DiagMatr1 matr, std::string indent=defaultMatrIndent);
void print_matrix(DiagMatr2 matr, std::string indent=defaultMatrIndent);
void print_matrix(DiagMatr  matr, std::string indent=defaultMatrIndent);
void print_matrix(FullStateDiagMatr matr, std::string indent=defaultMatrIndent);



/*
 * TABLE PRINTING
 */

void print_table(std::string title, std::vector<std::tuple<std::string, std::string>>   rows, std::string indent=defaultTableIndent);
void print_table(std::string title, std::vector<std::tuple<std::string, long long int>> rows, std::string indent=defaultTableIndent);



/*
 * PAULI PRINTING
 */

void print_pauliStr(PauliStr str, int numQubits);

void print_pauliStrSumInfo(qindex numTerms, qindex numBytes);

void print_pauliStrSum(PauliStrSum sum, std::string indent=defaultMatrIndent);



#endif // PRINTER_HPP