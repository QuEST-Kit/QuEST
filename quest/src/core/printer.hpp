/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef PRINTER_HPP
#define PRINTER_HPP

#include "quest/include/types.h"
#include "quest/include/matrices.h"
#include "quest/include/paulis.h"

#include <tuple>
#include <vector>
#include <complex>
#include <string>
#include <sstream>
#include <type_traits>

using std::tuple;
using std::string;
using std::vector;
using std::complex;



/*
 * USER-CONFIGURABLE GLOBALS
 */

void printer_setMaxNumPrintedItems(qindex num);



/*
 * TYPE NAME STRINGS
 */

string printer_getQrealType();

string printer_getQcompType();

string printer_getQindexType();

string printer_getFloatPrecisionFlag();



/*
 * STRING CASTS
 */


// type T can be int, qindex, float, double, long double.
// this more specific printer_toStr() definition is called when passing complex
template <typename T>
string printer_toStr(complex<T> num);

// type T can be any non-complex type recognised by ostringstream!
template<typename T> 
string printer_toStr(T expr) {

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
    extern string eq;
    extern string mu;
    extern string by;
    extern string pn;
    extern string pg;
    extern string pm;
    extern string bt;
    extern string na;
    extern string un;
}



/*
 * DEFAULT INDENTS OF PRINTERS
 *
 * which doesn't need to be publicly accessible but we want them as
 * default arguments to the public signatures, so whatever. Declared
 * static to avoid symbol duplication by repeated header inclusion.
 */

static string defaultMatrIndent  = "    ";
static string defaultTableIndent = "  ";



/*
 * PRIMITIVE PRINTING
 */

void print(const char* str);
void print(string str);
void print(qcomp num);



/*
 * MATRIX PRINTING
 */

void print_matrixInfo(string nameStr, int numQubits, qindex dim, size_t elemMem, size_t otherMem, int numNodes=1);

void print_matrix(CompMatr1 matr, string indent=defaultMatrIndent);
void print_matrix(CompMatr2 matr, string indent=defaultMatrIndent);
void print_matrix(CompMatr  matr, string indent=defaultMatrIndent);
void print_matrix(DiagMatr1 matr, string indent=defaultMatrIndent);
void print_matrix(DiagMatr2 matr, string indent=defaultMatrIndent);
void print_matrix(DiagMatr  matr, string indent=defaultMatrIndent);
void print_matrix(FullStateDiagMatr matr, string indent=defaultMatrIndent);



/*
 * TABLE PRINTING
 */

void print_table(string title, vector<tuple<string, string>>   rows, string indent=defaultTableIndent);
void print_table(string title, vector<tuple<string, long long int>> rows, string indent=defaultTableIndent);



/*
 * PAULI PRINTING
 */

void print_pauliStr(PauliStr str, int numQubits);

void print_pauliStrSumInfo(qindex numTerms, qindex numBytes);

void print_pauliStrSum(PauliStrSum sum, string indent=defaultMatrIndent);



#endif // PRINTER_HPP