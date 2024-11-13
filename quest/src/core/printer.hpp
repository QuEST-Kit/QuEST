/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef PRINTER_HPP
#define PRINTER_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"
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

void printer_setMaxNumPrintedScalars(qindex numRows, qindex numCols);



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
 * SUBSTRING PREPARATION
 */

string printer_getMemoryWithUnitStr(size_t numBytes);



/*
 * SUBSTRINGS USED BY REPORTERS
 *
 * which are padded with 1 whitespace either side
 */

namespace printer_substrings { 
    extern string eq; // =
    extern string pn; // per node
    extern string pg; // per gpu
    extern string ig; // in gpu
    extern string pm; // per machine
    extern string bt; // 2^
    extern string na; // N/A
    extern string un; // unknown
}



/*
 * DEFAULT INDENTS OF PRINTERS
 *
 * which doesn't need to be publicly accessible but we want them as
 * default arguments to the public signatures, so whatever. Declared
 * static to avoid symbol duplication by repeated header inclusion.
 */

static string defaultMatrIndent  = "    ";
static string defaultQuregIndent = "    ";
static string defaultKrausIndent = "  ";
static string defaultTableIndent = "  ";



/*
 * PRIMITIVE PRINTING
 */

void print(const char* str);
void print(string str);
void print(qcomp num);



/*
 * STRUCT HEADER PRINTING
 */

void print_header(CompMatr1 m,  size_t numBytes);
void print_header(CompMatr2 m,  size_t numBytes);
void print_header(CompMatr  m,  size_t numBytes);
void print_header(DiagMatr1 m,  size_t numBytes);
void print_header(DiagMatr2 m,  size_t numBytes);
void print_header(DiagMatr  m,  size_t numBytes);
void print_header(SuperOp  op,  size_t numBytes);
void print_header(KrausMap map, size_t numBytes);
void print_header(PauliStrSum sum, size_t numBytes);

void print_header(FullStateDiagMatr m, size_t numBytesPerNode);
void print_header(Qureg qureg, size_t numBytesPerNode);



/*
 * STRUCT ELEMENT PRINTING
 */

void print_elems(CompMatr1 matr, string indent=defaultMatrIndent);
void print_elems(CompMatr2 matr, string indent=defaultMatrIndent);
void print_elems(CompMatr  matr, string indent=defaultMatrIndent);
void print_elems(DiagMatr1 matr, string indent=defaultMatrIndent);
void print_elems(DiagMatr2 matr, string indent=defaultMatrIndent);
void print_elems(DiagMatr  matr, string indent=defaultMatrIndent);
void print_elems(KrausMap map,   string indent=defaultMatrIndent);
void print_elems(SuperOp op,     string indent=defaultMatrIndent);
void print_elems(Qureg qureg,    string indent=defaultMatrIndent);
void print_elems(PauliStrSum sum, string indent=defaultMatrIndent);
void print_elems(PauliStr str, int numQubits);
void print_elems(FullStateDiagMatr qureg, string indent=defaultMatrIndent);



/*
 * TABLE PRINTING
 */

void print_table(string title, vector<tuple<string, string       >> rows, string indent=defaultTableIndent);
void print_table(string title, vector<tuple<string, long long int>> rows, string indent=defaultTableIndent);



#endif // PRINTER_HPP