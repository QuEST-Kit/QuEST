/** @file
 * Functions for formatting and outputting strings, used
 * primarily by reporters (e.g. reportQureg). 
 * 
 * @author Tyson Jones
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

// beware that files including this header receive all these
// namespace items; a worthwhile evil to keep this readable
using std::tuple;
using std::string;
using std::vector;
using std::complex;



/*
 * USER-CONFIGURABLE GLOBALS
 */

void printer_setMaxNumPrintedScalars(qindex numRows, qindex numCols);

void printer_setMaxNumPrintedSigFig(int numSigFigs);

void printer_setNumTrailingNewlines(int numNewlines);

int printer_getNumTrailingNewlines();

void printer_setPauliChars(string newChars);

void printer_setPauliStrFormat(int flag);



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
    // are automatically converted to scientific notation when necessary. Beware
    // that the user configured significant figures are not reflected here.

    std::ostringstream buffer;
    buffer << expr;
    return buffer.str();
}

// explicit qreal version of above, affected by user-set significant figures
string printer_toStr(qreal num);



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
 * TRAILING NEWLINE PRINTING
 */

void print_newlines();
void print_oneFewerNewlines();



/*
 * PRIMITIVE PRINTING
 */

// never prints any newlines (except within supplied str)
void print(const char* str);
void print(string str);
void print(qcomp num);



/*
 * STRUCT HEADER PRINTING
 */

// always prints "label:\n"
void print_label(string label);

// always prints a single trailing newline
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

// always prints a SINGLE trailing newline, so never supports run-on
// of multiple matrices/multi-line prints; caller must ergo validate
// that user's num-newlines > 0 (else we violatd it) and add lines
// if necessary. Note the default indentation is almost never overriden
void print_elems(CompMatr1,         string indent=defaultMatrIndent);
void print_elems(CompMatr2,         string indent=defaultMatrIndent);
void print_elems(CompMatr ,         string indent=defaultMatrIndent);
void print_elems(DiagMatr1,         string indent=defaultMatrIndent);
void print_elems(DiagMatr2,         string indent=defaultMatrIndent);
void print_elems(DiagMatr,          string indent=defaultMatrIndent);
void print_elems(KrausMap,          string indent=defaultMatrIndent);
void print_elems(SuperOp,           string indent=defaultMatrIndent);
void print_elems(Qureg,             string indent=defaultMatrIndent);
void print_elems(PauliStrSum,       string indent=defaultMatrIndent);
void print_elems(FullStateDiagMatr, string indent=defaultMatrIndent);

// always prints NO trailing newlines
void print_elemsWithoutNewline(PauliStr, string indent="");



/*
 * TABLE PRINTING
 */

void print_table(string title, vector<tuple<string, string>>        rows, string indent=defaultTableIndent);
void print_table(string title, vector<tuple<string, long long int>> rows, string indent=defaultTableIndent);
void print_table(string title, string emptyPlaceholder,                   string indent=defaultTableIndent);



#endif // PRINTER_HPP