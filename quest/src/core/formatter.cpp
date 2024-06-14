/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#include "quest/include/types.h"
#include "quest/include/structures.h"

#include "quest/src/core/formatter.hpp"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <tuple>
#include <string>


// I like to live dangerously; and I am willing to debate with the
// most dogmatic C++ enthusiast about the unreadability of namespaces
using namespace std;


// aesthetic constants
const int MIN_SPACE_BETWEEN_MATRIX_COLS = 2;
const int MIN_SPACE_BETWEEN_TABLE_COLS = 5;

const char TABLE_SPACE_CHAR = '.';
const char MATRIX_SPACE_CHAR = ' ';



/*
 * TYPE NAME STRINGS
 */


// try safely check for de-mangler library
#if defined __has_include
    #if __has_include (<cxxabi.h>)
        #include <cxxabi.h>
        #define DEMANGLE_TYPE true
    #else
        #define DEMANGLE_TYPE false
    #endif

// fall-back to dangerous compiler-specific check (we know MSVC doesn't mangle)
#elif ! defined _MSC_VER

    #warning "Attempting to include compiler-specific library <cxxabi.h>"

    #include <cxxabi.h>
    #define DEMANGLE_TYPE true
#else
    #define DEMANGLE_TYPE false
#endif


// macros for printing multi-word macros (necessary indirection)
#define GET_STR_INTERNAL(x) #x
#define GET_STR(x) GET_STR_INTERNAL(x)


template<typename T> 
string getTypeName(T) { 

    // get possibly-mangled type name
    string name = typeid(T).name();

    // non-MSVC compilers can de-mangle
    #ifdef DEMANGLE_TYPE
        name = abi::__cxa_demangle(name.c_str(), nullptr, nullptr, nullptr);
    #endif

    return name;
}


string form_getQcompType() {

    qcomp x;
    return getTypeName(x);
}

string form_getQrealType() {

    // more portable than getTypeName
    return GET_STR( FLOAT_TYPE );
}

string form_getQindexType() {

    // more portable than getTypeName
    return GET_STR( INDEX_TYPE );
}

string form_getFloatPrecisionFlag() {

    return GET_STR( FLOAT_PRECISION );
}



/*
 * COMPLEX STRINGS
 *
 * which can be any precision, and ergo different to that of qcomp
 */

template <typename T>
string floatToStr(T num, bool hideSign=false) {

    // write to stream (instead of calling to_string()) to auto-use scientific notation
    ostringstream buffer;

    // ensure +- is not shown if forced
    if (hideSign) {
        buffer << noshowpos;
        if (num < 0)
            num *= -1;
    }

    // return 0, 1, .1, 1.2e-5
    buffer << num;
    return buffer.str();
}

template <typename T>
string compToStr(complex<T> num) {

    // precise 0 is rendered as a real integer
    if (real(num) == 0 && imag(num) == 0)
        return "0";

    // real complex-floats neglect imag component
    if (imag(num) == 0)
        return floatToStr(real(num));

    // imaginary complex-floats neglect real component
    string realStr = (real(num) == 0)? "" : floatToStr(real(num));

    // -1i is abbreviated to -i
    if (imag(num) == -1)
        return realStr + "-i";

    // +1i is abbreviated to +i, although the sign is dropped if there's no preceeding real component
    if (imag(num) == 1)
        return realStr + ((real(num) == 0)? "i" : "+i");

    // get imag component string but without sign
    string imagStr = floatToStr(imag(num), true);

    // scientific-notation components always get wrapped in paranthesis
    if (imagStr.find('e') != string::npos)
        imagStr = '(' + imagStr + ')';

    // negative imag components always get - sign prefix
    if (imag(num) < 0)
        imagStr = '-' + imagStr;

    // positive imag components only get + sign prefix if they follow a real component
    if (imag(num) > 0 && real(num) != 0)
        imagStr = '+' + imagStr;

    // all imag components end in 'i'
    imagStr += 'i';

    return realStr + imagStr;
}



/*
 * MATRIX PRINTING
 */


template <class T> 
void printMatrixElems(T matr, string indent) {

    // determine max width of each column
    vector<int> maxColWidths(matr.numRows, 0);
    for (qindex c=0; c<matr.numRows; c++) {
        for (qindex r=0; r<matr.numRows; r++) {
            int width = compToStr(matr.elems[r][c]).length();
            if (width > maxColWidths[c])
                maxColWidths[c] = width;
        }
    }

    // print each row, aligning columns
    for (qindex r=0; r<matr.numRows; r++) {
        cout << indent;

        for (qindex c=0; c<matr.numRows; c++)
            cout << left
                << setw(maxColWidths[c] + MIN_SPACE_BETWEEN_MATRIX_COLS) 
                << setfill(MATRIX_SPACE_CHAR)
                << compToStr(matr.elems[r][c]);

        cout << endl;
    }
}


void form_printMatrix(CompMatr1 matr, string indent) {

    printMatrixElems(matr, indent);
}

void form_printMatrix(CompMatr2 matr, string indent) {

    printMatrixElems(matr, indent);
}

void form_printMatrix(CompMatrN matr, string indent) {

    printMatrixElems(matr, indent);
}



/*
 * SUBSTRINGS USED BY REPORTERS
 */

namespace form_substrings { 
    string eq = " = ";
    string mu = " x ";
    string by = " bytes";
    string pn = " per node";
    string pg = " per gpu";
    string pm = " per machine";
    string bt = "2^";
    string na = "N/A";
    string un = "unknown";
}



/*
 * TABLE PRINTING
 */


void form_printTable(string title, vector<tuple<string, string>> rows, string indent) {

    // find max-width of left column
    int maxWidth = 0;
    for (auto const& [key, value] : rows)
        if (key.length() > maxWidth)
            maxWidth = key.length();

    // print table title (indented)
    cout << indent << "[" << title << "]" << endl;

    // pad left column (which is doubly indented) to align right column
    for (auto const& [key, value] : rows)
        cout 
            << indent << indent << left 
            << setw(maxWidth + MIN_SPACE_BETWEEN_TABLE_COLS) << setfill(TABLE_SPACE_CHAR) 
            << key << value << endl;
}


void form_printTable(string title, vector<tuple<string, long long int>> rows, string indent) {

    // convert all values to strings
    vector<tuple<string, string>> casted;
    for (auto const& [key, value] : rows)
        casted.push_back({key, to_string(value)});

    form_printTable(title, casted, indent);
}
