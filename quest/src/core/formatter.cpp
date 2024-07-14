/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#include "quest/include/types.h"
#include "quest/include/matrices.h"

#include "quest/src/core/formatter.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"

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
using namespace form_substrings;


// aesthetic constants
const int MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS = 2;
const int MIN_SPACE_BETWEEN_DIAG_MATRIX_COLS = 1;
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


// type T can be anything in principle, although it's currently only used for qcomp
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
 * NUMBER STRINGIFYING
 *
 * which is precision and type agnostic - i.e. we can stringify non-qcomp
 * and non-qreal numbers given by users or the backend
 */


// type T can be any number, e.g. int, qindex, float, double, long double
template <typename T>
string floatToStr(T num, bool hideSign=false) {

    // write to stream (instead of calling to_string()) to auto-use scientific notation
    ostringstream buffer;

    // ensure +- is not shown if forced hidden
    if (hideSign) {
        buffer << noshowpos;
        if (num < 0)
            num *= -1;
    }

    // return 0, 1, .1, 1.2e-5
    buffer << num;
    return buffer.str();
}


// type T can be precision decimal (independent of qreal) i.e. float, double, long double
template <typename T>
string form_str(complex<T> num) {

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


// explicitly instantiate all publicly passable types
template string form_str<int>(complex<int> num);
template string form_str<long int>(complex<long int> num);
template string form_str<long long int>(complex<long long int> num);

template string form_str<unsigned>(complex<unsigned> num);
template string form_str<long unsigned>(complex<long unsigned> num);
template string form_str<long long unsigned>(complex<long long unsigned> num);

template string form_str<float>(complex<float> num);
template string form_str<double>(complex<double> num);
template string form_str<long double>(complex<long double> num);



/*
 * MATRIX INFO PRINTING
 *
 * which prints information about the matrix type itself, rather than its elements
 */


void form_printMatrixInfo(std::string nameStr, int numQubits, qindex dim, size_t elemMem, size_t otherMem, int numNodes) {

    // only root node prints
    if (comm_getRank() != 0)
        return;

    // please don't tell anyone how I live
    bool isDiag = nameStr.find("Diag") != std::string::npos;

    // prepare substrings
    std::string sepStr  = ", ";
    std::string qbStr   = form_str(numQubits) + " qubit" + ((numQubits>1)? "s":"");
    std::string dimStr  = form_str(dim) + (isDiag? "" : mu + form_str(dim)) + " qcomps";
    std::string memStr  = form_str(elemMem) + " + " + form_str(otherMem) + by;

    if (numNodes > 1) {
        dimStr += " over " + form_str(numNodes) + " nodes";
        memStr += " per node";
    }

    // print e.g. FullStateDiagMatr (8 qubits, 256 qcomps over 4 nodes, 1024 + 32 bytes per node):
    std::cout 
        << nameStr << " (" 
        << qbStr   << sepStr 
        << dimStr  << sepStr 
        << memStr  << "):" 
        << std::endl;
}



/*
 * DENSE MATRIX PRINTING
 *
 * which print only the matrix elements, indented
 */


// type T can be qcomp(*)[] or qcomp**
template <typename T> 
void printDenseMatrixElems(T elems, qindex dim, string indent) {
    
    // determine max width of each column
    vector<int> maxColWidths(dim, 0);
    for (qindex c=0; c<dim; c++) {
        for (qindex r=0; r<dim; r++) {
            int width = form_str(elems[r][c]).length();
            if (width > maxColWidths[c])
                maxColWidths[c] = width;
        }
    }

    // print each row, aligning columns
    for (qindex r=0; r<dim; r++) {
        cout << indent;

        for (qindex c=0; c<dim; c++)
            cout << left
                << setw(maxColWidths[c] + MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS) 
                << setfill(MATRIX_SPACE_CHAR)
                << form_str(elems[r][c]);

        cout << endl;
    }
}

// type T can be qcomp(*)[] or qcomp**
template <typename T> 
void rootPrintDenseMatrixElems(T elems, qindex dim, string indent) {
    if (comm_isRootNode())
        printDenseMatrixElems(elems, dim, indent);
}


void form_printMatrix(CompMatr1 matr, string indent) {
    rootPrintDenseMatrixElems(matr.elems, matr.numRows, indent);
}
void form_printMatrix(CompMatr2 matr, string indent) {
    rootPrintDenseMatrixElems(matr.elems, matr.numRows, indent);
}
void form_printMatrix(CompMatr matr, string indent) {
    rootPrintDenseMatrixElems(matr.cpuElems, matr.numRows, indent);
}



/*
 * DIAGONAL MATRIX PRINTING
 *
 * which print only the matrix diagonals, indented
 */


// diagonals don't need templating because arrays decay to pointers, yay!
void printDiagMatrixElems(qcomp* elems, qindex dim, string indent) {

    // determine max width of each column
    vector<int> maxColWidths(dim, 0);
    for (qindex c=0; c<dim; c++) {
        int width = form_str(elems[c]).length();
        if (width > maxColWidths[c])
            maxColWidths[c] = width;
    }

    // print each "virtual" row, aligning columns
    for (qindex r=0; r<dim; r++) {
        cout << indent;

        // print nothing in every column except diagonal
        for (qindex c=0; c<dim; c++)
            cout << left
                << setw(maxColWidths[c] + MIN_SPACE_BETWEEN_DIAG_MATRIX_COLS) 
                << setfill(MATRIX_SPACE_CHAR)
                << ((r==c)? form_str(elems[c]) : "");

        cout << endl;
    }
}

void rootPrintDiagMatrixElems(qcomp* elems, qindex dim, string indent) {
    if (comm_isRootNode())
        printDiagMatrixElems(elems, dim, indent);
}


void form_printMatrix(DiagMatr1 matr, string indent) {
    rootPrintDiagMatrixElems(matr.elems, matr.numElems, indent);
}
void form_printMatrix(DiagMatr2 matr, string indent) {
    rootPrintDiagMatrixElems(matr.elems, matr.numElems, indent);
}
void form_printMatrix(DiagMatr matr, string indent) {
    rootPrintDiagMatrixElems(matr.cpuElems, matr.numElems, indent);
}

void form_printMatrix(FullStateDiagMatr matr, string indent) {

    // non-distributed edge-case is trivial; root node prints all
    if (!matr.isDistributed) {
        rootPrintDiagMatrixElems(matr.cpuElems, matr.numElemsPerNode, indent);
        return;
    }

    // but when distributed, each node prints its partition
    comm_sync();
    int numNodes = comm_getNumNodes();

    for (int r=0; r<numNodes; r++) {

        if (r == comm_getRank()) {
            std::cout << defaultTableIndent << "[rank " << r << "]" << std::endl;
            printDiagMatrixElems(matr.cpuElems, matr.numElemsPerNode, indent + defaultTableIndent);
        }

        comm_sync();
    }
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
