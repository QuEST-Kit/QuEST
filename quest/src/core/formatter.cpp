/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#include "quest/include/types.h"
#include "quest/include/matrices.h"

#include "quest/src/core/formatter.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

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



/*
 * AESTHETIC CONSTANTS
 */


const int MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS = 2;
const int SET_SPACE_BETWEEN_DIAG_MATRIX_COLS = 2;
const int MIN_SPACE_BETWEEN_TABLE_COLS = 5;

const char TABLE_SPACE_CHAR = '.';
const char MATRIX_SPACE_CHAR = ' ';
const string MATRIX_VDOTS_CHAR = "⋮";   // unicode too wide for char literal
const string MATRIX_DDOTS_CHAR = "⋱";
const string MATRIX_HDOTS_CHAR = "…";



/*
 * USER-CONFIGURABLE GLOBALS
 */


// for diagonal matrices, this is the total number of printed elements,
// but for dense matrices, it is the number in the first column, as well
// as the number of printed columns (starting from the left)
qindex maxNumPrintedMatrixElems = (1 << 5);

void form_setMaxNumPrintedMatrixElems(qindex num) {

    maxNumPrintedMatrixElems = num;
}



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


void form_printMatrixInfo(string nameStr, int numQubits, qindex dim, size_t elemMem, size_t otherMem, int numNodes) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // please don't tell anyone how I live
    bool isDiag = nameStr.find("Diag") != string::npos;

    // prepare substrings
    string sepStr  = ", ";
    string qbStr   = form_str(numQubits) + " qubit" + ((numQubits>1)? "s":"");
    string dimStr  = form_str(dim) + (isDiag? "" : mu + form_str(dim)) + " qcomps";
    string memStr  = form_str(elemMem) + " + " + form_str(otherMem) + by;

    if (numNodes > 1) {
        dimStr += " over " + form_str(numNodes) + " nodes";
        memStr += " per node";
    }

    // print e.g. FullStateDiagMatr (8 qubits, 256 qcomps over 4 nodes, 1024 + 32 bytes per node):
    cout 
        << nameStr << " (" 
        << qbStr   << sepStr 
        << dimStr  << sepStr 
        << memStr  << "):" 
        << endl;
}



/*
 * DENSE MATRIX PRINTING
 *
 * which print only the matrix elements, indented. We truncate matrices which 
 * contain more rows than maxNumPrintedMatrixElems, by printing only their
 * upper-left and lower-left quadrants. This somewhat unusual style of truncation
 * avoids having to align unicode characters and all the resulting headaches.
 */


// type T can be qcomp(*)[] or qcomp**
template <typename T> 
void rootPrintTruncatedDenseMatrixElems(T elems, qindex dim, string indent) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // determine how to truncate the reported matrix
    qindex numTopElems = maxNumPrintedMatrixElems / 2; // floors
    qindex numBotElems = maxNumPrintedMatrixElems - numTopElems; // includes remainder
    qindex numReported = numTopElems + numBotElems;

    // prevent truncation if the matrix is too small, or user has disabled
    bool isTruncated = (
        (maxNumPrintedMatrixElems != 0) &&
        (dim > (numTopElems + numBotElems)));
    if (!isTruncated) {
        numReported = dim;
        numTopElems = dim;
        numBotElems = 0;
    }

    // determine max width of each reported column...
    vector<int> maxColWidths(numReported, 0);
    for (qindex c=0; c<numReported; c++) {

        // by considering the top-most rows
        for (qindex r=0; r<numTopElems; r++) {
            int width = form_str(elems[r][c]).length();
            if (width > maxColWidths[c])
                maxColWidths[c] = width;
        }

        // and the bottom-most rows
        for (qindex i=0; i<numBotElems; i++) {
            int r = dim - i - 1;
            int width = form_str(elems[r][c]).length();
            if (width > maxColWidths[c])
                maxColWidths[c] = width;
        }
    }

    // print top-most rows, aligning columns
    for (qindex r=0; r<numTopElems; r++) {
        cout << indent;

        for (qindex c=0; c<numReported; c++)
            cout << left
                << setw(maxColWidths[c] + MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS) 
                << setfill(MATRIX_SPACE_CHAR)
                << form_str(elems[r][c]);

        // print a trailing ellipsis after the top-most row
        if (isTruncated && r == 0)
            cout << MATRIX_HDOTS_CHAR;

        cout << endl;
    }

    // we are finished if there was no bottom elems (because matrix was not truncated)
    if (!isTruncated)
        return;
    
    // otherwise, separate the top-left and bottom-left quadrants by an ellipsis and blank row,
    // which we align with the center of the left-most column (just to be a little pretty)
    cout << indent;
    for (int i=0; i<maxColWidths[0]/2; i++)
        cout << MATRIX_SPACE_CHAR;
    cout << MATRIX_VDOTS_CHAR << endl;

    // print the remaining bottom rows
    for (qindex r=dim-numBotElems; r<dim; r++) {
        cout << indent;

        for (qindex c=0; c<numReported; c++)
            cout << left
                << setw(maxColWidths[c] + MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS) 
                << setfill(MATRIX_SPACE_CHAR)
                << form_str(elems[r][c]);

        // print a trailing ellipsis after the bottom-most row
        if (r == dim-1)
            cout << MATRIX_HDOTS_CHAR;

        cout << endl;
    }
}

void form_printMatrix(CompMatr1 matr, string indent) {
    rootPrintTruncatedDenseMatrixElems(matr.elems, matr.numRows, indent);
}
void form_printMatrix(CompMatr2 matr, string indent) {
    rootPrintTruncatedDenseMatrixElems(matr.elems, matr.numRows, indent);
}
void form_printMatrix(CompMatr matr, string indent) {
    rootPrintTruncatedDenseMatrixElems(matr.cpuElems, matr.numRows, indent);
}



/*
 * LOCAL DIAGONAL MATRIX PRINTING
 *
 * which print only the matrix diagonals, indented. We don't have to muck around with
 * templating because arrays decay to pointers, yay! We truncate matrices which 
 * contain more diagonals than maxNumPrintedMatrixElems, by printing only their
 * top-most and lowest partitions, separated by ellipsis.
 */


void rootPrintAllDiagonalElems(qcomp* elems, qindex len, qindex indentOffset, string indent) {

    if (!comm_isRootNode())
        return;

    for (qindex i=0; i<len; i++)
        cout 
            << indent
            << string((i + indentOffset) * SET_SPACE_BETWEEN_DIAG_MATRIX_COLS, MATRIX_SPACE_CHAR)
            << form_str(elems[i])
            << endl;
}


void rootPrintDiagonalDots(qindex elemInd, string indent) {

    if (!comm_isRootNode())
        return;

    cout
        << indent 
        << string(elemInd * SET_SPACE_BETWEEN_DIAG_MATRIX_COLS, MATRIX_SPACE_CHAR)
        << MATRIX_DDOTS_CHAR 
        << endl;
}


void rootPrintTruncatedDiagonalMatrixElems(qcomp* elems, qindex dim, string indent) {

    if (!comm_isRootNode())
        return;

    // determine how to truncate the reported matrix
    qindex numTopElems = maxNumPrintedMatrixElems / 2; // floors
    qindex numBotElems = maxNumPrintedMatrixElems - numTopElems; // includes remainder

    // prevent truncation if the matrix is too small, or user has disabled
    bool isTruncated = (
        (maxNumPrintedMatrixElems != 0) &&
        (dim > (numTopElems + numBotElems)));
    if (!isTruncated) {
        numTopElems = dim;
        numBotElems = 0;
    }

    // print each row at a regular indentation, regardless of element width
    rootPrintAllDiagonalElems(elems, numTopElems, 0, indent);

    // we are finished if there was no bottom elems (because matrix was not truncated)
    if (!isTruncated)
        return;

    // otherwise, separate the top-left and bottom-left partitions by a diagonal ellipsis,
    rootPrintDiagonalDots(numTopElems, indent);

    // print remaining elements, keeping consistent indentation, adjusting for above ellipsis
    rootPrintAllDiagonalElems(&elems[dim-numBotElems], numBotElems, numTopElems+1, indent);
}


void form_printMatrix(DiagMatr1 matr, string indent) {
    rootPrintTruncatedDiagonalMatrixElems(matr.elems, matr.numElems, indent);
}
void form_printMatrix(DiagMatr2 matr, string indent) {
    rootPrintTruncatedDiagonalMatrixElems(matr.elems, matr.numElems, indent);
}
void form_printMatrix(DiagMatr matr, string indent) {
    rootPrintTruncatedDiagonalMatrixElems(matr.cpuElems, matr.numElems, indent);
}



/*
 * DISTRIBUTED DIAGONAL MATRIX PRINTING
 *
 * which complicates the truncation printing logic
 */


void rootPrintRank(int rank) {

    if (comm_isRootNode())
        cout << defaultTableIndent << "[rank " << rank << "]" << endl;
}

void rootPrintRanks(int startRank, int endRankIncl) {

    if (comm_isRootNode())
        cout << defaultTableIndent << "[ranks " << startRank << "-" << endRankIncl << "]" << endl;
}


void form_printMatrix(FullStateDiagMatr matr, string indent) {

    // non-distributed edge-case is trivial; root node prints all, handling truncation
    if (!matr.isDistributed) {
        rootPrintTruncatedDiagonalMatrixElems(matr.cpuElems, matr.numElemsPerNode, indent);
        return;
    }

    int thisRank = comm_getRank();
    int numRanks = comm_getNumNodes();

    // prevent truncation if user has disabled by imitating a sufficiently large truncation threshold
    qindex maxPrintedElems = (maxNumPrintedMatrixElems == 0)?
        matr.numElems :
        maxNumPrintedMatrixElems;

    // when distributed, multiple nodes may print their elements depending on the matrix truncation
    int numNodesWorth = maxPrintedElems / matr.numElemsPerNode; // floors

    // these nodes are divided into "full" nodes which will print all their local elements...
    int numFullReportingNodes = 2 * (numNodesWorth/2); // floors
    int numTopFullReportingNodes = numFullReportingNodes / 2; // floors
    int numBotFullReportingNodes = numFullReportingNodes - numTopFullReportingNodes;

    // and two "partial" nodes which print only some of their elements
    int topPartialRank = numTopFullReportingNodes;
    int botPartialRank = numRanks - numBotFullReportingNodes - 1;

    // although it's possible these "partial" nodes end up printing all or none of their elements
    qindex numPartialElems = maxPrintedElems - (numFullReportingNodes * matr.numElemsPerNode);
    qindex numTopPartialElems = numPartialElems / 2; // floors
    qindex numBotPartialElems = numPartialElems - numTopPartialElems;

    // determine which, if any, nodes are "occluded" (do not print at all)
    int topLastPrintingRank  = (numTopPartialElems > 0)? topPartialRank : topPartialRank - 1;
    int botFirstPrintingRank = (numBotPartialElems > 0)? botPartialRank : botPartialRank + 1;
    int firstOccludedRank = topLastPrintingRank + 1;
    int lastOccludedRank  = botFirstPrintingRank - 1;
    bool anyRanksOccluded = lastOccludedRank > firstOccludedRank;

    // nodes print by sending their elements to the root node's buffer, which root then prints.
    // We will allocate the minimum size buffer possible (without incurring more rounds of
    // communication) to avoid astonishing the user with a big allocation. We also only alloc
    // the buffer on the root node, to save memory when each machine runs multiple MPI nodes.
    std::vector<qcomp> rootBuff;

    // if there are any nodes which get all their elements printed...
    if (numTopFullReportingNodes > 0) {

        // the first is always root, which is special; it can print without communication
        rootPrintRank(0);
        rootPrintAllDiagonalElems(matr.cpuElems, matr.numElemsPerNode, 0, indent);

        // if there are more full nodes to print, root allocs its buffer space
        if (numTopFullReportingNodes > 1)
            if (comm_isRootNode(thisRank))
                rootBuff.reserve(matr.numElemsPerNode); 

        // then all top full nodes (except root) send their elems in-turn to root, which prints on their behalf
        for (int r=1; r<numTopFullReportingNodes; r++) {
            rootPrintRank(r);
            comm_sendAmpsToRoot(r, matr.cpuElems, rootBuff.data(), matr.numElemsPerNode);
            rootPrintAllDiagonalElems(rootBuff.data(), matr.numElemsPerNode, 0, indent);
        }
    }

    // if the top "partial node" has any elements to print...
    if (numTopPartialElems > 0) {
        rootPrintRank(topPartialRank);

        // and it happens to be the root node, then printing requires no communication
        if (comm_isRootNode(topPartialRank))
            rootPrintAllDiagonalElems(matr.cpuElems, numTopPartialElems, 0, indent);

        else {
            // otherwise, we ensure root buffer has space (there may have been no top full nodes)
            if (comm_isRootNode(thisRank))
                rootBuff.reserve(numTopPartialElems);

            // send a subset (or all) of the top partial node's amps to root, which prints
            comm_sendAmpsToRoot(topPartialRank, matr.cpuElems, rootBuff.data(), numTopPartialElems);
            rootPrintAllDiagonalElems(rootBuff.data(), numTopPartialElems, 0, indent);
        }

        // if not all elems were printed, add an ellipsis
        if (numTopPartialElems < matr.numElemsPerNode)
            rootPrintDiagonalDots(numTopPartialElems, indent);
    }

    // if any nodes are occluded (they don't have any elems printed), root prints them as ellipsis
    if (anyRanksOccluded) {
        rootPrintRanks(firstOccludedRank, lastOccludedRank);
        rootPrintDiagonalDots(0, indent);
    }

    // if the bottom "partial node" has any elements to print...
    if (numBotPartialElems > 0) {
        rootPrintRank(botPartialRank);

        // ensure root buffer has space (bottom partial can have more elements than top due to truncation indivisibility)
        if (comm_isRootNode(thisRank))
            rootBuff.reserve(numBotPartialElems);

        // if not all elems will be printed, prefix with ellipsis, which will shift subsequently printed elems
        bool isPartial = numBotPartialElems < matr.numElemsPerNode;
        if (isPartial)
            rootPrintDiagonalDots(0, indent);

        // send a subset (or all) of the bottom partial node's amps to root, which prints
        qcomp* elems = &matr.cpuElems[matr.numElemsPerNode - numBotPartialElems];
        comm_sendAmpsToRoot(botPartialRank, elems, rootBuff.data(), numBotPartialElems);
        rootPrintAllDiagonalElems(rootBuff.data(), numBotPartialElems, (int) isPartial, indent);
    }

    // finally, print any bottom "full" nodes
    if (numBotFullReportingNodes > 0) {

        // root pedantically ensures necessary buffer space
        if (comm_isRootNode(thisRank))
            rootBuff.reserve(matr.numElemsPerNode);

        // all bottom full nodes send their elems in-turn to root, which prints on their behalf
        for (int i=0; i<numBotFullReportingNodes; i++) {
            int r = numRanks - numBotFullReportingNodes + i;
            rootPrintRank(r);
            comm_sendAmpsToRoot(r, matr.cpuElems, rootBuff.data(), matr.numElemsPerNode);
            rootPrintAllDiagonalElems(rootBuff.data(), matr.numElemsPerNode, 0, indent);
        }
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
