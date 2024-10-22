/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#include "quest/include/types.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"
#include "quest/include/paulis.h"

#include "quest/src/core/printer.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <type_traits> 
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <tuple>
#include <string>

using namespace printer_substrings;

using std::cout;
using std::endl;
using std::left;
using std::setw;
using std::tuple;
using std::string;
using std::vector;
using std::setfill;
using std::complex;



/*
 * AESTHETIC CONSTANTS
 */


const int MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS = 2;
const int SET_SPACE_BETWEEN_DIAG_MATRIX_COLS = 2;
const int MIN_SPACE_BETWEEN_TABLE_COLS = 5;
const int MIN_SPACE_BETWEEN_COEFF_AND_PAULI_STR = 2;

// must be kept as single char; use numbers above to change spacing
const char TABLE_SPACE_CHAR = '.';
const char MATRIX_SPACE_CHAR = ' ';
const char PAULI_STR_SPACE_CHAR = ' ';

const string VDOTS_CHAR = "⋮";   // unicode too wide for char literal
const string DDOTS_CHAR = "⋱";
const string HDOTS_CHAR = "…";



/*
 * SUBSTRINGS USED BY REPORTERS
 */

namespace printer_substrings { 
    string eq = " = ";
    string mu = " x ";
    string pl = " + ";
    string by = " bytes";
    string pn = " per node";
    string pg = " per gpu";
    string pm = " per machine";
    string bt = "2^";
    string na = "N/A";
    string un = "unknown";
}



/*
 * USER-CONFIGURABLE GLOBALS
 */

// controls the truncation of reported matrices and pauli string sums.
// for diagonal matrices, this is the total number of printed elements,
// but for dense matrices, it is the number in the first column, as well
// as the number of printed columns (starting from the left). this is
// user-configurable and can be set to 0 to disable truncation
qindex maxNumPrintedItems = (1 << 5);
bool isTruncationEnabled  = (maxNumPrintedItems != 0);

void printer_setMaxNumPrintedItems(qindex num) {
    maxNumPrintedItems  = num;
    isTruncationEnabled = (num != 0);
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

// otherwise just give up; reporters might print mangled types
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


string printer_getQcompType() {

    qcomp x;
    return getTypeName(x);
}

string printer_getQrealType() {

    // more portable than getTypeName
    return GET_STR( FLOAT_TYPE );
}

string printer_getQindexType() {

    // more portable than getTypeName
    return GET_STR( INDEX_TYPE );
}

string printer_getFloatPrecisionFlag() {

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
    std::stringstream buffer;

    // ensure +- is not shown if forced hidden
    if (hideSign) {
        buffer << std::noshowpos;
        if (num < 0)
            num *= -1;
    }

    // return 0, 1, 0.1, 1.2e-5, -1e+5
    buffer << num;
    return buffer.str();
}


// type T can be precision decimal (independent of qreal) i.e. float, double, long double
template <typename T>
string printer_toStr(complex<T> num) {

    // precise 0 is rendered as a real integer
    if (real(num) == 0 && imag(num) == 0)
        return "0";

    // real complex-floats neglect imag component
    if (imag(num) == 0)
        return floatToStr(real(num));

    // imaginary complex-floats neglect real component entirely (instead of 0+3i)
    string realStr = (real(num) == 0)? "" : floatToStr(real(num));

    // scientific-notation real component (when followed by any imaginary component) gets brackets
    if (realStr.find('e') != string::npos)
        realStr = '(' + realStr + ')';

    // -1i is abbreviated to -i
    if constexpr (std::is_signed_v<T>)
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
template string printer_toStr<int>(complex<int> num);
template string printer_toStr<long int>(complex<long int> num);
template string printer_toStr<long long int>(complex<long long int> num);

template string printer_toStr<unsigned>(complex<unsigned> num);
template string printer_toStr<long unsigned>(complex<long unsigned> num);
template string printer_toStr<long long unsigned>(complex<long long unsigned> num);

template string printer_toStr<float>(complex<float> num);
template string printer_toStr<double>(complex<double> num);
template string printer_toStr<long double>(complex<long double> num);


// alias as toStr() just for internal brevity
// (this seems backward; ordinarily we would define toStr() as
// the templated inner-function and define concretely-typed public
// printer_toStr() overloads. This current implementation, where
// we expose the templated printer_toStr() in the header, avoids
// polluting the header with the above precision-agnostic 
// instantiations (which would become overloads). Stinky!)

template <class T> string toStr(T num) {

    return printer_toStr(num);
}



/*
 * PRIMITIVE PRINTING
 *
 * used by other files which need to directly print
 * unformatted strings (from root node only)
 */


void print(const char* str) {
    if (comm_isRootNode())
        cout << str << endl;
}

void print(string str) {
    print(str.data());
}

void print(qcomp num) {
    print(toStr(num));
}



/*
 * INTERNAL FUNCTIONS
 *
 * used by the below public functions for preparing substrings 
 */


string getNumQubitsStr(int numQubits) {

    // e.g. "1 qubit" or "2 qubits"
    return toStr(numQubits) + " qubit" + ((numQubits>1)? "s":"");
}


string getNumMatrixElemsStr(qindex dim, bool isDiag, int numNodes) {

    // we do not bother handling the dim=1 non-plural case, because it
    // never occurs! We always receive power-of-2 dimension matrices

    // e.g. "2 qcomps" or "2x2 qcomps"
    string out = toStr(dim) + (isDiag? "" : mu + toStr(dim)) + " qcomps";

    // FullStateDiagMatr may be distributed over >1 numNodes (for other matrix types, numNodes=1 always)
    if (numNodes > 1)
        out += " over " + toStr(numNodes) + " nodes";

    return out;
}


string getMemoryCostsStr(size_t numMainBytes, size_t numOverheadBytes, int numNodes) {

    // e.g. "1024 + 32 bytes"
    string out = toStr(numMainBytes) + pl + toStr(numOverheadBytes) + by;

    // FullStateDiagMatr may be distributed over >1 numNodes (for other matrix types, numNodes=1 always)
    if (numNodes > 1)
        out += " per node";

    return out;
}



/*
 * MATRIX INFO PRINTING
 *
 * which prints information about the matrix type itself, rather than its elements
 */


void print_matrixInfo(string nameStr, int numQubits, qindex dim, size_t elemMem, size_t otherMem, int numNodes, bool hasGpuMem) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // please don't tell anyone how I live
    // edit: Tyson from 2 months in the future has found this and is appalled, but is somewhat
    // sated by the potential bug (a matrix being incorrectly assessed as diagonal if its name
    // changes) causing only a minor printed/typeset error (e.g. showing "5x5" instead of "5") 
    bool isDiag = nameStr.find("Diag") != string::npos;

    // print e.g. FullStateDiagMatr (8 qubits, 256 qcomps over 4 nodes, 1024 + 32 bytes per node):
    cout 
        << nameStr 
        << " (" 
        << getNumQubitsStr(numQubits)
        << ", " 
        << getNumMatrixElemsStr(dim, isDiag, numNodes)
        << ", "  
        << getMemoryCostsStr(elemMem, otherMem, numNodes)
        << (hasGpuMem? (" in GPU") : "") 
        << "):"  // colon, since caller will hereafter print matrix elements
        << endl;
}



/*
 * DENSE MATRIX PRINTING
 *
 * which print only the matrix elements, indented. We truncate matrices which 
 * contain more rows than maxNumPrintedItems, by printing only their
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
    qindex numTopElems = maxNumPrintedItems / 2; // floors
    qindex numBotElems = maxNumPrintedItems - numTopElems; // includes remainder
    qindex numReported = numTopElems + numBotElems;

    // prevent truncation if the matrix is too small, or user has disabled
    bool isTruncated = isTruncationEnabled && dim > (numTopElems + numBotElems);
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
            int width = toStr(elems[r][c]).length();
            if (width > maxColWidths[c])
                maxColWidths[c] = width;
        }

        // and the bottom-most rows
        for (qindex i=0; i<numBotElems; i++) {
            int r = dim - i - 1;
            int width = toStr(elems[r][c]).length();
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
                << toStr(elems[r][c]);

        // print a trailing ellipsis after the top-most row
        if (isTruncated && r == 0)
            cout << HDOTS_CHAR;

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
    cout << VDOTS_CHAR << endl;

    // print the remaining bottom rows
    for (qindex r=dim-numBotElems; r<dim; r++) {
        cout << indent;

        for (qindex c=0; c<numReported; c++)
            cout << left
                << setw(maxColWidths[c] + MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS) 
                << setfill(MATRIX_SPACE_CHAR)
                << toStr(elems[r][c]);

        // print a trailing ellipsis after the bottom-most row
        if (r == dim-1)
            cout << HDOTS_CHAR;

        cout << endl;
    }
}

void print_matrix(CompMatr1 matr, string indent) {
    rootPrintTruncatedDenseMatrixElems(matr.elems, matr.numRows, indent);
}
void print_matrix(CompMatr2 matr, string indent) {
    rootPrintTruncatedDenseMatrixElems(matr.elems, matr.numRows, indent);
}
void print_matrix(CompMatr matr, string indent) {
    rootPrintTruncatedDenseMatrixElems(matr.cpuElems, matr.numRows, indent);
}



/*
 * LOCAL DIAGONAL MATRIX PRINTING
 *
 * which print only the matrix diagonals, indented. We don't have to muck around with
 * templating because arrays decay to pointers, yay! We truncate matrices which 
 * contain more diagonals than maxNumPrintedItems, by printing only their
 * top-most and lowest partitions, separated by ellipsis.
 */


void rootPrintAllDiagonalElems(qcomp* elems, qindex len, qindex indentOffset, string indent) {

    if (!comm_isRootNode())
        return;

    for (qindex i=0; i<len; i++)
        cout 
            << indent
            << string((i + indentOffset) * SET_SPACE_BETWEEN_DIAG_MATRIX_COLS, MATRIX_SPACE_CHAR)
            << toStr(elems[i])
            << endl;
}


void rootPrintDiagonalDots(qindex elemInd, string indent) {

    if (!comm_isRootNode())
        return;

    cout
        << indent 
        << string(elemInd * SET_SPACE_BETWEEN_DIAG_MATRIX_COLS, MATRIX_SPACE_CHAR)
        << DDOTS_CHAR 
        << endl;
}


void rootPrintTruncatedDiagonalMatrixElems(qcomp* elems, qindex dim, string indent) {

    if (!comm_isRootNode())
        return;

    // determine how to truncate the reported matrix
    qindex numTopElems = maxNumPrintedItems / 2; // floors
    qindex numBotElems = maxNumPrintedItems - numTopElems; // includes remainder

    // prevent truncation if the matrix is too small, or user has disabled
    bool isTruncated = isTruncationEnabled && dim > (numTopElems + numBotElems);
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


void print_matrix(DiagMatr1 matr, string indent) {
    rootPrintTruncatedDiagonalMatrixElems(matr.elems, matr.numElems, indent);
}
void print_matrix(DiagMatr2 matr, string indent) {
    rootPrintTruncatedDiagonalMatrixElems(matr.elems, matr.numElems, indent);
}
void print_matrix(DiagMatr matr, string indent) {
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


void print_matrix(FullStateDiagMatr matr, string indent) {

    // non-distributed edge-case is trivial; root node prints all, handling truncation
    if (!matr.isDistributed) {
        rootPrintTruncatedDiagonalMatrixElems(matr.cpuElems, matr.numElemsPerNode, indent);
        return;
    }

    int thisRank = comm_getRank();
    int numRanks = comm_getNumNodes();

    // prevent truncation if user has disabled by imitating a sufficiently large truncation threshold
    qindex maxPrintedElems = matr.numElems;
    if (isTruncationEnabled && maxNumPrintedItems < maxPrintedElems)
        maxPrintedElems = maxNumPrintedItems;

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
 * CHANNEL PRINTING
 */


void print_superOpInfo(SuperOp op) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // determine memory costs
    size_t elemMem = mem_getLocalSuperOpMemoryRequired(op.numQubits);
    size_t structMem = sizeof(op);

    // print e.g. SuperOp (1 qubit, 4 x 4 qcomps, 256 + 40 bytes)
    cout 
        << "SuperOp" 
        << " (" 
        << getNumQubitsStr(op.numQubits)
        << ", " 
        << getNumMatrixElemsStr(op.numRows, false, 1)  // superoperator is a non-diagonal non-distributed matrix
        << ", "  
        << getMemoryCostsStr(elemMem, structMem, 1) 
        << "):"  // colon, since caller will hereafter print superoperator elements
        << endl;
}

void print_superOp(SuperOp op, string indent) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    rootPrintTruncatedDenseMatrixElems(op.cpuElems, op.numRows, indent);
}

void print_krausMapInfo(KrausMap map) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // determine memory costs (gauranteed not to overflow)
    size_t krausMem = mem_getLocalMatrixMemoryRequired(map.numQubits, true, 1) * map.numMatrices;
    size_t superMem = mem_getLocalSuperOpMemoryRequired(map.superop.numQubits);
    size_t strucMem = sizeof(map);




    // print e.g. 
    cout 
        << "KrausMap" 
        << " (" 
        << getNumQubitsStr(map.numQubits)
        << ", " 
        << toStr(map.numMatrices) + " " + ((map.numMatrices>1)? "matrices" : "matrix")
        << ", " 
        << getNumMatrixElemsStr(map.numRows, false, 1) 
        << ", "  
        << toStr(krausMem) << " + " << getMemoryCostsStr(superMem, strucMem, 1) 
        << "):"  // colon, since caller will thereafter print Kraus matrices (via print_krausMap())
        << endl;
}

void print_krausMap(KrausMap map, string indexIndent, string matrixIndent) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // print every Kraus operator in-turn (regardless of reporter truncation settings)
    for (int n=0; n<map.numMatrices; n++) {

        // print the Kraus map index
        cout << indexIndent << "[" << toStr(n) << "]" << endl;

        // then its corresponding matrix (which may be truncated)
        rootPrintTruncatedDenseMatrixElems(map.matrices[n], map.numRows, matrixIndent);
    }
}



/*
 * TABLE PRINTING
 */


void print_table(string title, vector<tuple<string, string>> rows, string indent) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // find max-width of left column
    size_t maxWidth = 0;
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


void print_table(string title, vector<tuple<string, long long int>> rows, string indent) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // convert all values to strings
    vector<tuple<string, string>> casted;
    for (auto const& [key, value] : rows)
        casted.push_back({key, std::to_string(value)});

    print_table(title, casted, indent);
}



/*
 * PAULI STRING AND SUM REPORTERS
 */


// we'll make use of these internal functions from paulis.cpp
extern int paulis_getPauliAt(PauliStr str, int ind);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str);


string getPauliStrAsString(PauliStr str, int maxNumQubits) {

    string labels = "IXYZ";

    // ugly but adequate
    string out = "";
    for (int i=maxNumQubits-1; i>=0; i--)
        out += labels[paulis_getPauliAt(str, i)];
    
    return out;
}


void print_pauliStr(PauliStr str, int numQubits) {

    // only root node ever prints
    if (!comm_isRootNode())
        return;

    cout << getPauliStrAsString(str, numQubits) << endl;
}


void print_pauliStrSumInfo(qindex numTerms, qindex numBytes) {

    // only root node ever prints
    if (!comm_isRootNode())
        return;

    // print e.g. PauliStrSum (3 terms, 24 bytes)
    cout 
        << "PauliStrSum (" 
        << numTerms << " term" << ((numTerms>1)? "s":"") << ", "
        << numBytes << by 
        << "):"  // colon, since caller will hereafter print each pauli string
        << endl;
}


int getWidthOfWidestElem(qcomp* elems, qindex numElems) {

    int maxWidth = 0;

    for (qindex i=0; i<numElems; i++) {
        int width = toStr(elems[i]).size();
        if (width > maxWidth)
            maxWidth = width;
    }

    return maxWidth;
}


int getIndexOfLeftmostPauliAmongStrings(PauliStr* strings, qindex numStrings) {

    int maxInd = 0;

    for (qindex i=0; i<numStrings; i++) {
        int ind = paulis_getIndOfLefmostNonIdentityPauli(strings[i]);
        if (ind > maxInd)
            maxInd = ind;
    }

    return maxInd;
}


void printPauliStrSumSubset(PauliStrSum sum, qindex startInd, qindex endIndExcl, int coeffWidth, int maxNumPaulis, string indent) {

    // each line prints as e.g. -3+2i XIXXI
    for (int i=startInd; i<endIndExcl; i++)
        cout 
            << indent << left 
            << setw(coeffWidth + MIN_SPACE_BETWEEN_COEFF_AND_PAULI_STR) 
            << setfill(PAULI_STR_SPACE_CHAR) 
            << toStr(sum.coeffs[i]) 
            << getPauliStrAsString(sum.strings[i], maxNumPaulis)
            << endl;
}


void print_pauliStrSum(PauliStrSum sum, string indent) {
    
    // only root node ever prints
    if (!comm_isRootNode())
        return;

    // determine how to truncate the reported pauli string
    qindex numTopElems = maxNumPrintedItems / 2; // floors
    qindex numBotElems = maxNumPrintedItems - numTopElems; // includes remainder
    qindex botElemsInd = sum.numTerms - numBotElems;

    // prevent truncation if the matrix is too small, or user has disabled
    bool isTruncated = isTruncationEnabled && sum.numTerms > maxNumPrintedItems;
    if (!isTruncated) {
        numTopElems = sum.numTerms;
        numBotElems = 0;
        // botElemsInd is disregarded
    }

    // choose coeffs padding to make consistent column alignment
    int maxTopCoeffWidth = getWidthOfWidestElem(sum.coeffs, numTopElems);
    int maxBotCoeffWidth = getWidthOfWidestElem(&(sum.coeffs[botElemsInd]), numBotElems);
    int maxCoeffWidth = std::max(maxTopCoeffWidth, maxBotCoeffWidth);

    // choose the fixed number of Paulis to show (the min necessary to fit all non-truncated strings)
    int maxTopNumPaulis = 1 + getIndexOfLeftmostPauliAmongStrings(sum.strings, numTopElems);
    int maxBotNumPaulis = 1 + getIndexOfLeftmostPauliAmongStrings(&(sum.strings[botElemsInd]), numBotElems);
    int maxNumPaulis = std::max(maxTopNumPaulis,maxBotNumPaulis);

    // print all top elems
    printPauliStrSumSubset(sum, 0, numTopElems, maxCoeffWidth, maxNumPaulis, indent);

    // if we just printed all elems, there's nothing left to do
    if (!isTruncated)
        return;

    // otherwise print ellipsis between coeff and string columns
    cout << indent << string(maxCoeffWidth+1, PAULI_STR_SPACE_CHAR) << VDOTS_CHAR << endl;

    // print all bottom elems
    printPauliStrSumSubset(sum, botElemsInd, sum.numTerms, maxCoeffWidth, maxNumPaulis, indent);
}
