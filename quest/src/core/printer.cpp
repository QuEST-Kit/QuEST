/** @file
 * Functions for formatting and outputting strings, used
 * primarily by reporters (e.g. reportQureg). A substantial
 * amount of logic is dedicated to making outputs look
 * extra pretty, managing distributed objects, and 
 * aesthetically truncating large matrices.
 * 
 * @author Tyson Jones
 * @author Erich Essmann (improved OS agnosticism, patched mem-leak)
 */

#include "quest/include/qureg.h"
#include "quest/include/types.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"
#include "quest/include/paulis.h"

#include "quest/src/core/printer.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <type_traits> 
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include <new>


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
 * PRIVATE AESTHETIC CONSTANTS
 */


const int MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS = 2;
const int SET_SPACE_BETWEEN_DIAG_MATRIX_COLS = 2;
const int MIN_SPACE_BETWEEN_TABLE_COLS = 5;
const int SET_SPACE_BETWEEN_KET_AND_RANK = 2;

// must be kept as single char; use numbers above to change spacing
const char TABLE_SPACE_CHAR = '.';
const char MATRIX_SPACE_CHAR = ' ';
const char PAULI_SPACE_CHAR = ' ';

// should not contain newline
const string LABEL_SUFFIX = ":";

const string VDOTS_CHAR = "⋮";   // unicode too wide for char literal
const string DDOTS_CHAR = "⋱";
const string HDOTS_CHAR = "…";

const string HEADER_OPEN_BRACKET = " (";
const string HEADER_CLOSE_BRACKET = "):";
const string HEADER_DELIMITER = ", ";

const string IMAGINARY_UNIT = "i";



/*
 * PUBLIC AESTHETIC CONSTATNS
 */


namespace printer_substrings { 
    string eq = " = ";
    string pn = " per node";
    string pg = " per gpu";
    string ig = " in gpu";
    string pm = " per machine";
    string bt = "2^";
    string na = "N/A";
    string un = "unknown";
}



/*
 * USER-CONFIGURABLE GLOBALS
 */


// user truncation of printed quregs, matrices, pauli-strings, etc
qindex global_maxNumPrintedRows = (1 << 5);
qindex global_maxNumPrintedCols = 4;

void printer_setMaxNumPrintedScalars(qindex numRows, qindex numCols) {

    global_maxNumPrintedRows = numRows;
    global_maxNumPrintedCols = numCols;
}


// user configuration of printed scalars
int global_maxNumPrintedSigFigs = 5;

void printer_setMaxNumPrintedSigFig(int numSigFigs) {

    global_maxNumPrintedSigFigs = numSigFigs;
}


// user configuration of trailing newlines after reporters

int global_numTrailingNewlines = 2;

void printer_setNumTrailingNewlines(int numNewlines) {

    global_numTrailingNewlines = numNewlines;
}

int printer_getNumTrailingNewlines() {

    return global_numTrailingNewlines;
}


// user specified characters to represent Pauli operators

string global_pauliChars = "IXYZ";

void printer_setPauliChars(string newChars) {

    global_pauliChars = newChars;
}


// user specified Pauli string format

int global_pauliStrFormatFlag = 0;

void printer_setPauliStrFormat(int flag) {

    /// @todo
    /// this could be changed to an enum, even though the
    /// API will likely always receive an int (for C and
    /// C++ agnosticism and simplicity - I hate polluting
    /// an API with custom types and constants!).

    global_pauliStrFormatFlag = flag;
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

    // #warning command safe in non-MSVC compiler
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


inline std::string demangleTypeName(const char* mangledName) {
#if defined(__GNUC__) || defined(__clang__)
    int status = 0;

    // __cxa_demangle returns a malloc'd string, so we wrap it in a unique_ptr
    // for automatic cleanup. (The custom deleter calls free().)
    std::unique_ptr<char, void(*)(void*)> demangled(
        abi::__cxa_demangle(mangledName, nullptr, nullptr, &status),
        std::free
    );

    if (status == 0 && demangled) {
        return demangled.get();
    } else {
        // fall back to the mangled name
        return mangledName;
    }
#else
    // e.g. MSVC or unknown compiler: no standard demangler
    return mangledName;
#endif
}


// type T can be anything in principle, although it's currently only used for qcomp
template <typename T>
std::string getTypeName(T _unused) {
    // For MSVC, typeid(T).name() typically returns something like "class Foo"
    // or "struct Foo", but it's still not exactly "Foo".
    // For GCC/Clang, you get a raw "mangled" name, e.g. "N3FooE".
    const char* rawName = typeid(T).name();

    // We'll try to demangle if we can:
    return demangleTypeName(rawName);
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
string floatToStr(T num, bool hideSign=false, int overrideSigFigs=-1) {

    // write to stream (instead of calling to_string()) to auto-use scientific notation
    std::stringstream buffer;

    // ensure +- is not shown if forced hidden
    if (hideSign) {
        buffer << std::noshowpos;
        if (num < 0)
            num *= -1;
    }

    // impose user-set significant figures, unless caller overrode
    int sigFigs = (overrideSigFigs != -1)? overrideSigFigs : global_maxNumPrintedSigFigs;
    buffer << std::setprecision(sigFigs);

    // get string of e.g. 0, 1, 0.1, 1.2e-05, -1e+05
    buffer << num;
    string out = buffer.str();

    // remove superflous 0 prefix in scientific notation exponent
    size_t ePos = out.find('e');
    if (ePos != string::npos && (out[ePos + 2] == '0'))
        out.erase(ePos + 2, 1);

    // permit superflous + prefix of scientific notatio exponent, e.g. 2e+6

    return out;
}


// type T can be precision decimal (independent of qreal) i.e. float, double, long double
template <typename T>
string printer_toStr(complex<T> num) {

    // precise 0 is rendered as a real integer
    if (std::real(num) == 0 && std::imag(num) == 0)
        return "0";

    // real complex-floats neglect imag component
    if (std::imag(num) == 0)
        return floatToStr(std::real(num));

    // imaginary complex-floats neglect real component entirely (instead of 0+3i)
    string realStr = (std::real(num) == 0)? "" : floatToStr(std::real(num));

    // scientific-notation real component (when followed by any imaginary component) gets brackets
    if (realStr.find('e') != string::npos)
        realStr = '(' + realStr + ')';

    // -1i is abbreviated to -i
    if constexpr (std::is_signed_v<T>)
        if (std::imag(num) == -1)
            return realStr + "-" + IMAGINARY_UNIT;

    // +1i is abbreviated to +i, although the sign is dropped if there's no preceeding real component
    if (std::imag(num) == 1)
        return realStr + ((std::real(num) == 0)? IMAGINARY_UNIT : ("+" + IMAGINARY_UNIT));

    // get imag component string but without sign
    string imagStr = floatToStr(std::imag(num), true);

    // scientific-notation components always get wrapped in paranthesis
    if (imagStr.find('e') != string::npos)
        imagStr = '(' + imagStr + ')';

    // negative imag components always get - sign prefix
    if (std::imag(num) < 0)
        imagStr = '-' + imagStr;

    // positive imag components only get + sign prefix if they follow a real component
    if (std::imag(num) > 0 && std::real(num) != 0)
        imagStr = '+' + imagStr;

    // all imag components end in 'i'
    imagStr += IMAGINARY_UNIT;

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


// explicit qreal overload so that real sig-figs can be changed
string printer_toStr(qreal num) {

    // uses user-set significant figures
    return floatToStr(num);
}


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

    if (!comm_isRootNode())
        return;

    // never prints a trailing newline
    cout << str;
}

void print(string str) {
    print(str.data());
}

void print(qcomp num) {
    print(toStr(num));
}



/*
 * NEWLINE PRINTING
 */


void print_newlines() {
    assert_printerGivenNonNegativeNumNewlines();

    print(string(global_numTrailingNewlines, '\n'));
}

void print_oneFewerNewlines() {
    assert_printerGivenPositiveNumNewlines();

    print(string(global_numTrailingNewlines - 1, '\n'));
}



/*
 * STRINGIFYING FUNCTIONS
 *
 * some of which are exposed to other files, while others
 * are private convenience functions to prepare substrings
 */


string printer_getMemoryWithUnitStr(size_t numBytes) {

    // retain "bytes" instead of "B" since more recognisable
    vector<string> units = {"bytes", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"};

    // unit sizes (in bytes) may overflow; but then so too will numBytes
    vector<size_t> sizes(units.size());
    for (size_t i=0; i<units.size(); i++)
        sizes[i] = powerOf2(10 * i); // = 1024^i

    // find the biggest unit for which numBytes/sizes[ind] > 1
    int ind = 1;
    while (numBytes > sizes[ind])
        ind++;
    ind--;

    // express numBytes in terms of new unit, forcefully rounding to 3 sig-figs max,
    // except when the chosen unit is bytes (then we permit all 4 digits)
    qreal frac = numBytes / static_cast<qreal>(sizes[ind]);
    return floatToStr(frac, false, (ind==0)? 4 : 3) + " " + units[ind];
}


string getTableTitleStr(string title) {

    return "[" + title + "]";
}


string getRankStr(int rank) {

    return "[node " + toStr(rank) + "]";
}


string getKrausOperatorIndStr(qindex ind) {

    return "[matrix " + toStr(ind) + "]";
}


vector<string> getBasisKets(qindex startInd, qindex numInds) {

    vector<string> kets(numInds);

    for (qindex i=0; i<numInds; i++)
        kets[i] = std::string("|") + toStr(i+startInd) + std::string("⟩");

    return kets;
}


string getNumQubitsStr(int numQubits) {

    // e.g. "1 qubit" or "2 qubits"
    return toStr(numQubits) + " qubit" + (numQubits>1? "s":"");
}


string getQuregTypeStr(Qureg qureg) {

    // never pluralise qubit here because of following word
    return toStr(qureg.numQubits) + " qubit " + (qureg.isDensityMatrix? "density matrix" : "statevector");
}


string getProductStr(qindex a, qindex b) {

    // e.g. "4x3"
    return toStr(a) + "x" + toStr(b);
}


string getNumMatrixElemsStr(qindex dim, bool isDiag, int numNodes) {

    // we do not bother handling the dim=1 non-plural case, because it
    // never occurs! We always receive power-of-2 dimension matrices

    // e.g. "2 qcomps" or "2x2 qcomps"
    string out = (isDiag? toStr(dim) : getProductStr(dim, dim)) + " qcomps"; // never non-plural

    // FullStateDiagMatr may be distributed over >1 numNodes (for other matrix types, numNodes=1 always)
    if (numNodes > 1)
        out += " over " + toStr(numNodes) + " nodes";

    return out;
}


string getMemoryCostsStr(size_t numBytesPerNode, bool isDistrib, bool isGpu) {

    using namespace printer_substrings;

    // when data is GPU-accelerated, we report all memory as being
    // in the GPU; this isn't exactly true, since the GPU VRAM will
    // not contain struct fields nor KrausMap matrices. However,
    // the difference is tiny and unimportant, and shrinks exponentially
    // (or quadratically, for KrausMap) with increasing number of qubits

    string mem = printer_getMemoryWithUnitStr(numBytesPerNode);

    if (isDistrib && isGpu)
        return mem + pg;
    if (isDistrib)
        return mem + pn;
    if (isGpu)
        return mem + ig;
    
    return mem;
}



/*
 * INTERNAL CONVENIENCE TYPES
 *
 * used for compactly passing arguments to internal
 * functions which format and print 1D and 2D matrices
 */


using qcompmatr = vector<vector<qcomp>>;
using stringmatr = vector<vector<string>>;


stringmatr getMatrixOfQcompStrings(qcompmatr in) {

    if (in.empty())
        return stringmatr(0);

    size_t numRows = in.size();
    size_t numCols = in[0].size();
    stringmatr out(numRows, vector<string>(numCols));

    for (size_t r=0; r<numRows; r++)
        for (size_t c=0; c<numCols; c++)
            out[r][c] = toStr(in[r][c]);

    return out;
}



/*
 * PRINTING MATRICES DIVIDED INTO QUADRANTS
 *
 * where the quadrants are those resulting
 * from horizontally and vertically truncating
 * the matrices. The functions in this section
 * are compatible with both 1D and 2D data
 * (vectors become single-column matrices)
 */


size_t getMaxWidthOfColumn(stringmatr matr, size_t col) {

    size_t maxWidth = 0;

    for (size_t r=0; r<matr.size(); r++) {
        size_t width = matr[r][col].size();
        if (width > maxWidth)
            maxWidth = width;
    }

    return maxWidth;
}


vector<size_t> getMaxWidthOfColumns(stringmatr upper, stringmatr lower) {

    if (upper.empty())
        return {};

    size_t numCols = upper[0].size(); // matches lower, unless lower empty
    vector<size_t> maxWidths(numCols, 0);

    for (size_t c=0; c<numCols; c++)
        maxWidths[c] = std::max(
            getMaxWidthOfColumn(upper, c), 
            getMaxWidthOfColumn(lower, c)); // 0 if lower empty

    return maxWidths;
}


void expandMaxWidthsAccordingToLabels(vector<size_t> &widths, vector<string> &labels) {

    if (labels.empty())
        return;

    for (size_t c=0; c<widths.size(); c++) {
        size_t labelWidth = labels[c].size();
        if (labelWidth > widths[c])
            widths[c] = labelWidth;
    }
}


void printPerRowIndent(string baseIndent, size_t row, string indentPerRow) {

    cout << baseIndent;

    // K. I. S. S.
    for (size_t r=0; r<row; r++)
        cout << indentPerRow;
}


void printRowInTwoQuadrants(
    vector<string> leftRow, vector<size_t> leftColWidths,
    vector<string> rightRow, vector<size_t> rightColWidths
) {
    // print left portion of row
    for (size_t c=0; c<leftRow.size(); c++)
        cout << left << setw(leftColWidths[c] + MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS) << setfill(MATRIX_SPACE_CHAR)
             << leftRow[c];

    // finish without printing trailing newline if there is no right portion
    if (rightRow.empty())
        return;

    // draw ellipsis; left-spacer already applied by previous amp
    string spacer = string(MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS, MATRIX_SPACE_CHAR);
    cout << HDOTS_CHAR << spacer;

    // print right portion
    for (size_t c=0; c<rightRow.size(); c++)
        cout << left << setw(rightColWidths[c] + MIN_SPACE_BETWEEN_DENSE_MATRIX_COLS) << setfill(MATRIX_SPACE_CHAR)
             << rightRow[c];

    // do not print a trailing new line; caller may wish to append more characters
}


void printContiguousRowsInTwoQuadrants(
    stringmatr left, vector<size_t> leftColWidths, 
    stringmatr right, vector<size_t> rightColWidths, 
    vector<string> rowLabels,
    string indent, string indentPerRow, size_t rowOffsetForIndent
) {
    for (size_t r=0; r<left.size(); r++) {
        printPerRowIndent(indent, r + rowOffsetForIndent, indentPerRow);

        auto rightRow = (right.empty())? vector<string>{} : right[r];
        printRowInTwoQuadrants(left[r], leftColWidths, rightRow, rightColWidths);

        // print optional row label, using trailing space left by above row print
        cout << (rowLabels.empty()? "" : rowLabels[r]);
        cout << endl;
    }
}


void printMatrixInFourQuadrants(
    qcompmatr upperLeft, qcompmatr upperRight, 
    qcompmatr lowerLeft, qcompmatr lowerRight, 
    vector<string> leftColLabels,  vector<string> rightColLabels, // shown above the first row
    vector<string> upperRowLabels, vector<string> lowerRowLabels, // shown right of the last column
    string baseIndent, string indentPerRow, string verticalEllipsis
) {
    // this function prints a matrix which has been prior divided into one, two or four
    // quadrants, resulting from reporter truncation of the original matrix. It can be
    // used to print truncated (or non-truncated) non-square matrices, horizontal and
    // vertical vectors, diagonal matrices, and lists of scalars (e.g. Pauli coefficients)

    // the below preconditions are assumed but not enforced:
    // - it is valid for lowerLeft to be empty, implying upperLeft contains entire columns;
    //   in that scenario, lowerRight must also be empty
    // - it is valid for upperRight to be empty, implying upperLeft contains entire rows;
    //   in that scenario, lowerRight must also be empty
    // - it is valid for upperRight=lowerLeft=lowerRight to be empty, implying upperLeft
    //   contains the entire matrix
    // - it is valid for qcompmatr={} (empty), but it is INVALID to be {{}} (reported non-empty)
    // - it is INVALID for any qcompmatr to have an inconsistent number of columns across its rows
    // - it is INVALID for upperLeft and lowerLeft to differ in #columns, though #rows can vary;
    //   the same applies to upperRight and lowerRight
    // - it is INVALID for upperLeft and upperRight to differ in #rows, though #columns can vary;
    //   the same applies to lowerLeft and lowerRight
    // - it is INVALID for upperLeft to ever be empty
    // - it is valid for leftColLabels=rightColLabels to be empty, but otherwise...
    //   - leftColLabels must match the #cols in upperLeft (can contain empty strings)
    //   - rightColLabels must similarly match the #cols in upperRight
    // - it is valid for upperRight=lowerRight to be empty, while lowerLeft is not empty;
    //   this would be the case when printing a vector. Non-distributed vectors like this
    //   are fine, although distributed vectors should use a separate, bespoke routine since
    //   this function cannot display labeled ranks spanning vertically

    if (!comm_isRootNode())
        return;

    // prepare element strings
    stringmatr upperLeftStrings  = getMatrixOfQcompStrings(upperLeft);
    stringmatr upperRightStrings = getMatrixOfQcompStrings(upperRight);
    stringmatr lowerLeftStrings  = getMatrixOfQcompStrings(lowerLeft);
    stringmatr lowerRightStrings = getMatrixOfQcompStrings(lowerRight);

    // decide column widths
    vector<size_t> leftColWidths  = getMaxWidthOfColumns(upperLeftStrings,  lowerLeftStrings);
    vector<size_t> rightColWidths = getMaxWidthOfColumns(upperRightStrings, lowerRightStrings);

    // including consideration of column labels
    expandMaxWidthsAccordingToLabels(leftColWidths,  leftColLabels);
    expandMaxWidthsAccordingToLabels(rightColWidths, rightColLabels);

    // optionally print all column labels
    if (!leftColLabels.empty()) {
        cout << baseIndent;
        printRowInTwoQuadrants(leftColLabels, leftColWidths, rightColLabels, rightColWidths);
        cout << endl;
    }

    // print the upper rows
    printContiguousRowsInTwoQuadrants(
        upperLeftStrings, leftColWidths, 
        upperRightStrings, rightColWidths, 
        upperRowLabels,
        baseIndent, indentPerRow, 0);

    // finish immediately if there are no lower rows (matrix wasn't vertically truncated)
    if (lowerLeftStrings.empty())
        return;

    // otherwise, draw ellipsis, aligned (approx) in the midpoint of leftmost column (which may be indented)
    printPerRowIndent(baseIndent, upperLeft.size(), indentPerRow);
    cout << string(leftColWidths[0]/2, MATRIX_SPACE_CHAR) << verticalEllipsis << endl;

    // print the lower rows (continuing indentation from upper)
    printContiguousRowsInTwoQuadrants(
        lowerLeftStrings, leftColWidths, 
        lowerRightStrings, rightColWidths, 
        lowerRowLabels,
        baseIndent, indentPerRow, upperLeft.size() + 1); // +1 for ellipsis
}



/*
 * FUNCTIONS TO DIVIDE A MATRIX INTO QUADRANTS
 *
 * as per the user's printing truncation thresholds.
 * These functions decide how to divide a matrix or
 * vector into one, two or four quadrants (or for
 * vectors, at most two quadrants) for subsequent
 * pretty printing
 */


struct MatrixQuadrantInds {

    qindex numUpperLeftRows=0, numUpperRightRows=0; 
    qindex numLowerLeftRows=0, numLowerRightRows=0;

    qindex numUpperLeftCols=0, numUpperRightCols=0;
    qindex numLowerLeftCols=0, numLowerRightCols=0;

    qindex rightStartCol=0;
    qindex lowerStartRow=0;

    // these are always fixed at 0, but included for clarity
    qindex leftStartCol = 0;
    qindex upperStartRow = 0;
};


MatrixQuadrantInds getTruncatedMatrixQuadrantInds(qindex numRows, qindex numCols) {

    MatrixQuadrantInds inds;

    // find maximum size of matrix quadrants according to user-truncatins.
    // Choose right & lower first so that When num=odd, extra left & upper elem is shown
    qindex maxNumRightCols  = global_maxNumPrintedCols / 2; // floors
    qindex maxNumLeftCols = global_maxNumPrintedCols - maxNumRightCols;
    qindex maxNumLowerRows = global_maxNumPrintedRows / 2; // floors
    qindex maxNumUpperRows = global_maxNumPrintedRows - maxNumLowerRows;

    // may be ignored (and will unimportantly underflow when not truncating)
    inds.rightStartCol = numCols - maxNumRightCols;
    inds.lowerStartRow = numRows - maxNumLowerRows;

    // decide among four possible quadrant population configurations;
    // this code can be significantly shortened but we leave it explicit
    // for clarity, and so that *NumRows=0 <=> *NumCols=0

    if (numRows <= global_maxNumPrintedRows && numCols <= global_maxNumPrintedCols) {

        // upper left contains entire matrix
        inds.numUpperLeftRows = numRows;  inds.numUpperLeftCols = numCols;

    } else if (numRows <= global_maxNumPrintedRows) {

        // matrix divided into upper left and upper right (bottoms empty)
        inds.numUpperLeftRows  = numRows;  inds.numUpperLeftCols  = maxNumLeftCols;
        inds.numUpperRightRows = numRows;  inds.numUpperRightCols = maxNumRightCols;
    
    } else if (numCols <= global_maxNumPrintedCols) {

        // matrix divided into upper left and lower left (rights empty)
        inds.numUpperLeftRows = maxNumUpperRows;  inds.numUpperLeftCols = numCols;
        inds.numLowerLeftRows = maxNumLowerRows;  inds.numLowerLeftCols = numCols;

    } else {
        
        // matrix divided into into four quadrants
        inds.numUpperLeftRows  = maxNumUpperRows;  inds.numUpperLeftCols  = maxNumLeftCols;
        inds.numUpperRightRows = maxNumUpperRows;  inds.numUpperRightCols = maxNumRightCols;
        inds.numLowerLeftRows  = maxNumLowerRows;  inds.numLowerLeftCols  = maxNumLeftCols;
        inds.numLowerRightRows = maxNumLowerRows;  inds.numLowerRightCols = maxNumRightCols;
    }

    return inds;
}



/*
 * FUNCTIONS TO POPULATE THE QUADRANTS
 *
 * the sizes of which are decided by the above functions.
 * The below functions invoke communication and GPU-to-CPU 
 * copying as necessary to populate the quadrants on the 
 * root node (that which subsequently prints them)
 */


void allocateMatrixQuadrants(MatrixQuadrantInds inds, qcompmatr &ul, qcompmatr &ur, qcompmatr &ll, qcompmatr &lr) {

    util_tryAllocMatrix(ul, inds.numUpperLeftRows,  inds.numUpperLeftCols,  error_printerFailedToAllocTempMemory);
    util_tryAllocMatrix(ur, inds.numUpperRightRows, inds.numUpperRightCols, error_printerFailedToAllocTempMemory); // may be zero-size
    util_tryAllocMatrix(ll, inds.numLowerLeftRows,  inds.numLowerLeftCols,  error_printerFailedToAllocTempMemory); // may be zero-size
    util_tryAllocMatrix(lr, inds.numLowerRightRows, inds.numLowerRightCols, error_printerFailedToAllocTempMemory); // may be zero-size
}


void populateSingleColumnQcompmatr(qcompmatr &matr, qindex startInd, PauliStrSum sum) {

    // matr is single-column, and has been prior resized such that 
    // even before population, matr.size() is non-zero unless it is
    // intended to stay empty
    size_t numCoeffs = matr.size(); // numRows
    if (numCoeffs == 0)
        return;

    // serially copy the CPU-only, non-distributed coefficients
    for (size_t i=0; i<numCoeffs; i++)
        matr[i][0] = sum.coeffs[startInd + i];
}


// T can be any 1D type, i.e. Qureg (.isDensity=0), FullStateDiagMatr, DiagMatr1, DiagMatr2, DiagMatr
template <class T>
void populateSingleColumnQcompmatr(qcompmatr &matr, qindex startInd, T obj) {

    // matr is single-column, and has been prior resized such that 
    // even before population, matr.size() is non-zero unless it is
    // intended to stay empty
    size_t numAmps = matr.size(); // numRows
    if (numAmps == 0)
        return;

    // prepare temp memory where adjacent columns are adjacent elems (unlike matr)
    vector<qcomp> column;
    util_tryAllocVector(column, numAmps, error_printerFailedToAllocTempMemory);

    // populate temp memory, potentially invoking CPU-GPU copying and communication...
    if constexpr (util_isFullStateDiagMatr<T>()) {
        localiser_fullstatediagmatr_getElems(column.data(), obj, startInd, numAmps);
    } else if constexpr (util_isQuregType<T>()) {
        localiser_statevec_getAmps(column.data(), obj, startInd, numAmps);

    // or invoking nothing...
    } else if constexpr (util_isFixedSizeMatrixType<T>()) {
        column = vector<qcomp>(obj.elems, obj.elems + numAmps);

    // or invoking merely a GPU-CPU copy 
    } else {
        column = vector<qcomp>(obj.cpuElems, obj.cpuElems + numAmps);
        auto gpuPtr = util_getGpuMemPtr(obj);
        if (mem_isAllocated(gpuPtr))
            gpu_copyGpuToCpu(gpuPtr, column.data(), numAmps);
    }

    // overwrite matr; serial iteration is fine since we assume few elems printed
    for (size_t i=0; i<numAmps; i++)
        matr[i][0] = column[i];
}


void populateManyColumnQcompmatr(qcompmatr &matr, qindex startRow, qindex startCol, Qureg qureg) {

    // matr has been prior resized such that even before
    // population, matr.size() is non-zero unless it is
    // intended to stay empty
    if (matr.size() == 0)
        return;

    // prepare list of pointers to matr rows, compatible with localiser 2D input format
    vector<qcomp*> ptrs;
    util_tryAllocVector(ptrs, matr.size(), error_printerFailedToAllocTempMemory);
    for (size_t r=0; r<matr.size(); r++)
        ptrs[r] = matr[r].data();

    // overwrite matr (throug ptrs), which may trigger communication and GPU-to-CPU copying.
    // note this achieves node consensus; every node populates their matr, which is much more
    // demanding than only the root note, however we expect this is not a big deal; live printing
    // in distributed settings is exceptionally rare and will print small, tractable subsets
    localiser_densmatr_getAmps(ptrs.data(), qureg, startRow, startCol, matr.size(), matr[0].capacity());
}


// T can be qcomp** or qcomp(*)[]
template <typename T>
void populateManyColumnQcompmatrFromPtr(qcompmatr &matr, qindex startRow, qindex startCol, T elems) {

    if (matr.size() == 0)
        return;

    // serially copy elements; performance is no issue here
    for (size_t r=0; r<matr.size(); r++)
        for (size_t c=0; c<matr[r].size(); c++)
            matr[r][c] = elems[r+startRow][c+startCol];
}


// T can be 2D matrix type, i.e CompMatr, CompMatr1, CompMatr2, SuperOp
template <typename T>
void populateManyColumnQcompmatr(qcompmatr &matr, qindex startRow, qindex startCol, T obj) {

    if constexpr (util_isFixedSizeMatrixType<T>())
        populateManyColumnQcompmatrFromPtr(matr, startRow, startCol, obj.elems);

    else {
        // even though matrix GPU memory should be unchanged, we copy it
        // over to CPU to overwrite any changes the user may have done
        // to the CPU matrix; this is so that the displayed matrix is
        // definitely consistent with the simulated backend. Because
        // we expect square matrices to be small, we do not bother copying
        // the specifically displayed sub-matrix (like we do for Qureg and
        // FullStateDiagMatr); instead, we just copy over the whole matrix.
        if constexpr (util_isHeapMatrixType<T>()) {
            auto gpuPtr = util_getGpuMemPtr(obj);
            if (mem_isAllocated(gpuPtr))
                gpu_copyGpuToCpu(obj);
        }

        populateManyColumnQcompmatrFromPtr(matr, startRow, startCol, obj.cpuElems);
    }
}


void populateMatrixQuadrants(MatrixQuadrantInds inds, qcompmatr &ul, qcompmatr &ur, qcompmatr &ll, qcompmatr &lr, Qureg obj) {

    // dense matrices can populate as many as 4 quadrants
    if (obj.isDensityMatrix) {
        populateManyColumnQcompmatr(ul, inds.upperStartRow, inds.leftStartCol,  obj);
        populateManyColumnQcompmatr(ur, inds.upperStartRow, inds.rightStartCol, obj);
        populateManyColumnQcompmatr(ll, inds.lowerStartRow, inds.leftStartCol,  obj);
        populateManyColumnQcompmatr(lr, inds.lowerStartRow, inds.rightStartCol, obj);

    // but vectors always become 1 or 2 quadrants (the left, as a potentially split column)
    } else {
        populateSingleColumnQcompmatr(ul, inds.upperStartRow, obj);
        populateSingleColumnQcompmatr(ll, inds.lowerStartRow, obj);
    }
}


// T can be all 1D and 2D types, i.e. CompMatr/1/2, DiagMatr1/2, FullStateDiagMatr, SuperOp,
// and can furthermore be PauliStrSum (in order to print the real coefficients).
// Importantly, this excludes KrausMap which is handled explicitly, and annoyingly, Qureg,
// which must be handled seperately above (because it must runtime branch, not compile-time)
template <class T>
void populateMatrixQuadrants(MatrixQuadrantInds inds, qcompmatr &ul, qcompmatr &ur, qcompmatr &ll, qcompmatr &lr, T obj) {

    // dense matrices can populate as many as 4 quadrants
    if constexpr (util_isDenseMatrixType<T>()) {
        populateManyColumnQcompmatr(ul, inds.upperStartRow, inds.leftStartCol,  obj);
        populateManyColumnQcompmatr(ur, inds.upperStartRow, inds.rightStartCol, obj);
        populateManyColumnQcompmatr(ll, inds.lowerStartRow, inds.leftStartCol,  obj);
        populateManyColumnQcompmatr(lr, inds.lowerStartRow, inds.rightStartCol, obj);

    // but vectors always become 1 or 2 quadrants (the left, as a potentially split column)
    } else {
        populateSingleColumnQcompmatr(ul, inds.upperStartRow, obj);
        populateSingleColumnQcompmatr(ll, inds.lowerStartRow, obj);
    }
}



/*
 * PRINTING DENSE MATRICES
 *
 * These functions print only the elements, 
 * indented, and optionally rank labels, e.g.
 *    [rank 0] [rank 1]
 *     3.1     3.2i
 *     -3.     5.3
 */


void setColumnLabelsToRanks(vector<string> &labels, int &lastRank, qindex colIndOffset, Qureg qureg) {

    // we can afford to serially enumerate printed columns, since expected few
    for (size_t i=0; i<labels.size(); i++) {
        int colRank = util_getRankContainingColumn(qureg, i + colIndOffset);
        labels[i] = (colRank > lastRank)? getRankStr(colRank) : "";

        // lastRank is modified, even for caller
        lastRank = colRank;
    }
}


void populateDensityMatrixColumnLabels(vector<string> &leftColLabels, vector<string> &rightColLabels, Qureg qureg, MatrixQuadrantInds inds) {

    leftColLabels.resize(inds.numUpperLeftCols);
    rightColLabels.resize(inds.numUpperRightCols);

    int lastRank = -1;
    setColumnLabelsToRanks(leftColLabels,  lastRank, inds.leftStartCol,  qureg);
    setColumnLabelsToRanks(rightColLabels, lastRank, inds.rightStartCol, qureg);
}


// T can be 2D types, e.g. Qureg (.isDensityMatrix=1), CompMatr/1/2, SuperOp.
// importantly, it excludes KrausMap which is handled separately
template <class T> 
void printDenseSquareMatrix(T obj, string indent) {

    // determine the full matrix dimension
    qindex numRows;
    if constexpr (util_isDenseMatrixType<T>())
        numRows = obj.numRows;
    else
        numRows = powerOf2(obj.numQubits); // this would be wrong for SuperOp

    // work out which global amps are going to be printed, divided into quadrants
    MatrixQuadrantInds inds = getTruncatedMatrixQuadrantInds(numRows, numRows);

    // populate quadrants (which may involve a GPU to CPU copy and/or communication)
    qcompmatr ul,ur,ll,lr;
    allocateMatrixQuadrants(inds, ul,ur,ll,lr);
    populateMatrixQuadrants(inds, ul,ur,ll,lr, obj);

    // columns are only labelled if given a distributed Qureg (they become sender node)
    vector<string> leftColLabels, rightColLabels;
    if constexpr (util_isQuregType<T>())
        if (obj.isDistributed)
            populateDensityMatrixColumnLabels(leftColLabels, rightColLabels, obj, inds);

    // do not label rows nor successively indent them
    vector<string> rowLabels = {};
    string indentPerRow = "";

    printMatrixInFourQuadrants(
        ul, ur, ll, lr,
        leftColLabels, rightColLabels, rowLabels, rowLabels, 
        indent, indentPerRow, VDOTS_CHAR);
}


void print_elems(CompMatr1 obj, string indent) { printDenseSquareMatrix(obj, indent); }
void print_elems(CompMatr2 obj, string indent) { printDenseSquareMatrix(obj, indent); }
void print_elems(CompMatr  obj, string indent) { printDenseSquareMatrix(obj, indent); }
void print_elems(SuperOp   obj, string indent) { printDenseSquareMatrix(obj, indent); }

void print_elems(KrausMap map, string indent) {

    if (!comm_isRootNode())
        return;

    // we ignore the superoperator, and instead print every CPU-only Kraus map in-turn...
    for (qindex i=0; i<map.numMatrices; i++) {

        // via wrapping its 2D memory in a spoofed CompMatr and re-using print_elems(CompMatr)
        CompMatr spoof;
        spoof.numQubits = map.numQubits;
        spoof.numRows   = map.numRows;

        // heap fields are not created since not consulted
        spoof.isApproxUnitary   = nullptr;
        spoof.isApproxHermitian = nullptr;
        spoof.wasGpuSynced      = nullptr;

        spoof.cpuElems = map.matrices[i];
        spoof.cpuElemsFlat = nullptr; // not consulted
        spoof.gpuElemsFlat = nullptr; // consulted; printer checks GPU-alloc via non-null

        cout << indent << getKrausOperatorIndStr(i) << endl; // single indent
        print_elems(spoof, indent + indent); // double indent
    }
}

// print_elems(Qureg.isDensityMatrix=1) handled separately below



/*
 * PRINTING VECTORS AND DIAGONAL MATRICES
 *
 * These functions print only the elements, 
 * indented, and optionally rank and ket labels, e.g.
 *     3.1  [rank 0]
 *       3.4i
 *         5.6  [rank 1]
 *           0.5i
 */


// T can be Qureg (.isDistributed=1, .isDensityMatrix=0) or FullStateDiagMatr
template <class T>
void setRowLabelsToRanks(vector<string> &labels, int &lastRank, qindex rowIndOffset, T obj) {

    // we can afford to serially enumerate printed rows, since expected few
    for (size_t i=0; i<labels.size(); i++) {

        // show rank only when distinct from the previous amp's rank
        int rank = util_getRankContainingIndex(obj, i + rowIndOffset);
        if (rank == lastRank)
            continue;

        // add extra spacing if the label already contains a ket
        if (!labels[i].empty())
            labels[i] += string(SET_SPACE_BETWEEN_KET_AND_RANK, MATRIX_SPACE_CHAR);

        labels[i] += getRankStr(rank);

        // lastRank is modified, even for caller
        lastRank = rank;
    }
}


// T can be Qureg (.isDistributed=1, .isDensityMatrix=0) or FullStateDiagMatr
template <class T>
void populateDistributedVectorRowLabels(vector<string> &upperRowLabels, vector<string> &lowerRowLabels, T obj, MatrixQuadrantInds inds) {

    // these might already contain kets; resize preserves existing elems
    upperRowLabels.resize(inds.numUpperLeftRows);
    lowerRowLabels.resize(inds.numLowerLeftRows);

    // to ensure rank aligns nicely, pad existing labels with the widest elem (which is the last one)
    int len = (lowerRowLabels.empty()? upperRowLabels : lowerRowLabels).back().length();
    for (auto &label : upperRowLabels) label += string(len - label.length(), MATRIX_SPACE_CHAR);
    for (auto &label : lowerRowLabels) label += string(len - label.length(), MATRIX_SPACE_CHAR);

    int lastRank = -1;
    setRowLabelsToRanks(upperRowLabels, lastRank, inds.upperStartRow, obj);
    setRowLabelsToRanks(lowerRowLabels, lastRank, inds.lowerStartRow, obj);
}


// T can be any 1D type, e.g. Qureg (.isDensityMatrix=0), FullStateDiagMatr, DiagMatr/1/2
template <class T>
void printVector(T obj, string indent) {

    // divide vector into one of two column vectors (ul and ll; ur and lr stay empty)
    qcompmatr ul,ur,ll,lr;
    qindex numAmps = powerOf2(obj.numQubits);
    MatrixQuadrantInds inds = getTruncatedMatrixQuadrantInds(numAmps, 1);
    allocateMatrixQuadrants(inds, ul,ur,ll,lr);
    populateMatrixQuadrants(inds, ul,ur,ll,lr, obj);

    // rows are only labelled if we're printing a Qureg...
    vector<string> upperRowLabels = {};
    vector<string> lowerRowLabels = {};
    if constexpr (util_isQuregType<T>()) {
        upperRowLabels = getBasisKets(inds.upperStartRow, inds.numUpperLeftRows);
        lowerRowLabels = getBasisKets(inds.lowerStartRow, inds.numLowerLeftRows);
    }

    // and/or if the type is distributed
    if constexpr (util_isDistributableType<T>())
        if (obj.isDistributed)
            populateDistributedVectorRowLabels(upperRowLabels, lowerRowLabels, obj, inds);

    // columns are never labelled
    vector<string> colLabels = {};

    // only successively indent each row when vector represets a diagonal matrix
    bool indentEachRow = ! util_isQuregType<T>(); // since else diagonal
    string indentPerRow = (indentEachRow)? string(SET_SPACE_BETWEEN_DIAG_MATRIX_COLS, MATRIX_SPACE_CHAR) : ""; 
    string ellipsisChar = (indentEachRow)? DDOTS_CHAR : VDOTS_CHAR;

    // will print the vector in one or two quadrants (right two quadrants are unpopulated)
    printMatrixInFourQuadrants(
        ul, ur, ll, lr,
        colLabels, colLabels, upperRowLabels, lowerRowLabels, 
        indent, indentPerRow, ellipsisChar);
}


void print_elems(DiagMatr1 obj, string indent) { printVector(obj, indent); }
void print_elems(DiagMatr2 obj, string indent) { printVector(obj, indent); }
void print_elems(DiagMatr  obj, string indent) { printVector(obj, indent); }
void print_elems(FullStateDiagMatr obj, string indent) { printVector(obj, indent); }

void print_elems(Qureg qureg, string indent) {

    (qureg.isDensityMatrix)?
        printDenseSquareMatrix(qureg, indent):
        printVector(qureg, indent);
}



/*
 * PRINTING PAULI STRING SUMS
 *
 * printing only the elements (i.e. coefficients and
 * strings), with no meta-data. e.g.
 *    0.123  XXIXX
 *    1.23i  XYZXZ
 *    -1-6i  IIIII
 */


// we'll make use of these internal functions from paulis.cpp
extern int paulis_getPauliAt(PauliStr str, int ind);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStr* strings, qindex numStrings);


string getPauliStrAsAllQubitsString(PauliStr str, int numPaulis) {

    // avoid repeated allocations in below string concatenation
    string out = "";
    out.reserve(numPaulis);

    // ugly but adequate - like me (call me)
    for (int i=numPaulis-1; i>=0; i--) {
        int code = paulis_getPauliAt(str, i); // 0123
        out += global_pauliChars[code];       // IXYZ unless user-overriden
    }
    
    return out;
}


string getPauliStrAsIndexString(PauliStr str, int numPaulis) {

    string out = "";

    // prematurely optimise (hehe) by avoiding repeated allocations
    // induced by below string concatenation, since we can bound the
    // string length; each Pauli is 1 char, each index is max 2 digits
    // (<64) and each inter-Pauli space is 1 char, so <=4 chars per Pauli.
    int maxLen = numPaulis * 4; // <= 256 always
    out.reserve(maxLen);

    for (int i=0; i<numPaulis; i++) {
        int code = paulis_getPauliAt(str, i);

        // don't report identities
        if (code == 0)
            continue;

        // e.g. X12 
        out += global_pauliChars[code] + toStr(i) + PAULI_SPACE_CHAR;
    }

    // out includes a trailing PAULI_SPACE_CHAR which is fine,
    // since nothing else is ever post-printed on the same line
    return out;
}


string getPauliStrAsString(PauliStr str, int numPauliChars=0) {

    // numPauliChars (only used by all-qubits reporting format)
    // allows the caller to pad the printed PauliStr with left I,
    // useful when reporting a truncated PauliStrSum to ensure
    // the top and bottom partitions are aligned and equal width.
    // When not passed, it defaults to minimum padding
    if (numPauliChars == 0)
        numPauliChars = 1 + paulis_getIndOfLefmostNonIdentityPauli(str);

    /// @todo
    /// to be absolutely pedantic/defensively-designed, we should
    /// assert numPauliChars >= 1+paulis_getIndOfLefmostNonIdentityPauli
    /// to ensure a bug never occludes a PauliStr. Very low priority!

    switch(global_pauliStrFormatFlag) {
        case 0: return getPauliStrAsAllQubitsString(str, numPauliChars);
        case 1: return getPauliStrAsIndexString(str, numPauliChars);
    }

    /// @todo
    /// to be defensively-designed, we should throw a runtime error
    /// here because global_pauliStrFormatFlag was mutilated. This is
    /// a very low priority; so we just return a suspicious string!
    return "(PauliStr stringifying failed - please alert the QuEST devs)";
}


void print_elemsWithoutNewline(PauliStr str, string indent) {

    // only root node ever prints
    if (!comm_isRootNode())
        return;

    cout << indent << getPauliStrAsString(str);

    // beware that NO trailing newlines are printed so subsequent
    // calls (without calling print_newlines()) will run-on
}


void print_elems(PauliStrSum sum, string indent) {

    // divide the list of sum coefficients into one or two quadrants (ul and ll; ur and lr stay empty)
    qcompmatr ul,ur,ll,lr;
    MatrixQuadrantInds inds = getTruncatedMatrixQuadrantInds(sum.numTerms, 1);
    allocateMatrixQuadrants(inds, ul,ur,ll,lr);
    populateMatrixQuadrants(inds, ul,ur,ll,lr, sum);

    // find the widest Pauli string, which informs how wide all printed Pauli strings are
    int upperWidth = 1 + paulis_getIndOfLefmostNonIdentityPauli(sum.strings, ul.size());
    int lowerWidth = (ll.empty())? 0 : 1 + paulis_getIndOfLefmostNonIdentityPauli(&sum.strings[inds.lowerStartRow], ll.size());
    int width = std::max(upperWidth, lowerWidth);

    // set row labels to be the Pauli strings
    vector<string> upperLabels(ul.size());
    vector<string> lowerLabels(ll.size());
    for (qindex i=0; i<inds.numUpperLeftRows; i++)
        upperLabels[i] = getPauliStrAsString(sum.strings[i], width);
    for (qindex i=0; i<inds.numLowerLeftRows; i++)
        lowerLabels[i] = getPauliStrAsString(sum.strings[i + inds.lowerStartRow], width);

    // columns are not labelled
    vector<string> colLabels = {};

    // will print the vector in one or two quadrants (right two quadrants are unpopulated)
    printMatrixInFourQuadrants(
        ul, ur, ll, lr,
        colLabels, colLabels, upperLabels, lowerLabels, 
        indent, "", VDOTS_CHAR);
}



/*
 * PRINTING STRUCT HEADERS
 *
 * which is the first line of text printed by reporters
 * (like reportCompMatr) which details struct fields/info,
 * preceeding the printing of the data/elements therein. E.g.
 *   PauliStrSum (3 terms, 104 bytes):
 */


void print_label(string label) {

    // always contains trailing newline (platform agnostic)
    print(label + LABEL_SUFFIX + "\n");
}


void printHeader(string name, vector<string> notes) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    cout << name << HEADER_OPEN_BRACKET;

    // safely assume notes.size() > 0
    for (size_t i=0; i<notes.size() - 1; i++)
        cout << notes[i] << HEADER_DELIMITER;

    cout << notes.back() << HEADER_CLOSE_BRACKET;
    cout << endl;
}


template <class T>
void printMatrixHeader(T matr, size_t numBytes) {

    bool isDiag = ! util_isDenseMatrixType<T>();
    bool isGpu = util_isGpuAcceleratedMatrix(matr);
    bool isDistrib = util_isDistributedMatrix(matr);
    int numNodes = (isDistrib)? comm_getNumNodes() : 1;

    // print e.g. FullStateDiagMatr (8 qubits, 256 qcomps over 4 nodes, 1056 bytes per GPU):
    printHeader(util_getMatrixTypeName<T>(), {
        getNumQubitsStr(matr.numQubits), 
        getNumMatrixElemsStr(util_getMatrixDim(matr), isDiag, numNodes),
        getMemoryCostsStr(numBytes, isDistrib, isGpu)});
}

void print_header(CompMatr1 m, size_t numBytes) { printMatrixHeader(m, numBytes); }
void print_header(CompMatr2 m, size_t numBytes) { printMatrixHeader(m, numBytes); }
void print_header(CompMatr  m, size_t numBytes) { printMatrixHeader(m, numBytes); }
void print_header(DiagMatr1 m, size_t numBytes) { printMatrixHeader(m, numBytes); }
void print_header(DiagMatr2 m, size_t numBytes) { printMatrixHeader(m, numBytes); }
void print_header(DiagMatr  m, size_t numBytes) { printMatrixHeader(m, numBytes); }
void print_header(FullStateDiagMatr m, size_t numBytes) { printMatrixHeader(m, numBytes); }


void print_header(SuperOp op, size_t numBytes) {

    // superoperator is a non-diagonal non-distributed matrix, maybe GPU
    bool isDiag = false;
    bool isDistrib = false;
    bool isGpu = getQuESTEnv().isGpuAccelerated;
    int numNodes = 1;

    // print e.g. SuperOp (1 qubit, 4x4 qcomps, 296 bytes in GPU)
    printHeader("SuperOp", {
        getNumQubitsStr(op.numQubits), 
        getNumMatrixElemsStr(op.numRows, isDiag, numNodes), 
        getMemoryCostsStr(numBytes, isDistrib, isGpu)});
}


void print_header(KrausMap map, size_t numBytes) {

    // KrausMap is non-distributed, but its superoperator is GPU-accel if the environment is
    bool isDistrib = false;
    bool isGpu = getQuESTEnv().isGpuAccelerated;

    // KrausMap uses custom substrings, in lieu of getNumMatrixElemsStr(), to report multiple matrices
    qindex mDim = map.numRows;
    qindex sDim = map.superop.numRows;
    string matrStr = toStr(map.numMatrices) + " " + getProductStr(mDim,mDim) + " " + (map.numMatrices>1? "matrices" : "matrix");
    string supStr =  "1 " + getProductStr(sDim,sDim) + " superoperator";

    // print e.g. KrausMap (1 qubit, 3 4x4 matrices, 1 16x16 superoperator, 1024 bytes in GPU)
    printHeader("KrausMap", {
        getNumQubitsStr(map.numQubits),
        matrStr, supStr,
        getMemoryCostsStr(numBytes, isDistrib, isGpu)});
}


void print_header(PauliStrSum sum, size_t numBytes) {

    // PauliStrSum coeffs are always local and CPU-only
    bool isDistrib = false;
    bool isGpu = false;

    printHeader("PauliStrSum", {
        toStr(sum.numTerms) + " term" + ((sum.numTerms>1)? "s":""),
        getMemoryCostsStr(numBytes, isDistrib, isGpu)});
}


void print_header(Qureg qureg, size_t numBytesPerNode) {

    // print e.g. Qureg (8 qubit statevector, 256 qcomps over 4 nodes, 1056 bytes per GPU):
    printHeader("Qureg", {
        getQuregTypeStr(qureg),
        getNumMatrixElemsStr(powerOf2(qureg.numQubits), ! qureg.isDensityMatrix, qureg.numNodes), // Qureg resembles matrix
        getMemoryCostsStr(numBytesPerNode, qureg.isDistributed, qureg.isGpuAccelerated)});
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
    for (auto const& [key, value] : rows) {
        if (key.length() > maxWidth)
            maxWidth = key.length();

        // pedantically suppressing unused-variable warning
        // (while still mirroring below enumeration for clarity)
        (void) value;
    }

    // print table title (indented)
    cout << indent << getTableTitleStr(title) << endl;

    // pad left column (which is doubly indented) to align right column
    for (auto const& [key, value] : rows)
        cout 
            << indent << indent << left 
            << setw(maxWidth + MIN_SPACE_BETWEEN_TABLE_COLS) << setfill(TABLE_SPACE_CHAR) 
            << key << value << endl;
}


void print_table(string title, vector<tuple<string, long long int>> rows, string indent) {

    // convert all values to strings
    vector<tuple<string, string>> casted;
    for (auto const& [key, value] : rows)
        casted.push_back({key, std::to_string(value)});

    print_table(title, casted, indent);
}


void print_table(string title, string emptyPlaceholder, string indent) {

    // only root node prints
    if (!comm_isRootNode())
        return;

    // print indented table title and placeholder
    cout << indent << getTableTitleStr(title) << endl;
    cout << indent << indent << emptyPlaceholder << endl;
}
