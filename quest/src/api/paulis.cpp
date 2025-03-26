/** @file
 * API functions for creating PauliStr and PauliStrSum,
 * and initialising and reporting them
 * 
 * @author Tyson Jones
 */

#include "quest/include/precision.h"
#include "quest/include/paulis.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/parser.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <array>

using std::string;
using std::vector;
using std::array;



/*
 * PRIVATE CONSTANTS
 */


static const int MAX_NUM_PAULIS_PER_MASK = sizeof(PAULI_MASK_TYPE) * 8 / 2;
static const int MAX_NUM_PAULIS_PER_STR  = MAX_NUM_PAULIS_PER_MASK * 2;



/*
 * PRIVATE UTILITIES
 */


int getPauliFromMaskAt(PAULI_MASK_TYPE mask, int ind) {

    return getTwoAdjacentBits(mask, 2*ind); // bits at (ind+1, ind)
}


bool didAnyAllocsFailOnAnyNode(PauliStrSum sum) {

    bool anyFail = (
        ! mem_isAllocated(sum.strings) || 
        ! mem_isAllocated(sum.coeffs)  || 
        ! mem_isAllocated(sum.isApproxHermitian) );
    
    if (comm_isInit())
        anyFail = comm_isTrueOnAllNodes(anyFail);

    return anyFail;
}


void freePauliStrSum(PauliStrSum sum) {

    // these do not need to be allocated (freeing nullptr is legal)
    cpu_deallocPauliStrings(sum.strings);
    cpu_deallocArray(sum.coeffs);
    util_deallocEpsilonSensitiveHeapFlag(sum.isApproxHermitian);
}


void freeAllMemoryIfAnyAllocsFailed(PauliStrSum sum) {

    // do nothing if everything allocated successfully between all nodes
    if (!didAnyAllocsFailOnAnyNode(sum))
        return;

    // otherwise free every successful allocation (freeing nullptr is legal)
    freePauliStrSum(sum);
}



/*
 * INTERNAL UTILITIES
 *
 * callable by other internal files but which are not exposed in the header
 * because we do not wish to make them visible to users. Ergo other internal
 * files must declare these functions as extern where needed. Yes, it's ugly :(
 */


int paulis_getPauliAt(PauliStr str, int ind) {

    return (ind < MAX_NUM_PAULIS_PER_MASK)?
        getPauliFromMaskAt(str.lowPaulis,  ind) :
        getPauliFromMaskAt(str.highPaulis, ind - MAX_NUM_PAULIS_PER_MASK);
}


int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str) {

    int ind   = (str.highPaulis == 0)? 0 : MAX_NUM_PAULIS_PER_MASK;
    auto mask = (str.highPaulis == 0)? str.lowPaulis : str.highPaulis;

    while (mask) {
        mask >>= 2;
        ind++;
    }

    return ind - 1;
}


int paulis_getIndOfLefmostNonIdentityPauli(PauliStr* strings, qindex numStrings) {

    int maxInd = 0;

    for (qindex i=0; i<numStrings; i++) {
        int ind = paulis_getIndOfLefmostNonIdentityPauli(strings[i]);
        if (ind > maxInd)
            maxInd = ind;
    }

    return maxInd;
}


int paulis_getIndOfLefmostNonIdentityPauli(PauliStrSum sum) {

    return paulis_getIndOfLefmostNonIdentityPauli(sum.strings, sum.numTerms);
}


bool paulis_containsXOrY(PauliStr str) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    for (int i=0; i<=maxInd; i++) {
        int pauli = paulis_getPauliAt(str, i);

        if (pauli == 1 || pauli == 2)
            return true;
    }

    return false;
}


bool paulis_containsXOrY(PauliStrSum sum) {

    for (qindex i=0; i<sum.numTerms; i++)
        if (paulis_containsXOrY(sum.strings[i]))
            return true;

    return false;
}


bool paulis_hasOddNumY(PauliStr str) {

    bool odd = false;

    for (int targ=0; targ < MAX_NUM_PAULIS_PER_STR; targ++) 
        if (paulis_getPauliAt(str, targ) == 2)
            odd = !odd;

    return odd;
}


int paulis_getPrefixZSign(Qureg qureg, vector<int> prefixZ) {

    int sign = 1;

    // each Z contributes +- 1
    for (int qubit : prefixZ)
        sign *= util_getRankBitOfQubit(qubit, qureg)? -1 : 1;

    return sign;
}


qcomp paulis_getPrefixPaulisElem(Qureg qureg, vector<int> prefixY, vector<int> prefixZ) {

    // each Z contributes +- 1
    qcomp elem = paulis_getPrefixZSign(qureg, prefixZ);

    // each Y contributes -+ i
    for (int qubit : prefixY)
        elem *= 1_i * (util_getRankBitOfQubit(qubit, qureg)? 1 : -1);

    return elem;
}


vector<int> paulis_getInds(PauliStr str) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    vector<int> inds(0);
    inds.reserve(maxInd+1);

    for (int i=0; i<=maxInd; i++)
        if (paulis_getPauliAt(str, i) != 0)
            inds.push_back(i);

    return inds;
}


array<vector<int>,3> paulis_getSeparateInds(PauliStr str, Qureg qureg) {

    vector<int> iXYZ = paulis_getInds(str);
    vector<int> iX, iY, iZ;

    vector<int>* ptrs[] = {&iX, &iY, &iZ};

    for (int i : iXYZ)
        ptrs[paulis_getPauliAt(str, i) - 1]->push_back(i);

    return {iX, iY, iZ};
}


PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift) {

    if (pauliShift <= 0 || pauliShift >= MAX_NUM_PAULIS_PER_MASK)
        error_pauliStrShiftedByIllegalAmount();

    int numBitsPerPauli = 2;
    int numMaskBits = numBitsPerPauli * MAX_NUM_PAULIS_PER_MASK;
    int bitShift    = numBitsPerPauli * pauliShift;

    // record the bits we will lose from lowPaulis, to move to highPaulis
    PAULI_MASK_TYPE lostBits = getBitsLeftOfIndex(str.lowPaulis, numMaskBits - bitShift - 1);

    // ensure we actually lose these bits from lowPaulis
    PAULI_MASK_TYPE lowerBits = getBitsRightOfIndex(str.lowPaulis, numMaskBits - bitShift) << bitShift;

    // and add them to highPaulis; we don't have to force lose upper bits of high paulis
    PAULI_MASK_TYPE upperBits = concatenateBits(str.highPaulis, lostBits, bitShift);

    return {
        .lowPaulis = lowerBits,
        .highPaulis = upperBits
    };
}


PauliStr paulis_getKetAndBraPauliStr(PauliStr str, Qureg qureg) {

    PauliStr shifted = paulis_getShiftedPauliStr(str, qureg.numQubits);
    
    return {
        .lowPaulis  = str.lowPaulis  | shifted.lowPaulis,
        .highPaulis = str.highPaulis | shifted.highPaulis
    };
}


PAULI_MASK_TYPE paulis_getKeyOfSameMixedAmpsGroup(PauliStr str) {

    PAULI_MASK_TYPE key = 0;

    // in theory, we can reduce the number of involved operations by bit-shifting
    // str left by 1, XOR'ing this with str, and retaining every 2nd bit, producing
    // e.g. key=0110 from str=IXYZ. However, this is an insignificant speedup which
    // risks sneaky bugs related to handling str's two masks.

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    for (int i=0; i<=maxInd; i++) {
        int pauli = paulis_getPauliAt(str, i);
        int isXY = (pauli == 1 || pauli == 2);
        key |= (isXY << i);
    }

    return key;
}



/*
 * PAULI STRING INITIALISATION
 *
 * some of which are exposed directly to C, and some of which are C++-only overloads
 */


extern "C" PauliStr getPauliStr(const char* paulis, int* indices, int numPaulis) {
    validate_newPauliStrParams(paulis, indices, numPaulis, MAX_NUM_PAULIS_PER_STR, __func__);

    // begin masks at all-identity 'I' = 0
    PAULI_MASK_TYPE lowPaulis = 0;
    PAULI_MASK_TYPE highPaulis = 0;

    // change targeted indices to the given Paulis
    for (int i=0; i<numPaulis; i++) {

        // cast single Pauli to full precision mask to enable below shifts
        auto pauli = (PAULI_MASK_TYPE) parser_getPauliIntFromChar(paulis[i]);

        // add the Pauli to either the lower or upper pauli masks
        if (indices[i] < MAX_NUM_PAULIS_PER_MASK)
            lowPaulis  |= pauli << (2*indices[i]);
        else
            highPaulis |= pauli << (2*(indices[i] - MAX_NUM_PAULIS_PER_MASK));
    }

    // return a new stack PauliStr instance, returning by copy
    return {
        .lowPaulis = lowPaulis,
        .highPaulis = highPaulis
    };
}


PauliStr getPauliStr(int* paulis, int* indices, int numPaulis) {
    validate_newPauliStrParams(paulis, indices, numPaulis, MAX_NUM_PAULIS_PER_STR, __func__);

    // validation ensures never causes stack overflow
    char pauliChars[MAX_NUM_PAULIS_PER_STR + 1]; // +1 for null-terminal

    // made a char array from the pauli codes
    for (int i=0; i<numPaulis; i++)
        pauliChars[i] = "IXYZ"[paulis[i]];

    // including the trailing null char, used to infer string end/length
    pauliChars[numPaulis] = '\0';

    return getPauliStr(pauliChars, indices, numPaulis);
}


extern "C" PauliStr _getPauliStrFromInts(int* paulis, int* indices, int numPaulis) {

    return getPauliStr(paulis, indices, numPaulis);
}


PauliStr getPauliStr(string paulis, int* indices, int numPaulis) {

    // additionally validate 'paulis' string has 'numPaulis' chars
    validate_newPauliStrNumChars(paulis.length(), numPaulis, __func__);

    return getPauliStr(paulis.data(), indices, numPaulis); // validates
}

PauliStr getPauliStr(string paulis, vector<int> indices) {

    // additionally validate 'paulis' string has 'numPaulis' chars
    validate_newPauliStrNumChars(paulis.length(), indices.size(), __func__);

    return getPauliStr(paulis.data(), indices.data(), indices.size()); // validates
}

PauliStr getPauliStr(string paulis) {

    // pedantically validate the string length isn't so long that it would stackoverflow a vector
    validate_newPauliStrNumPaulis(paulis.size(), MAX_NUM_PAULIS_PER_STR, __func__);

    // automatically target the lowest-index qubits, interpreting rightmost is least significant
    vector<int> indices(paulis.size());
    for (size_t i=0; i<paulis.size(); i++)
        indices[i] = paulis.size() - 1 - i;

    return getPauliStr(paulis, indices); // validates
}



/*
 * PAULI STRING SUM CREATION
 *
 * some of which are exposed directly to C, and some of which are C++-only overloads
 */


extern "C" PauliStrSum createPauliStrSum(PauliStr* strings, qcomp* coeffs, qindex numTerms) {

    // note we do not require nor impose the strings to be unique
    validate_newPauliStrSumParams(numTerms, __func__);

    // create struct
    PauliStrSum out = {
        .numTerms = numTerms,
        .strings = cpu_allocPauliStrings(numTerms),                // nullptr if failed
        .coeffs  = cpu_allocArray(numTerms),                       // nullptr if failed
        .isApproxHermitian = util_allocEpsilonSensitiveHeapFlag(), // nullptr if failed
    };

    // if either alloc failed, clear both before validation to avoid leak
    freeAllMemoryIfAnyAllocsFailed(out);
    validate_newPauliStrSumAllocs(out, numTerms*sizeof(PauliStr), numTerms*sizeof(qcomp), __func__);

    // otherwise copy given data into new heap structure, and set initial flags
    cpu_copyPauliStrSum(out, strings, coeffs);
    util_setFlagToUnknown(out.isApproxHermitian);

    return out;
}

PauliStrSum createPauliStrSum(vector<PauliStr> strings, vector<qcomp> coeffs) {

    // additionally validate 'strings' and 'coeffs' are the same length
    validate_newPauliStrSumMatchingListLens(strings.size(), coeffs.size(), __func__);

    return createPauliStrSum(strings.data(), coeffs.data(), coeffs.size()); // validates
}


extern "C" PauliStrSum createInlinePauliStrSum(const char* str) {

    // str must be null-terminated
    return createInlinePauliStrSum(string(str));
}

PauliStrSum createInlinePauliStrSum(string str) {

    bool rightIsLeastSig = true;
    return parser_validateAndParsePauliStrSum(str, rightIsLeastSig, __func__);
}


extern "C" PauliStrSum createPauliStrSumFromFile(const char* fn) {

    // fn must be null-terminated
    return createPauliStrSumFromFile(string(fn));
}

PauliStrSum createPauliStrSumFromFile(string fn) {
    validate_canReadFile(fn, __func__);

    // all distributed nodes will simultaneously read the file (that's fine)
    string str = parser_loadFile(fn);

    bool rightIsLeastSig = true;
    return parser_validateAndParsePauliStrSum(str, rightIsLeastSig, __func__);
}


extern "C" PauliStrSum createPauliStrSumFromReversedFile(const char* fn) {

    // fn must be null-terminated
    return createPauliStrSumFromReversedFile(string(fn));
}

PauliStrSum createPauliStrSumFromReversedFile(string fn) {
    validate_canReadFile(fn, __func__);

    // all distributed nodes will simultaneously read the file (that's fine)
    string str = parser_loadFile(fn);

    bool rightIsLeastSig = false;
    return parser_validateAndParsePauliStrSum(str, rightIsLeastSig, __func__);
}



/*
 * DESTROYERS
 */


extern "C" void destroyPauliStrSum(PauliStrSum sum) {
    validate_pauliStrSumFields(sum, __func__);

    freePauliStrSum(sum);
}



/*
 * API REPORTERS
 */


extern "C" void reportPauliStr(PauliStr str) {

    // no header, so no indentation
    string indent = "";
    print_elemsWithoutNewline(str, indent);

    // print all user-set newlines (including none)
    print_newlines();
}


extern "C" void reportPauliStrSum(PauliStrSum sum) {
    validate_pauliStrSumFields(sum, __func__);
    validate_numReportedNewlinesAboveZero(__func__);

    // calculate memory usage
    qindex numStrBytes   = sum.numTerms * sizeof *sum.strings;
    qindex numCoeffBytes = sum.numTerms * sizeof *sum.coeffs;
    qindex numStrucBytes = sizeof(sum);

    // we don't bother checking for overflow since total memory scales
    // linearly with user input parameters, unlike Qureg and matrices.
    qindex numTotalBytes = numStrBytes + numCoeffBytes + numStrucBytes;

    print_header(sum, numTotalBytes);
    print_elems(sum);
    
    // exclude mandatory newline above
    print_oneFewerNewlines();
}
