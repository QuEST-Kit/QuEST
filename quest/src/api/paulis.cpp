/** @file
 * API functions for creating PauliStr and PauliStrSum,
 * and initialising and reporting them
 * 
 * @author Tyson Jones
 */

#include "quest/include/precision.h"
#include "quest/include/paulis.h"

#include "quest/src/core/paulilogic.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/parser.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"

#include <vector>
#include <string>

using std::string;
using std::vector;



/*
 * PRIVATE UTILITIES
 */


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

    // return a new stack PauliStr instance (avoiding C++20 initialiser)
    PauliStr out;
    out.lowPaulis = lowPaulis;
    out.highPaulis = highPaulis;
    return out;
}


PauliStr getPauliStr(int* paulis, int* indices, int numPaulis) {
    validate_newPauliStrParams(paulis, indices, numPaulis, MAX_NUM_PAULIS_PER_STR, __func__);

    // validation ensures never causes stack overflow
    char pauliChars[MAX_NUM_PAULIS_PER_STR + 1]; // +1 for null-terminal

    // make a char array from the pauli codes, using an arbitrary 
    // choice of the Pauli characters accepted by the API (like IXYZ)
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

    // prepare output PauliStrSum (avoiding C++20 designated initialiser)
    PauliStrSum out;
    out.numTerms = numTerms;
    out.strings = cpu_allocPauliStrings(numTerms);                // nullptr if failed
    out.coeffs  = cpu_allocArray(numTerms);                       // nullptr if failed
    out.isApproxHermitian = util_allocEpsilonSensitiveHeapFlag(); // nullptr if failed

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
