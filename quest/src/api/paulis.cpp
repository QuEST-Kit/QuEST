/** @file
 * API definitions for creating and managing Pauli strings,
 * and weighted sums thereof.
 */

#include "quest/include/precision.h"
#include "quest/include/paulis.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/parser.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <iostream>
#include <vector>
#include <string>

using std::string;
using std::vector;



/*
 * PRIVATE CONSTANTS
 */


static const int MAX_NUM_PAULIS_PER_MASK = sizeof(PAULI_MASK_TYPE) * 8 / 2;
static const int MAX_NUM_PAULIS_PER_STR  = MAX_NUM_PAULIS_PER_MASK * 2;



/*
 * PRIVATE UTILITIES
 */


int getPauliFromMaskAt(PAULI_MASK_TYPE mask, int ind) {

    // get adjacent 2 bits at (ind+1, ind)
    return (mask >> (2*ind)) & 3;
}


bool didAnyAllocsFailOnAnyNode(PauliStrSum sum) {

    bool anyFail = (sum.strings == NULL) || (sum.coeffs == NULL);
    if (comm_isInit())
        anyFail = comm_isTrueOnAllNodes(anyFail);

    return anyFail;
}


void freeAllMemoryIfAnyAllocsFailed(PauliStrSum sum) {

    // do nothing if everything allocated successfully between all nodes
    if (!didAnyAllocsFailOnAnyNode(sum))
        return;

    // otherwise free every successful allocation (freeing NULL is legal)
    free(sum.strings);
    free(sum.coeffs);
}



/*
 * INTERNAL UTILITIES
 *
 * callable by other internal files but which are not exposed in the header.
 * Ergo other files must declare these functions as extern where needed.
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


vector<int> paulis_getSortedIndsOfNonIdentityPaulis(PauliStr str) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    vector<int> inds(0);
    inds.reserve(maxInd+1);

    for (int i=0; i<=maxInd; i++)
        if (paulis_getPauliAt(str, i) != 0)
            inds.push_back(i);

    return inds;
}


bool paulis_containsXOrY(PauliStr str) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(str);

    for (int i=0; i<maxInd; i++) {
        int pauli = paulis_getPauliAt(str, i);

        if (pauli == 1 || pauli == 2)
            return true;
    }

    return false;
}


vector<int> paulis_getTargsWithEitherPaulis(vector<int> targs, PauliStr str, int pauliA, int pauliB) {

    vector<int> subsetTargs(0);  subsetTargs.reserve(targs.size());

    for (int targ : targs) {
        int pauli = paulis_getPauliAt(str, targ);
        if (pauli == pauliA || pauli == pauliB)
            subsetTargs.push_back(targ);
    }

    return subsetTargs;
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
    return (PauliStr) {
        .lowPaulis = lowPaulis,
        .highPaulis = highPaulis
    };

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

    // automatically target the bottom-most qubits, preserving rightmost is least significant
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

    // memory needed, as is reported to user after alloc failure
    qindex numBytesStrings = numTerms * sizeof(PauliStr);
    qindex numBytesCoeffs  = numTerms * sizeof(qcomp);

    // create struct
    PauliStrSum out = {
        .strings = (PauliStr*) malloc(numBytesStrings), // NULL if failed
        .coeffs  = (qcomp*)    malloc(numBytesCoeffs),  // NULL if failed
        .numTerms = numTerms
    };

    // if either alloc failed, clear both before validation to avoid leak
    freeAllMemoryIfAnyAllocsFailed(out);
    validate_newPauliStrSumAllocs(out, numBytesStrings, numBytesCoeffs, __func__);

    // serially copy data over to new heap memory
    for (int i=0; i<numTerms; i++) {
        out.strings[i] = strings[i];
        out.coeffs[i] = coeffs[i];
    }

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

    bool rightIsLeastSig = true;
    string str = parser_loadFile(fn);
    return parser_validateAndParsePauliStrSum(str, rightIsLeastSig, __func__);
}


extern "C" PauliStrSum createPauliStrSumFromReversedFile(const char* fn) {

    // fn must be null-terminated
    return createPauliStrSumFromReversedFile(string(fn));
}

PauliStrSum createPauliStrSumFromReversedFile(string fn) {
    validate_canReadFile(fn, __func__);

    bool rightIsLeastSig = false;
    string str = parser_loadFile(fn);
    return parser_validateAndParsePauliStrSum(str, rightIsLeastSig, __func__);
}



/*
 * DESTROYERS
 */


extern "C" void destroyPauliStrSum(PauliStrSum sum) {
    validate_pauliStrSumFields(sum, __func__);

    free(sum.strings);
    free(sum.coeffs);
}



/*
 * API REPORTERS
 */


extern "C" void reportPauliStr(PauliStr str) {

    // avoid printing leftmost superfluous I operators
    int numPaulis = 1 + paulis_getIndOfLefmostNonIdentityPauli(str);
    print_pauliStr(str, numPaulis);
}


extern "C" void reportPauliStrSum(PauliStrSum str) {
    validate_pauliStrSumFields(str, __func__);

    // calculate memory usage
    qindex numStrBytes   = str.numTerms * sizeof *str.strings;
    qindex numCoeffBytes = str.numTerms * sizeof *str.coeffs;
    qindex numStrucBytes = sizeof(str);

    // we don't bother checking for overflow since total memory scales
    // linearly with user input parameters, unlike Qureg and matrices.
    qindex numTotalBytes = numStrBytes + numCoeffBytes + numStrucBytes;

    print_pauliStrSumInfo(str.numTerms, numTotalBytes);
    print_pauliStrSum(str);
}
