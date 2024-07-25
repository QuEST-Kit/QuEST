/** @file
 * API definitions for creating and managing Pauli strings,
 * and weighted sums thereof.
 */

#include "quest/include/precision.h"
#include "quest/include/paulis.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/parser.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <iostream>
#include <vector>
#include <string>



/*
 * PRIVATE CONSTANTS
 */


static const int MAX_NUM_PAULIS_PER_MASK = sizeof(PAULI_MASK_TYPE) * 8 / 2;
static const int MAX_NUM_PAULIS_PER_STR  = MAX_NUM_PAULIS_PER_MASK * 2;



/*
 * PRIVATE UTILITIES
 */


int getPauliAt(PAULI_MASK_TYPE mask, int ind) {

    // get adjacent 2 bits at (ind+1, ind)
    return (mask >> (2*ind)) & 3;
}

int getPauliAt(PauliStr str, int ind) {

    return (ind < MAX_NUM_PAULIS_PER_MASK)?
        getPauliAt(str.lowPaulis,  ind) :
        getPauliAt(str.highPaulis, ind - MAX_NUM_PAULIS_PER_MASK);
}


int getIndOfLefmostPauli(PauliStr str) {

    int ind   = (str.highPaulis == 0)? 0 : MAX_NUM_PAULIS_PER_MASK;
    auto mask = (str.highPaulis == 0)? str.lowPaulis : str.highPaulis;

    while (mask) {
        mask >>= 2;
        ind++;
    }

    return ind - 1;
}


/*
 * PAULI STRING INITIALISATION
 *
 * some of which are exposed directly to C, and some of which are C++-only overloads
 */


extern "C" PauliStr getPauliStr(char* paulis, int* indices, int numPaulis) {
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

PauliStr getPauliStr(std::string paulis, int* indices, int numPaulis) {

    // additionally validate 'paulis' string has 'numPaulis' chars
    validate_newPauliStrNumChars(paulis.length(), numPaulis, __func__);

    return getPauliStr(paulis.data(), indices, numPaulis); // validates
}

PauliStr getPauliStr(std::string paulis, std::vector<int> indices) {

    // additionally validate 'paulis' string has 'numPaulis' chars
    validate_newPauliStrNumChars(paulis.length(), indices.size(), __func__);

    return getPauliStr(paulis.data(), indices.data(), indices.size()); // validates
}

PauliStr getPauliStr(std::string paulis) {

    // pedantically validate the string length isn't so long that it would stackoverflow a vector
    validate_newPauliStrNumPaulis(paulis.size(), MAX_NUM_PAULIS_PER_STR, __func__);

    // automatically target the bottom-most qubits, preserving rightmost is least significant
    std::vector<int> indices(paulis.size());
    for (int i=0; i<paulis.size(); i++)
        indices[i] = paulis.size() - 1 - i;

    return getPauliStr(paulis, indices); // validates
}



/*
 * API REPORTERS
 */


extern "C" void reportPauliStr(PauliStr str) {

    // avoid printing leftmost superfluous I operators
    int numPaulis = 1 + getIndOfLefmostPauli(str);
    form_printPauliStr(str, numPaulis);
}


