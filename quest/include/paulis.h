/** @file
 * Signatures of API Pauli string and sum data structures, and their getters and setters, 
 * as well as reporting utilities.
 */

#ifndef PAULIS_H
#define PAULIS_H

#include "quest/include/precision.h"
#include "quest/include/types.h"

// C++ gets string and vector initialiser overloads
#ifdef __cplusplus
    #include <string>
    #include <vector>
#endif



/*
 * PAULI STRUCTS
 *
 * which are visible to both C and C++, and don't require demangling.
 */


typedef struct {

    // represent Pauli strings as base-4 numerals, split into their
    // upper and lower halves (as max-bit unsigned integers). This
    // imposes a strict upperbound on the number of stored Paulis.
    PAULI_MASK_TYPE lowPaulis;
    PAULI_MASK_TYPE highPaulis;

} PauliStr;


/*
 * PAULI STRING CREATION
 */


#ifdef __cplusplus

    // C++ users can access the base C method
    extern "C" PauliStr getPauliStr(char* paulis, int* indices, int numPaulis);

    // but also get overloads to accept natural C++ string types (like literals)
    PauliStr getPauliStr(std::string paulis, int* indices, int numPaulis);

    // and an additional overload to use vectors for brevity
    PauliStr getPauliStr(std::string paulis, std::vector<int> indices);

    // and an overload assuming indices={n,...,2,1,0}, used internally by parsers
    PauliStr getPauliStr(std::string paulis);

    // inline macro (included for consistency with C) calls the string & vector overload
    #define getInlinePauliStr(str, ...) \
        getPauliStr(str, __VA_ARGS__)

#else

    // C only supports passing a char array or string literal with a specified number of Paulis
    PauliStr getPauliStr(char* paulis, int* indices, int numPaulis);

    // inline macro exploits the compile-time size of a string literal, and enables array
    // literals without the C99 inline temporary array syntax; we further give the array
    // an explici size (instead of just (int[])) to gaurantee it contains at leasts as 
    // many elements as claimed, avoiding seg-faults if the user provides too few indices 
    #define getInlinePauliStr(str, ...) \
        getPauliStr((str), (int[sizeof(str)-1]) __VA_ARGS__, sizeof(str)-1)

#endif



/*
 * REPORTERS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif

    void reportPauliStr(PauliStr str);

// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // PAULIS_H