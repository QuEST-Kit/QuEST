/** @file
 * Definitions of PauliStr and PauliStrSum,
 * their initialisers, and reporting utilities.
 * 
 * @author Tyson Jones
 * 
 * @defgroup paulis Paulis
 * @ingroup api
 * @brief Data structures for representing Pauli strings and their weighted sums.
 * @{
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
 * unlike some other headers, we here intermix the C and C++-only
 * signatures, grouping them semantically & by their doc groups
 */



/*
 * PAULI STRUCTS
 *
 * which are visible to both C and C++, and don't require demangling.
 * PauliStr contain only stack primitives, while PauliStrSum
 * contains dynamic heap pointers. Notice that PauliStr has non-const
 * members, because users will typically store large collections of
 * them (like in a PauliStrSum) so we wish to retain copy overwriting.
 */


/** 
 * @defgroup paulis_structs Structs
 * @brief Data structures for representing tensors and weighted sums of Pauli operators
 * @{
 */


/// @notyetdoced
typedef struct {

    // represent Pauli strings as base-4 numerals, split into their
    // upper and lower halves (as max-bit unsigned integers). This
    // imposes a strict upperbound on the number of stored Paulis.
    PAULI_MASK_TYPE lowPaulis;
    PAULI_MASK_TYPE highPaulis;

} PauliStr;


/// @notyetdoced
typedef struct {

    qindex numTerms;

    // arbitrarily-sized collection of Pauli strings and their
    // coefficients are stored in heap memory.
    PauliStr* strings;
    qcomp* coeffs;

    // whether the sum constitutes a Hermitian operator (0, 1, or -1 to indicate unknown),
    // which is lazily evaluated when a function validates Hermiticity them. The flag is 
    // stored in heap so even copies of structs are mutable, but the pointer is immutable;
    // otherwise, the field of a user's struct could never be modified because of pass-by-copy.
    int* isApproxHermitian;

} PauliStrSum;


/** @} */



// we define the remaining doc groups in advance, since their signatures are
// more naturally grouped in an implementation-specific way below. Note the
// above structs were not doc'd this way (which would be more consistent)
// because it inexplicably causes Doxygen to duplicate their section at the
// top-level under Paulis (rather than under Structs). Bizarre! The order
// of declaration below will match the order shown in the html doc.
/** 
 * @defgroup paulis_create Constructors
 * @brief Functions for creating and initialising Pauli data structures.
 * 
 * @defgroup paulis_destroy Destructors
 * @brief Functions for destroying existing Pauli data structures.
 * 
 * @defgroup paulis_reporters Reporters
 * @brief Functions for printing Pauli data structures.
 */



/*
 * PAULI STRING CREATION
 */


// base method is C and C++ compatible
#ifdef __cplusplus
extern "C" {
#endif

    /** @ingroup paulis_create
     * @notyetdoced
     * 
     * @see
     * - reportPauliStr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStr getPauliStr(const char* paulis, int* indices, int numPaulis);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus

    // C++ users can access the above C method, along with direct overloads 
    // to accept integers (in lieu of chars), natural C++ string types
    // (like literals), and C++ vector types for brevity. Furthermore, C++
    // gets an overload which accepts only a string (no additional args)
    // which is used internally by parsers, and exposed for user-convenience.
    // note that C++ does NOT get an overload where the pauli codes (integers) are
    // passed as a vector; this is because integer literal initialiser lists like
    // {0,3,1} are valid std::string instances, causing overload ambiguity. Blegh!


    /** @ingroup paulis_create
     * @notyetdoced
     * 
     * @see
     * - reportPauliStr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStr getPauliStr(int* paulis, int* indices, int numPaulis);


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - getPauliStr()
     * - reportPauliStr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStr getPauliStr(std::string paulis, int* indices, int numPaulis);


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportPauliStr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStr getPauliStr(std::string paulis, std::vector<int> indices);


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportPauliStr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStr getPauliStr(std::string paulis);


    // never needs to be doc'd
    /// @private
    /// @neverdoced
    #define getInlinePauliStr(str, ...) \
        getPauliStr(str, __VA_ARGS__)


#else

    // C supports passing a char array or string literal with a specified number of Paulis,
    // or an overload accepting ints, achieved using a C11 _Generic. C also gets an inline 
    // macro which exploits the compile-time size of a string literal, and enables array
    // literals without the C99 inline temporary array syntax; we further give the array
    // an explici size (instead of just (int[])) to gaurantee it contains at leasts as 
    // many elements as claimed, avoiding seg-faults if the user provides too few indices 


    /// @ingroup paulis_create
    /// @private
    PauliStr _getPauliStrFromInts(int* paulis, int* indices, int numPaulis);


    // documented above (identical signatures to C)
    /// @neverdoced
    #define getPauliStr(paulis, ...) \
        _Generic((paulis), \
            int*    : _getPauliStrFromInts, \
            default : getPauliStr \
        )(paulis, __VA_ARGS__) 


    // documented below
    /// @neverdoced
    #define getInlinePauliStr(str, ...) \
        getPauliStr((str), (int[sizeof(str)-1]) __VA_ARGS__, sizeof(str)-1)

    // spoofing above macro as function to doc
    #if 0

        /** @ingroup paulis_create
         * @notyetdoced
         * @macrodoc
         * 
         * @see
         * - reportPauliStr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) and
         *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp)examples
         */
        PauliStr getInlinePauliStr(const char* paulis, { list });

    #endif


#endif



/*
 * PAULI STRING SUM CREATION
 */


// base methods are C and C++ compatible
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup paulis_create
     * @notyetdoced
     * 
     * @see
     * - reportPauliStrSum()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createPauliStrSum(PauliStr* strings, qcomp* coeffs, qindex numTerms);


    /** @ingroup paulis_create
     * @notyetdoced
     * 
     * @see
     * - reportPauliStrSum()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createInlinePauliStrSum(const char* str);


    /** @ingroup paulis_create
     * @notyetdoced
     * 
     * @see
     * - reportPauliStrSum()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createPauliStrSumFromFile(const char* fn);


    /** @ingroup paulis_create
     * @notyetdoced
     * 
     * @see
     * - reportPauliStrSum()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createPauliStrSumFromReversedFile(const char* fn);


#ifdef __cplusplus
}
#endif


// C++ users get additional overloads
#ifdef __cplusplus


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createPauliStrSum()
     * - reportPauliStrSum()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createPauliStrSum(std::vector<PauliStr> strings, std::vector<qcomp> coeffs);


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createInlinePauliStrSum()
     * - reportPauliStrSum()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createInlinePauliStrSum(std::string str);


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createPauliStrSumFromFile()
     * - reportPauliStrSum()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createPauliStrSumFromFile(std::string fn);


    /** @ingroup paulis_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createPauliStrSumFromReversedFile()
     * - reportPauliStrSum()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_paulis.cpp) examples
     */
    PauliStrSum createPauliStrSumFromReversedFile(std::string fn);


#endif



/*
 * PAULI STRING SUM DESTRUCTION
 */


// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


    /// @ingroup paulis_destroy
    /// @notyetdoced
    void destroyPauliStrSum(PauliStrSum sum);


// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * REPORTERS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup paulis_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_paulis.cpp) examples
     */
    void reportPauliStr(PauliStr str);


    /** @ingroup paulis_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_paulis.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_paulis.cpp) examples
     */
    void reportPauliStrSum(PauliStrSum str);


// end de-mangler
#ifdef __cplusplus
}
#endif



#endif // PAULIS_H

/** @} */ // (end file-wide doxygen defgroup)
