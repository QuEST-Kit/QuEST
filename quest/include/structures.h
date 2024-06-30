/** @file
 * Signatures of API data structures like gate matrices. 
 * Note QuESTEnv and Qureg structs have their own signatures 
 * in environment.h and qureg.h respectively.
 * 
 * This file uses extensive preprocessor trickery to achieve platform agnostic,
 * C and C++ compatible, type agnostic, getters and setters of complex matrices.
 * First, we define C and C++ "explicit" functions like getCompMatr1FromArr()
 * (though C must wrap these with wrappers.h due to bad qcomp interoperability).
 * These explicitly disambiguate the input types in the function name.
 * Next, we define getCompMatr1() as a generic C macro and C++ overloads, so
 * that users can agnostically pass 2D pointers, arrays of pointers, or 2D arrays,
 * including compound literals like (qcomp[2][2]) {{...}}. So far, only the C++ 
 * overloads accept literals without the declaration like {{...}} via std::vector.
 * So, we next define getLiteralCompMatr1() as a C and C++ (for consistency) macro 
 * which accepts direct, concise literals like getLiteralCompMatr1({{...}}). Viola! 
 * 
 * We employ similar tricks to make setCompMatrN(), but define setCompMatrNFromArr()
 * in this header (exposed only to C) because it requires C++-incompatible VLAs;
 * this definition must invoke a bespoke validation function (gross). Also, our
 * implementation of setLiteralCompMatrN has to resort to preprocessor stringifying 
 * and runtime parsing - blegh!
 * 
 * All of these "intermediate" CompMatr initialisations are exposed to users.
 * Keep in mind that definitions herein will lead to symbol duplication unless 
 * specified as 'static inline', but macros are fine.
 * 
 * You're moving in a land of both shadow and substance, of things and ideas.
 * You've just crossed over into the Twilight Zone.
 */

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "quest/include/types.h"

// C++ gets vector initialiser overloads, whereas C gets a macro
#ifdef __cplusplus
    #include <vector>
#endif



/*
 * MATRIX STRUCTS
 *
 * which are visible to both C and C++, where qcomp resolves
 * to the native complex type. These are not de-mangled because
 * C++ structs are already C compatible.
 */


typedef struct CompMatr1
{
    int numQubits;
    qindex numRows;
    qcomp elems[2][2];

} CompMatr1;


typedef struct CompMatr2
{
    int numQubits;
    qindex numRows;
    qcomp elems[4][4];

} CompMatr2;


typedef struct CompMatrN
{
    int numQubits;
    qindex numRows;
    qcomp** elems;

    // row-flattened elems in GPU memory, allocated only
    // in GPU-enabled QuEST environments (regardless of Quregs)
    qcomp* gpuElems;

} CompMatrN;



/*
 * EXPLICIT FIXED-SIZE MATRIX INITIALISERS
 *
 * which are exposed directly only to C++ because they are incompatile with C binaries
 * due to returning qcomp (or fixed size arrays thereof) by-value through their structs.
 * The incompatibility arises because C++'s' std::complex and C's complex.h type are
 * distinct in the application binary interface (ABI), so a C binary cannot directly call
 * a C++ function which receives/returns a std::complex and interpret it as a complex.h type.
 * Equivalent and identical (to the user) C-compatible definitions of below's functions are
 * defined in wrappers.h, and wrap alternate C definitions which pass/receive C++ complex
 * through pointers, rather than by value, which is fine because the C & C++ types have
 * the same memory layout.
 */


#ifdef __cplusplus

    CompMatr1 getCompMatr1FromArr(qcomp in[2][2]);
    CompMatr1 getCompMatr1FromPtr(qcomp** in);

    CompMatr2 getCompMatr2FromArr(qcomp in[4][4]);
    CompMatr2 getCompMatr2FromPtr(qcomp** in);

#endif



/*
 * OVERLOADED FIXED-SIZE MATRIX INITIALISERS
 */


#ifdef __cplusplus

    // C++ uses overloads, accepting even vector initialiser lists

    CompMatr1 getCompMatr1(qcomp in[2][2]);
    CompMatr1 getCompMatr1(qcomp** in);
    CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in);

    CompMatr2 getCompMatr2(qcomp in[4][4]);
    CompMatr2 getCompMatr2(qcomp** in);
    CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in);

#else

    // C uses C11 compile-time type inspection, but literals must have compound syntax.
    // Explicit qcomp[2][2] type isn't necessary because it decays to qcomp(*)[2].
    // Use of __VA_ARGS__ is necessary to accept multiple-token compound literals.
    // Sadly we cannot use _Generic 'default' to catch unrecognised types at compile time.

    // Using type qcomp(*) below would erroneously invoke the qcomp(re,im) macro.
    // Preventing expansion using (qcomp)(*) leads to _Generic not recognising the type.
    // So, in desperation, we alias the type qcomp with a name that has no colliding macro.
    typedef qcomp qalias;

    #define getCompMatr1(...) _Generic((__VA_ARGS__), \
        qalias(*)[2] : getCompMatr1FromArr, \
        qcomp**      : getCompMatr1FromPtr  \
        )((__VA_ARGS__))

    #define getCompMatr2(...) _Generic((__VA_ARGS__), \
        qalias(*)[4] : getCompMatr2FromArr, \
        qcomp**      : getCompMatr2FromPtr  \
        )((__VA_ARGS__))

#endif



/*
 * LITERAL FIXED-SIZE MATRIX INITIALISERS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * compound literal syntax. We expose these macros to C++ too for API consistency.
 * although C++'s getCompMatr1 vector overload achieves the same thing, and cannot
 * use C-style temporary arrays.
 */

#ifdef __cplusplus

    // C++ merely invokes the std::vector initialiser overload

    #define getLiteralCompMatr1(...) \
        getCompMatr1(__VA_ARGS__)

    #define getLiteralCompMatr2(...) \
        getCompMatr2(__VA_ARGS__)

#else

    // C adds compound literal syntax to make a temporary array

    #define getLiteralCompMatr1(...) \
        getCompMatr1FromArr((qcomp[2][2]) __VA_ARGS__)

    #define getLiteralCompMatr2(...) \
        getCompMatr2FromArr((qcomp[4][4]) __VA_ARGS__)

#endif



/*
 * VARIABLE-SIZE MATRIX CONSTRUCTORS
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    CompMatrN createCompMatrN(int numQubits);

    void destroyCompMatrN(CompMatrN matrix);

    void syncCompMatrN(CompMatrN matr);

#ifdef __cplusplus
}
#endif



/*
 * EXPLICIT VARIABLE-SIZE MATRIX INITIALISERS
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    void setCompMatrNFromPtr(CompMatrN matr, qcomp** vals);

#ifdef __cplusplus
}
#endif


// permit only C compilers to have a VLA version (not supported nor needed by C++)
#ifndef __cplusplus

    // expose this function's bespoke validation
    extern void validate_setCompMatrNFromArr(CompMatrN out);

     // static inline to avoid header-symbol duplication
    static inline void setCompMatrNFromArr(CompMatrN matr, qcomp arr[matr.numRows][matr.numRows]) {

        // this function will allocate stack memory of size matr.numRows, but that field could
        // be invalid since matr hasn't been validated, so we must invoke bespoke validation
        validate_setCompMatrNFromArr(matr);

        // new ptrs array safely fits in stack, since it's sqrt-smaller than user's passed stack array
        qcomp* ptrs[matr.numRows];

        // collect pointers to each row of arr
        for (qindex r=0; r<matr.numRows; r++)
            ptrs[r] = arr[r];

        // array decays to qcomp**, and *FromPtr function re-performs validation (eh)
        setCompMatrNFromPtr(matr, ptrs);
    }

#endif



/*
 * OVERLOADED VARIABLE-SIZE MATRIX INITIALISERS
 */


#ifdef __cplusplus

    // C++ uses overloads, accepting even vector initialiser lists, but cannot ever accept 2D arrays

    void setCompMatrN(CompMatrN out, qcomp** in);

    void setCompMatrN(CompMatrN out, std::vector<std::vector<qcomp>> in);

#else

    // C uses C11 compile-time type inspection, but literals must have compound syntax.
    // Explicit qcomp[][] type isn't necessary because it decays to qcomp(*)[].
    // Use of __VA_ARGS__ is necessary to accept multiple-token compound literals.
    // Sadly we cannot use _Generic 'default' to catch unrecognised types at compile time.

    // Using type qcomp(*) below would erroneously invoke the qcomp(re,im) macro.
    // Preventing expansion using (qcomp)(*) leads to _Generic not recognising the type.
    // So, in desperation, we re-use the qcomp alias made by the previous _Generic use.

    #define setCompMatrN(matr, ...) _Generic((__VA_ARGS__),   \
        qalias(*)[] : setCompMatrNFromArr, \
        qcomp**     : setCompMatrNFromPtr  \
        )((matr), (__VA_ARGS__))

#endif



/*
 * STRING PARSING VARIABLE-SIZE MATRIX INITIALISERS
 */


// de-mangle for both C and C++ invocation
#ifdef __cplusplus
extern "C" {
#endif

    void setStringCompMatrN(CompMatrN matr, char* string);

#ifdef __cplusplus
}
#endif



/*
 * LITERAL VARIABLE-SIZE MATRIX INITIALISERS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s getCompMatr1 vector overload achieves the same thing
 */


#ifdef __cplusplus

    // C++ gets an explicit redirect to setCompMatrN(std::vector...), since it's faster than stringifying

    #define setLiteralCompMatrN(matr, ...) \
        setCompMatrN(matr, __VA_ARGS__)

#else 

    // C VLAs cannot be initialized, so we cannot use the same trick as used for the fixed-size
    // arrays, and i.e. invoke setCompMatrNFromArr() with inline array (qcomp(*)[matr.numRows]) __VA_ARGS__.
    // Instead, we must diabolically convert the expression to a compile-time string and parse it. Eep!

    #define setLiteralCompMatrN(matr, ...) \
        setStringCompMatrN(matr, #__VA_ARGS__)

#endif



/*
 * MATRIX REPORTERS
 */


// de-mangle so below are directly callable by C binary
#ifdef __cplusplus
extern "C" {
#endif

    void reportCompMatr1(CompMatr1 matrix);

    void reportCompMatr2(CompMatr2 matrix);

    void reportCompMatrN(CompMatrN matrix);

#ifdef __cplusplus
}
#endif



#endif // STRUCTURES_H