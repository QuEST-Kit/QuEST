/** @file
 * Signatures of API data structures like gate matrices. Note QuESTEnv and Qureg 
 * structs have their own signatures in environment.h and qureg.h respectively.
 * 
 * This file uses extensive preprocessor trickery to achieve platform agnostic,
 * C and C++ compatible, type agnostic, getters and setters of complex matrices.
 * First, we define C and C++ "explicit" functions like getCompMatr1FromArr();
 * these definitions are in this header because of qcomp interoperability issues
 * (see below). The function names explicitly disambiguate the input types.
 * Next, we define getCompMatr1() as a generic C macro and C++ overloads, so
 * that users can agnostically pass 2D pointers, arrays of pointers, or 2D arrays,
 * including compound literals like (qcomp[2][2]) {{...}}. So far, only the C++ 
 * overloads accept inline literals like {{...}} via std::vector; C requires
 * using the ugly and verbose compound literal syntax to create temporary arrays.
 * So, we next define getInlineCompMatr1() as a C and C++ (for consistency) macro 
 * which accepts direct, concise literals like getInlineCompMatr1({{...}}). Viola! 
 * 
 * We employ similar tricks to make setCompMatrN(), but define setCompMatrNFromArr()
 * in this header (exposed only to C) because it requires C++-incompatible VLAs;
 * this definition must invoke a bespoke validation function (gross).
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


typedef struct {

    // const to prevent user modification
    const int numQubits;
    const qindex numRows;

    // elems are not const so that users can modify them after initialisation
    qcomp elems[2][2];

} CompMatr1;


typedef struct {

    // const to prevent user modification
    const int numQubits;
    const qindex numRows;

    // elems are not const so that users can modify them after initialisation
    qcomp elems[4][4];

} CompMatr2;


typedef struct {

    // const to prevent user modification
    const int numQubits;
    const qindex numRows;

    // consts prevent overwriting mem address (of outer list, and each
    // row therein) but permit modifying each element
    qcomp* const* const elems;

    // row-flattened elems in GPU memory, allocated only
    // in GPU-enabled QuEST environments (regardless of Quregs)
    qcomp* const gpuElems;

} CompMatrN;



/*
 * EXPLICIT FIXED-SIZE MATRIX INITIALISERS
 *
 * which are defined here in the header because the 'qcomp' type is interpreted
 * distinctly by C++ (the backend) and C (user code). The C and C++ ABIs do not
 * agree on a complex type, so a qcomp (to C; a _Complex, and to C++; a std::complex)
 * cannot be directly passed between C and C++ compiled binaries; nor can a CompMatr1
 * struct which unwraps the qcomp[][] array. However, the C and C++ complex types have 
 * identical memory layouts, so pointers to qcomp types can safely be passed between
 * C and C++ binaries. Ordinarily we leverage this by defining all qcomp-handling API
 * functions in C++, and defining additional C-only wrappers in wrappers.h, which
 * pass only pointers to qcomp.
 * 
 * Alas, we cannot use this trick for CompMatr, because they are declared 'const';
 * we cannot modify them through pointers, nor should we try to address them. Ergo
 * we directly define these functions below (static inline to avoid symbol duplication),
 * initializing the const CompMatr in one line. The functions are separately interpreted 
 * by the C and C++ compilers, resolving to their individual native types.
 */


static inline CompMatr1 getCompMatr1FromArr(qcomp in[2][2]) {

    return (CompMatr1) {
        .numQubits = 1,
        .numRows = 2,
        .elems = {
            {in[0][0], in[0][1]}, 
            {in[1][0], in[1][1]}}
    };
}

static inline CompMatr1 getCompMatr1FromPtr(qcomp** in) {

    return (CompMatr1) {
        .numQubits = 1,
        .numRows = 2,
        .elems = {
            {in[0][0], in[0][1]}, 
            {in[1][0], in[1][1]}}
    };
}


static inline CompMatr2 getCompMatr2FromArr(qcomp in[4][4]) {

    return (CompMatr2) {
        .numQubits = 2,
        .numRows = 4,
        .elems = {
            {in[0][0], in[0][1], in[0][2], in[0][3]},
            {in[1][0], in[1][1], in[1][2], in[1][3]},
            {in[2][0], in[2][1], in[2][2], in[2][3]},
            {in[3][0], in[3][1], in[3][2], in[3][3]}}
    };
}

static inline CompMatr2 getCompMatr2FromPtr(qcomp** in) {

    return (CompMatr2) {
        .numQubits = 2,
        .numRows = 4,
        .elems = {
            {in[0][0], in[0][1], in[0][2], in[0][3]},
            {in[1][0], in[1][1], in[1][2], in[1][3]},
            {in[2][0], in[2][1], in[2][2], in[2][3]},
            {in[3][0], in[3][1], in[3][2], in[3][3]}}
    };
}



/*
 * OVERLOADED FIXED-SIZE MATRIX INITIALISERS
 *
 * which permit both C and C++ users to call getCompMatr1() and pass
 * arrays or pointers, without having to call the above specialised
 * functions. We are effectively using macros to extend C++'s
 * overloaded API to C, though C++ users can additionally pass vectors.
 */


#ifdef __cplusplus

    // C++ uses overloads, accepting even vector initialiser lists,
    // which are defined in structures.cpp.

    CompMatr1 getCompMatr1(qcomp in[2][2]);
    CompMatr1 getCompMatr1(qcomp** in);
    CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in);

    CompMatr2 getCompMatr2(qcomp in[4][4]);
    CompMatr2 getCompMatr2(qcomp** in);
    CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in);

#else

    // C uses a header macro with C11 compile-time type inspection to expand
    // (at compile-time, rather than during pre-processing) with one of
    // the above inlined definitions. Note:
    // - we cannot accept C++ vectors (duh) so direct {{...}} initialisation
    //   isn't possible; users have to use C99 compound literals instead,
    //   which we address with a subsequent definition of getInlineCompMatr().
    // - the _Generic does not require a map for the qcomp[2][2] type, because
    //   a passed qcomp[2][2] decays to qcomp(*)[2], indistinguishable from
    //   an array of pointers which is already in the map.
    // - we use __VA_ARGS__ to accept multiple-token compound literals,
    //   i.e. C99 temporary arrays of syntax (qcomp[][]) {{...}}
    // - we cannot use _Generic's 'default' to catch unrecognised types at compile
    //   time (although that would be lovely) because we must expand _Generic to a 
    //   function (not a macro; preprocessing is finished by the time _Generic
    //   evaluates). Such a function could only throw an error at runtime, which is
    //   worse than letting _Generic throw an error at compile-time.

    // Another problem: type qcomp(*) below would erroneously invoke the qcomp(re,im) macro.
    // Preventing expansion using (qcomp)(*) leads to _Generic not recognising the type.
    // So, in desperation, we alias the type qcomp with a name that has no colliding macro.
    typedef qcomp qalias;

    #define getCompMatr1(...) \
        _Generic((__VA_ARGS__), \
            qalias(*)[2] : getCompMatr1FromArr, \
            qcomp**      : getCompMatr1FromPtr  \
        )((__VA_ARGS__))

    #define getCompMatr2(...) \
        _Generic((__VA_ARGS__), \
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

    #define getInlineCompMatr1(...) \
        getCompMatr1(__VA_ARGS__)

    #define getInlineCompMatr2(...) \
        getCompMatr2(__VA_ARGS__)

#else

    // C adds compound literal syntax to make a temporary array

    #define getInlineCompMatr1(...) \
        getCompMatr1FromArr((qcomp[2][2]) __VA_ARGS__)

    #define getInlineCompMatr2(...) \
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

    // C uses a header macro with C11 compile-time type inspection to expand
    // (at compile-time, rather than during pre-processing) with one of
    // the above inlined definitions. Note:
    // - we cannot accept C++ vectors (duh) so direct {{...}} initialisation
    //   isn't possible; users have to use C99 compound literals instead,
    //   which we address with a subsequent definition of setInlineCompMatrN().
    // - the _Generic does not require a map for the qcomp[][] type, because
    //   a passed qcomp[][] decays to qcomp(*)[], indistinguishable from
    //   an array of pointers which is already in the map.
    // - we use __VA_ARGS__ to accept multiple-token compound literals,
    //   i.e. C99 temporary arrays of syntax (qcomp[][]) {{...}}
    // - we cannot use _Generic's 'default' to catch unrecognised types at compile
    //   time (although that would be lovely) because we must expand _Generic to a 
    //   function (not a macro; preprocessing is finished by the time _Generic
    //   evaluates). Such a function could only throw an error at runtime, which is
    //   worse than letting _Generic throw an error at compile-time.

    // Another problem: type qcomp(*) below would erroneously invoke the qcomp(re,im) macro,
    // so we re-use 'qalias = qcomp' as defined at getCompMatr1.

    #define setCompMatrN(matr, ...) \
        _Generic((__VA_ARGS__),   \
            qalias(*)[] : setCompMatrNFromArr, \
            qcomp**     : setCompMatrNFromPtr  \
        )((matr), (__VA_ARGS__))

#endif



/*
 * LITERAL VARIABLE-SIZE MATRIX INITIALISERS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s getCompMatr1 vector overload achieves the same thing
 */


#ifdef __cplusplus

    // C++ gets an explicit redirect to setCompMatrN(std::vector...), ignoring numQb (blegh)

    #define setInlineCompMatrN(matr, numQb, ...) \
        setCompMatrN(matr, __VA_ARGS__)

#else 

    // C creates a compile-time-sized temporary array via a compound literal. We sadly
    // cannot use (qcomp[matr.numRows][matr.numRows]) to preclude passing 'numQb' 
    // because VLAs cannot be initialised inline.

    #define setInlineCompMatrN(matr, numQb, ...) \
        setCompMatrNFromArr(matr, (qcomp[1<<numQb][1<<numQb]) __VA_ARGS__)

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