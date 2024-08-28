/** @file
 * Signatures of API matrix data structures, and their getters and setters, 
 * as well as reporting utilities.
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
 * We employ similar tricks to make setCompMatr(), but define setCompMatrFromArr()
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

#ifndef MATRICES_H
#define MATRICES_H

#include "quest/include/types.h"

// C++ gets vector initialiser overloads, whereas C gets a macro
#ifdef __cplusplus
    #include <vector>
#endif



/*
 * COMPLEX MATRIX STRUCTS
 *
 * which are visible to both C and C++, where qcomp resolves
 * to the native complex type. These are not de-mangled because
 * C++ structs are already C compatible. We define their fields
 * as const to prevent users mangling them.
 * 
 * The compile-time sized structs have field 'elems', while
 * dynamic-sized structs have separate 'cpuElems' and 'gpuElems',
 * for persistent GPU allocation, and ergo need syncing. Note
 * 'gpuElems' is always 1D.
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

    // unitarity determined at sync; 0 or 1, or -1 to indicate unknown (when validation disabled).
    // flag is stored in heap so even copies of structs are mutable, but pointer is immutable.
    // otherwise, the field of a user's struct could never be modified because of pass-by-copy.
    int* const isUnitary;

    // 2D CPU memory; not const, so users can overwrite addresses (e.g. with NULL)
    qcomp** cpuElems;

    // row-flattened elems in GPU memory, allocated only
    // and always in GPU-enabled QuEST environments
    qcomp* gpuElemsFlat;

} CompMatr;



/*
 * DIAGONAL MATRIX STRUCTS
 *
 * with all the same nuances as the CompMatr structs described above. 
 */


typedef struct {

    // const to prevent user modification
    const int numQubits;
    const qindex numElems;

    // elems are not const so that users can modify them after initialisation
    qcomp elems[2];

} DiagMatr1;


typedef struct {

    // const to prevent user modification
    const int numQubits;
    const qindex numElems;

    // elems are not const so that users can modify them after initialisation
    qcomp elems[4];

} DiagMatr2;


typedef struct {

    // const to prevent user modification
    const int numQubits;
    const qindex numElems;

    // unitarity determined at sync; 0 or 1, or -1 to indicate unknown (when validation disabled).
    // flag is stored in heap so even copies of structs are mutable, but pointer is immutable
    int* const isUnitary;

    // CPU memory; not const, so users can overwrite addresses (e.g. with NULL)
    qcomp* cpuElems;

    // GPU memory, allocated only and always in GPU-enabled QuEST environments
    qcomp* gpuElemsFlat;

} DiagMatr;



/*
 * DISTRIBUTED MATRIX STRUCTS
 */


typedef struct {

    // data deployment configuration
    const int isDistributed;

    // const to prevent user modification
    const int numQubits;
    const qindex numElems;

    // will equal numElems if distribution is disabled at runtime (e.g. via autodeployment)
    const qindex numElemsPerNode;

    // unitarity determined at sync; 0 or 1, or -1 to indicate unknown (when validation disabled).
    // flag is stored in heap so even copies of structs are mutable, but pointer is immutable
    int* const isUnitary;

    // CPU memory; not const, so users can overwrite addresses (e.g. with NULL)
    qcomp* cpuElems;

    // GPU memory, allocated only and always in GPU-enabled QuEST environments
    qcomp* gpuElemsFlat;

} FullStateDiagMatr;



/*
 * EXPLICIT FIZED-SIZE DENSE MATRIX INITIALISERS
 *
 * which are defined here in the header because the 'qcomp' type is interpreted
 * distinctly by C++ (the backend) and C (user code). The C and C++ ABIs do not
 * agree on a complex type, so a qcomp (to C; a _Complex, and to C++; a std::complex)
 * cannot be directly passed between C and C++ compiled binaries; nor can a CompMatr1
 * struct which unwraps the qcomp[][] array. However, the C and C++ complex types have 
 * identical memory layouts, so pointers to qcomp types can safely be passed between
 * C and C++ binaries. 
 * 
 * Ordinarily we leverage this by defining all qcomp-handling API
 * functions in C++, and defining additional C-only wrappers in wrappers.h, which pass 
 * only pointers to qcomp. Alas, we cannot use this trick here, because the CompMatr1/2 
 * fields are declared 'const'; we cannot modify them through pointers, nor should we 
 * try to address them. Ergo we directly define these functions below (static inline to 
 * avoid symbol duplication), initializing the struct in one line. These functions will be 
 * separately interpreted by the C and C++ compilers, resolving qcomp to their individual 
 * native complex types.
 * 
 * Note separate pointer and array definitions are not necessary for DiagMatr, since
 * it 1D array field decays to a pointer. Finally, we mention these methods need never
 * be directly called by the user, who should instead call the overloaded getCompMatr1()
 * which will automatically dispatch to the below methods - even in C, using generics!
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
 * which permit both C and C++ users to call getCompMatr1/2() and pass
 * arrays or pointers, without having to call the above specialised
 * functions. We are effectively using macros to extend C++'s
 * overloaded API to C, though C++ users can additionally pass vectors.
 * 
 * The situation is simpler for DiagMatr since 1D array references
 * decay to pointers. We still wish to enable C++ users to pass
 * vectors and their in-place initialiser lists, so we overload those.
 */


#ifdef __cplusplus

    // C++ uses overloads, accepting even vector initialiser lists,
    // which are defined in matrices.cpp.

    CompMatr1 getCompMatr1(qcomp in[2][2]);
    CompMatr1 getCompMatr1(qcomp** in);
    CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in);

    CompMatr2 getCompMatr2(qcomp in[4][4]);
    CompMatr2 getCompMatr2(qcomp** in);
    CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in);

#else

    // C uses a header macro with C11 generic inspection to expand the macro
    // (at compile-time, rather than during pre-processing) with one of
    // the above explicit definitions. Note:
    // - we cannot accept C++ vectors (duh) so direct {{...}} initialisation
    //   isn't possible; users have to use C99 compound literals instead, or
    //   our subsequent definition of getInlineCompMatr().
    // - we explicitly check for pointers (qcomp**), but we use default to catch
    //   all array types (qcomp[][n], and qcomp(*)[] due to Generic pointer decay
    //   in GCC). This is to avoid dealing with VLA (illegal) in the variable
    //   size matrices, to avoid having to alias qcomp (because qcomp(*)[] 
    //   out qcomp(re,im) macro), and to spare the user from the _Generic
    //   compilation error (a bit hard to read) when passing an incorrect type.
    // - we use __VA_ARGS__ to accept both prefined variables and multiple-token 
    //   compound literals, i.e. C99 temporary arrays of syntax (qcomp[][]) {{...}}
    // - we cannot use _Generic's 'default' to catch unrecognised types at compile
    //   time (although that would be lovely) because we must expand _Generic to a 
    //   function (not a macro; preprocessing is finished by the time _Generic
    //   evaluates). Such a function could only throw an error at runtime.

    #define getCompMatr1(...) \
        _Generic((__VA_ARGS__), \
            qcomp** : getCompMatr1FromPtr, \
            default : getCompMatr1FromArr \
        )((__VA_ARGS__))

    #define getCompMatr2(...) \
        _Generic((__VA_ARGS__), \
            qcomp** : getCompMatr2FromPtr, \
            default : getCompMatr2FromArr \
        )((__VA_ARGS__))

#endif


// we define getDiagMatr1/2 in this header for the same reasons that we did so for 
// getCompMatr1/2FromArr/Ptr(). That is, we cannot pass DiagMatr1/2 instances between
// C and C++ binaries, so we define these simple inline functions afresh for each.
// Conveniently, we don't need to overload or use Generics, since arrays decay to ptrs.

static inline DiagMatr1 getDiagMatr1(qcomp* in) {

    return (DiagMatr1) {
        .numQubits = 1,
        .numElems = 2,
        .elems = {in[0], in[1]}
    };
}

static inline DiagMatr2 getDiagMatr2(qcomp* in) {

    return (DiagMatr2) {
        .numQubits = 2,
        .numElems = 4,
        .elems = {in[0], in[1], in[2], in[3]}
    };
}

// C++ users can additionally initialise from vectors, enabling in-line initialisation

#ifdef __cplusplus

    DiagMatr1 getDiagMatr1(std::vector<qcomp> in);

    DiagMatr2 getDiagMatr2(std::vector<qcomp> in);

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


    #define getInlineDiagMatr1(...) \
        getDiagMatr1(__VA_ARGS__)

    #define getInlineDiagMatr2(...) \
        getDiagMatr2(__VA_ARGS__)

#else

    // C adds compound literal syntax to make a temporary array

    #define getInlineCompMatr1(...) \
        getCompMatr1FromArr((qcomp[2][2]) __VA_ARGS__)

    #define getInlineCompMatr2(...) \
        getCompMatr2FromArr((qcomp[4][4]) __VA_ARGS__)


    #define getInlineDiagMatr1(...) \
        getDiagMatr1((qcomp[2]) __VA_ARGS__)

    #define getInlineDiagMatr2(...) \
        getDiagMatr2((qcomp[4]) __VA_ARGS__)

#endif



/*
 * VARIABLE-SIZE MATRIX CONSTRUCTORS
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    CompMatr createCompMatr(int numQubits);

    DiagMatr createDiagMatr(int numQubits);

    FullStateDiagMatr createFullStateDiagMatr(int numQubits);

    FullStateDiagMatr createCustomFullStateDiagMatr(int numQubits, int useDistrib);


    void destroyCompMatr(CompMatr matrix);

    void destroyDiagMatr(DiagMatr matrix);

    void destroyFullStateDiagMatr(FullStateDiagMatr matrix);


    void syncCompMatr(CompMatr matr);

    void syncDiagMatr(DiagMatr matr);

    void syncFullStateDiagMatr(FullStateDiagMatr matr);

#ifdef __cplusplus
}
#endif



/*
 * EXPLICIT VARIABLE-SIZE MATRIX INITIALISERS
 *
 * not necessary for DiagMatr since it maintains a 1D array,
 * which automatically decays to a qcomp pointer
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    void setCompMatrFromPtr(CompMatr matr, qcomp** vals);

#ifdef __cplusplus
}
#endif


// permit only C compilers to have a VLA version (not supported nor needed by C++)
#ifndef __cplusplus

    // expose this function's bespoke validation
    extern void validate_setCompMatrFromArr(CompMatr out);

     // static inline to avoid header-symbol duplication
    static inline void setCompMatrFromArr(CompMatr matr, qcomp arr[matr.numRows][matr.numRows]) {

        // this function will allocate stack memory of size matr.numRows, but that field could
        // be invalid since matr hasn't been validated, so we must invoke bespoke validation
        validate_setCompMatrFromArr(matr);

        // new ptrs array safely fits in stack, since it's sqrt-smaller than user's passed stack array
        qcomp* ptrs[matr.numRows];

        // collect pointers to each row of arr
        for (qindex r=0; r<matr.numRows; r++)
            ptrs[r] = arr[r];

        // array decays to qcomp**, and *FromPtr function re-performs validation (eh)
        setCompMatrFromPtr(matr, ptrs);
    }

#endif



/*
 * OVERLOADED VARIABLE-SIZE MATRIX INITIALISERS
 */


// both C and C++ can safely pass pointers or arrays to diagonals (since arrays decay)

#ifdef __cplusplus
extern "C" {
#endif

    void setDiagMatr(DiagMatr out, qcomp* in);

    void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems);

#ifdef __cplusplus
}
#endif


// but setCompMatr() is overloaded because pointers and 2D arrays are distinct

#ifdef __cplusplus

    // C++ uses overloads, accepting even vector initialiser lists, but cannot ever accept 2D arrays

    void setCompMatr(CompMatr out, qcomp** in);

    void setCompMatr(CompMatr out, std::vector<std::vector<qcomp>> in);


    // we also give diagonal matrices vector overloads, to allow initialiser lists

    void setDiagMatr(DiagMatr out, std::vector<qcomp> in);

    void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, std::vector<qcomp> in);

#else

    // C uses C11 generics, resolved at compile-time (after preprocessing)

    #define setCompMatr(matr, ...) \
        _Generic((__VA_ARGS__), \
            qcomp** : setCompMatrFromPtr, \
            default : setCompMatrFromArr \
        )((matr), (__VA_ARGS__))

#endif



/*
 * LITERAL VARIABLE-SIZE MATRIX INITIALISERS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s vector overloads achieve the same thing
 */


#ifdef __cplusplus

    // C++ gets an explicit redirect to set*Matr(std::vector...), ignoring numQb and numElems (blegh)

    #define setInlineCompMatr(matr, numQb, ...) \
        setCompMatr(matr, __VA_ARGS__)

    #define setInlineDiagMatr(matr, numQb, ...) \
        setDiagMatr(matr, __VA_ARGS__)

    #define setInlineFullStateDiagMatr(matr, startInd, numElems, ...) \
        setFullStateDiagMatr(matr, startInd, __VA_ARGS__)

#else 

    // C creates a compile-time-sized temporary array via a compound literal. We sadly
    // cannot use (qcomp[matr.numRows][matr.numRows]) to preclude passing 'numQb' 
    // because VLAs cannot be initialised inline.

    #define setInlineCompMatr(matr, numQb, ...) \
        setCompMatrFromArr(matr, (qcomp[1<<numQb][1<<numQb]) __VA_ARGS__)


    // 1D array arguments fortunately decay to pointers 

    #define setInlineDiagMatr(matr, numQb, ...) \
        setDiagMatr(matr, (qcomp[1<<numQb]) __VA_ARGS__)

    #define setInlineFullStateDiagMatr(matr, startInd, numElems, ...) \
        setFullStateDiagMatr(matr, startInd, (qcomp[numElems]) __VA_ARGS__, numElems)

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

    void reportCompMatr(CompMatr matrix);


    void reportDiagMatr1(DiagMatr1 matrix);

    void reportDiagMatr2(DiagMatr2 matrix);

    void reportDiagMatr(DiagMatr matrix);


    void reportFullStateDiagMatr(FullStateDiagMatr matr);

#ifdef __cplusplus
}
#endif



#endif // MATRICES_H