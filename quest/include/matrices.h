/** @file
 * Definitions of all dense and diagonal matrices, their getters and setters,
 * as well as their reporting utilities. Note that Kraus maps are treated in
 * a bespoke file (krausmaps.h).
 * 
 * This file uses extensive preprocessor trickery to achieve overloaded,
 * platform agnostic, C and C++ compatible, precision agnostic, getters 
 * and setters of complex matrices. Read on to begin your adventure.
 */

#ifndef MATRICES_H
#define MATRICES_H

#include "quest/include/types.h"

// C++ gets vector initialiser overloads, whereas C gets a macro
#ifdef __cplusplus
    #include <vector>
#endif



/*
 * DENSE MATRIX STRUCTS
 *
 * which are visible to both C and C++, where qcomp resolves
 * to the native complex type. These are not de-mangled because
 * C++ structs are already C compatible. We define their fields
 * as const to prevent users mangling them.
 * 
 * The compile-time sized structs have field 'elems', while
 * dynamic-sized structs have separate 'cpuElems' and 
 * 'gpuElemsFlat', for persistent GPU allocation, and ergo need 
 * syncing. Note 'gpuElemsFlat' is always 1D (hence the name).
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

    // unitarity of the matrix (0, 1, or -1 to indicate unknown) which is lazily evaluated,
    // deferred until a function actually asserts unitarity, at which point it is computed
    // and the flag fixed until the user modifies the matrix (through sync() or setAmps() etc).
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

    // unitarity of the matrix (0, 1, or -1 to indicate unknown) which is lazily evaluated,
    // deferred until a function actually asserts unitarity, at which point it is computed
    // and the flag fixed until the user modifies the matrix (through sync() or setAmps() etc).
    // flag is stored in heap so even copies of structs are mutable, but pointer is immutable.
    // otherwise, the field of a user's struct could never be modified because of pass-by-copy.
    int* const isUnitary;

    // CPU memory; not const, so users can overwrite addresses (e.g. with NULL)
    qcomp* cpuElems;

    // GPU memory, allocated only and always in GPU-enabled QuEST environments
    qcomp* gpuElems;

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

    // unitarity of the matrix (0, 1, or -1 to indicate unknown) which is lazily evaluated,
    // deferred until a function actually asserts unitarity, at which point it is computed
    // and the flag fixed until the user modifies the matrix (through sync() or setAmps() etc).
    // flag is stored in heap so even copies of structs are mutable, but pointer is immutable.
    // otherwise, the field of a user's struct could never be modified because of pass-by-copy.
    int* const isUnitary;

    // CPU memory; not const, so users can overwrite addresses (e.g. with NULL)
    qcomp* cpuElems;

    // GPU memory, allocated only and always in GPU-enabled QuEST environments
    qcomp* gpuElems;

} FullStateDiagMatr;



/*
 * FIZED-SIZE MATRIX GETTERS VIA POINTERS
 *
 * which are defined here in the header because the 'qcomp' type is interpreted
 * distinctly by C++ (the backend) and C (user code). The C and C++ ABIs do not
 * agree on a complex type, so a qcomp (to C; a _Complex, and to C++; a std::complex)
 * cannot be directly passed between C and C++ compiled binaries; nor can a CompMatr1
 * struct which unwraps the qcomp[][] array. However, the C and C++ complex types have 
 * identical memory layouts, so pointers to qcomp types can safely be passed between
 * C and C++ binaries. Ordinarily we leverage this by defining all qcomp-handling API
 * functions in C++, and defining additional C-only wrappers in wrappers.h, which pass 
 * only pointers to qcomp. Alas, we cannot use this trick here, because the CompMatr1/2 
 * fields are declared 'const'; we cannot modify them through pointers, nor should we 
 * try to address them. Ergo we directly define these functions below (static inline to 
 * avoid symbol duplication), initializing the struct in one line. These functions will be 
 * separately interpreted by the C and C++ compilers, resolving qcomp to their individual 
 * native complex types.
 * 
 * These functions permit users to pass heap and stack pointers:
 *   - qcomp** ptr = malloc(...); getCompMatr1(ptr);
 *   - qcomp* ptrs[2]; getCompMatr1(ptrs);
 * in both C and C++. Because of 1D pointer decay, they also permit:
 *   - qcomp* ptr = malloc(...); getDiagMatr1(ptr);
 *   - qcomp arr[2]; getDiagMatr1(arr);
 */


static inline CompMatr1 getCompMatr1(qcomp** in) {

    return (CompMatr1) {
        .numQubits = 1,
        .numRows = 2,
        .elems = {
            {in[0][0], in[0][1]}, 
            {in[1][0], in[1][1]}}
    };
}

static inline CompMatr2 getCompMatr2(qcomp** in) {

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



/*
 * FIZED-SIZE MATRIX GETTERS VIA ARRAYS & VECTORS
 *
 * which define additional overloads for arrays, VLAs, C99 temporary arrays,
 * vectors and vector initialisation lists. This empowers C users to call:
 *   - qcomp arr[2][2]; getCompMatr1(arr);
 *   - int n=2; qcomp arr[n][n]; getCompMatr1(arr);
 *   - getCompMatr1( (qcomp[2][2]) {...} );
 * and C++ users call:
 *   - qcomp arr[2][2]; getCompMatr1(arr);
 *   - std::vector vec(2); getCompMatr1(vec);
 *   - getCompMatr1( {...} );
 * An unintended but harmless side-effect is the exposure of function 
 * getCompMatr1FromArr() to the user.
 */


// define the array overloads with a distinct name from the base
// C function - we will alias it with getCompMatr() using Generics

static inline CompMatr1 getCompMatr1FromArr(qcomp in[2][2]) {

    qcomp* rowPtrs[] = {in[0], in[1]};
    return getCompMatr1(rowPtrs);
}

static inline CompMatr2 getCompMatr2FromArr(qcomp in[4][4]) {

    qcomp* rowPtrs[] = {in[0], in[1], in[2], in[3]};
    return getCompMatr2(rowPtrs);
}


// no array overloads are necessary for getDiagMatr(), because
// a 1D array automatically decays to a pointer


#ifdef __cplusplus

    // C++ defines overloads which merely wrap getCompMatr1FromArr()

    static inline CompMatr1 getCompMatr1(qcomp in[2][2]) { return getCompMatr1FromArr(in); }
    static inline CompMatr2 getCompMatr2(qcomp in[4][4]) { return getCompMatr2FromArr(in); }


    // C++ also defines additional std::vector overloads (for convenience, and for inline initialisation).
    // these are defined in matrices.cpp because they invoke validation (checking vector sizes)

    CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in);
    CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in);

    DiagMatr1 getDiagMatr1(std::vector<qcomp> in);
    DiagMatr2 getDiagMatr2(std::vector<qcomp> in);

#else

    // C uses C11 Generics to effectively overload getCompMatr1/2 to accept both
    // pointers (as prior defined) and arrays (wrapping getCompMatr1FromArr()). Note:
    // - our macros below accept C99 variadic arguments so that users pass C99
    //   compound literals (e.g. (qcomp[]) {1,2}) in addition to existing ptrs.
    //   they cannot however exclude the (qcomp[]) syntax like C++ users enjoy, 
    //   which is why we will subsequently define a getInlineCompMatr1()
    // - Generics evaluate at compile-time (AFTER preprocessing) so their RHS
    //   expressions are limited; because of this, it is impossible to avoid
    //   defining the getCompMatr1FromArr() inner functions to avoid exposing them.
    // - our Generics explicitly check for pointer types (qcomp**), but we use default 
    //   to catch all array types (qcomp[][n], or qcomp(*)[] due to automatic Generic 
    //   pointer decay in GCC). This avoids us using qcomp(*)[] in the macro which
    //   would get expanded into our qcomp(re,im) macro and become invalid Generic
    //   syntax, unless we use a qcomp alias macro (which is gross). It also makes 
    //   the code more consistent with our variable-size CompMatr macros later in this 
    //   file, which cannot use VLA in Generics at all. And finally, it avoids the user
    //   having to see a Generic compilation error message when they pass an invalid
    //   type. 
    // - Generic expansion does not recurse, hence our macro safely has the same name
    //   (e.g. getCompMatr1) as the inner function, defining a true overload 
    // - we could not have _Generic's 'default' to catch unrecognised types at compile
    //   time to issue a custom message, because we must expand _Generic to a function 
    //   rather than a macro; preprocessing is finished by the time _Generic evaluates,
    //   so a macro would always be substituted before compilation and if it contained
    //   a compile-time error, it will always be triggered. A function error however
    //   would compile fine, but the error message would only be triggered at runtime
    //   when the user actually calls getCompMatr1() which is much worse than a slightly
    //   less clear compile-time error!
    
    #define getCompMatr1(...) \
        _Generic((__VA_ARGS__), \
            qcomp** : getCompMatr1, \
            default : getCompMatr1FromArr \
        )((__VA_ARGS__))

    #define getCompMatr2(...) \
        _Generic((__VA_ARGS__), \
            qcomp** : getCompMatr2, \
            default : getCompMatr2FromArr \
        )((__VA_ARGS__))

#endif



/*
 * FIXED-SIZE MATRIX GETTERS VIA LITERALS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * compound literal syntax. We expose these macros to C++ too for API consistency.
 * although C++'s getCompMatr1 vector overload achieves the same thing, and cannot
 * use C-style temporary arrays.
 * 
 * These empower C and C++ users to call
 *   - getInlineCompMatr1( {{1,2},{3,4}} )
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
 * VARIABLE-SIZE MATRIX SETTERS VIA POINTERS
 *
 * These functions permit users to pass heap and stack pointers:
 *   - qcomp** ptr = malloc(...); setCompMatr(m, ptr);
 *   - qcomp* ptrs[8]; setCompMatr(m, ptrs);
 * in both C and C++. By decay, they also permit arrays to diagonals:
 *   - qcomp* ptr = malloc(...); setDiagMatr(m, ptr);
 *   - qcomp arr[8]; setDiagMatr(m, arr); 
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    void setCompMatr(CompMatr matr, qcomp** vals);

    void setDiagMatr(DiagMatr out, qcomp* in);

    void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems);

#ifdef __cplusplus
}
#endif



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA ARRAYS & VECTORS
 *
 * which define additional overloads for arrays, VLAs, vectors and vector initialisation lists.
 * C users can call:
 *   - qcomp arr[8][8]; setCompMatr(m, arr);
 *   - int n=8; qcomp arr[n][n]; setCompMatr(m, arr);
 * and C++ users can call:
 *   - int n=8; std::vector vec(n); setCompMatr(vec);
 *   - setCompMatr( {...} );
 * An unintended but harmless side-effect is the exposure of functions setCompMatrFromArr() and 
 * validate_setCompMatrFromArr() to the user.
 */


#ifdef __cplusplus

    // C++ defines vector overloads, permitting inline initialisation

    void setCompMatr(CompMatr out, std::vector<std::vector<qcomp>> in);

    void setDiagMatr(DiagMatr out, std::vector<qcomp> in);

    void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, std::vector<qcomp> in);


    // C++ cannot accept 2D arrays at all, because it does not support C99 VLA. 
    // It can however accept 1D arrays (which decay to pointers) already to setDiagMatr()

#else

    // C first defines a bespoke functions receiving C99 VLAs, which we have to define here in
    // the header becauses the C++ source cannot use VLA, nor should we pass a 2D qcomp array
    // directly between C and C++ binaries (due to limited interoperability)

    extern void validate_matrixFields(CompMatr matr, const char* caller);

     // static inline to avoid header-symbol duplication
    static inline void setCompMatrFromArr(CompMatr matr, qcomp arr[matr.numRows][matr.numRows]) {

        // this function will allocate stack memory of size matr.numRows, but that field could
        // be invalid since matr hasn't been validated, so we must first invoke validation. Note
        // the caller will likely have called setCompMatr() or setInlineCompMatr(); we'll just
        // report the former for relative clarity.
        validate_matrixFields(matr, "setCompMatr");

        // new ptrs array safely fits in stack, since it's sqrt-smaller than user's passed stack array
        qcomp* ptrs[matr.numRows];

        // collect pointers to each row of arr
        for (qindex r=0; r<matr.numRows; r++)
            ptrs[r] = arr[r];

        // array decays to qcomp**, and *FromPtr function re-performs validation (eh)
        setCompMatr(matr, ptrs);
    }


    // C then overloads setCompMatr() to call the above VLA when given arrays, using C11 Generics.
    // See the doc of getCompMatr1() above for an explanation of Generic, and its nuances

    #define setCompMatr(matr, ...) \
        _Generic((__VA_ARGS__), \
            qcomp** : setCompMatr, \
            default : setCompMatrFromArr \
        )((matr), (__VA_ARGS__))


    // no need to define bespoke overload for diagonal matrices, because 1D arrays decay to pointers

#endif



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA LITERALS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s vector overloads achieve the same thing.
 * 
 * These empower C and C++ users to call
 *   - setCompMatr(m, {{1,2},{3,4}} )
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