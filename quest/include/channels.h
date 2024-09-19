/** @file
 * Data structures for representing arbitrary channels as
 * superoperators and Kraus maps, including their constructors, 
 * getters, setters and reporters. Note the functions to
 * actually simulate these channels are exposed in decoherence.h
 * 
 * Like matrices.h, this file makes extensive use of macros to
 * overload struct initialisers for user convenience. All macros
 * herein expand to single-line definitions for safety. Some 
 * intendedly private functions are necessarily exposed here to 
 * the user, and are prefixed with an underscore.
 */

#ifndef CHANNELS_H
#define CHANNELS_H

#include "quest/include/types.h"

// C++ gets vector initialiser overloads, whereas C gets a macro
#ifdef __cplusplus
    #include <vector>
#endif



/*
 * STRUCTS
 */


typedef struct {

    // fields are const to prevent user modification
    const int numQubits;
    const qindex numRows;
    
    qcomp** cpuElems;
    qcomp* gpuElemsFlat;

    // whether the user has ever synchronised memory to the GPU, which is performed automatically
    // when calling functions like setCompMatr(), but which requires manual invocation with
    // syncCompMatr() after manual modification of the cpuElem. Note this can only indicate whether
    // the matrix has EVER been synced; it cannot be used to detect whether manual modifications
    // made after an initial sync have been re-synched. This is a heap pointer to remain mutable.
    int* const wasGpuSynced;

} SuperOp;


typedef struct {

    // fields are const to prevent user modification
    const int numQubits;

    // representation of the map as a collection of Kraus operators, kept exclusively 
    // in CPU memory, and used only for CPTP validation and reporting the map
    const qindex numMatrices;
    const qindex numRows;
    qcomp*** matrices;

    // representation of the map as a single superoperator, used for simulation
    SuperOp superop;

    // CPTP-ness is determined at validation; 0 or 1, or -1 to indicate unknown. The flag is 
    // stored in heap so even copies of structs are mutable, but pointer itself is immutable.
    int* const isCPTP;

} KrausMap;


// we define no fixed-size versions (e.g. KrausMap1/2), unlike we did for CompMatr1/2
// and DiagMatr1/2. This is because the 2-qubit superoperator is 256 elements big, and
// seems inadvisably large to be passing-by-copy through the QuEST backend layers, and
// would need explicit GPU memory allocation at each invocation of mixKrausMap2() (it
// exceeds the max number of CUDA kernel args). Furthermore, repeatedly calling
// createKrausMap2() would repeatedly invoke ~256*16 flops to compute te superoperator,
// which may be an user-astonishing overhead (more astonishing than the API asymmetry).
// Finally, computing the fixed-size superoperators must be in the header (to avoid
// the issues of qcmop interoperability, just like for getCompMatr1) and could not call
// an inner function which wouldn't be user-exposed; so we would end up redefining the
// superoperator calculation THREE times!



/*
 * BASIC FUNCTIONS
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    KrausMap createKrausMap(int numQubits, int numOperators);

    void syncKrausMap(KrausMap map);

    void destroyKrausMap(KrausMap map);

    void reportKrausMap(KrausMap map);


    SuperOp createSuperOp(int numQubits);

    void syncSuperOp(SuperOp op);

    void destroySuperOp(SuperOp op);

    void reportSuperOp(SuperOp op);

#ifdef __cplusplus
}
#endif



/*
 * POINTER INITIALISERS
 *
 * which permit users to pass heap and stack pointers in both C and C++, e.g.
 *   - qcomp** ptr = malloc(...); setSuperOp(m, ptr);
 *   - qcomp* ptrs[16]; setSuperOp(m, ptrs);
 *   - qcomp*** ptr = malloc(...); setKrausMap(m, ptr);
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif

    void setKrausMap(KrausMap map, qcomp*** matrices);

    void setSuperOp(SuperOp op, qcomp** matrix);

#ifdef __cplusplus
}
#endif



/*
 * ARRAY, VECTOR, MATRIX INITIALISERS
 *
 * which define additional overloads for arrays, VLAs, vectors and vector initialisation lists.
 * They permit C users to additionally call e.g.
 *   - qcomp arr[16][16]; setSuperOp(m, arr);
 *   - int n=16; qcomp arr[n][n]; setSuperOp(m, arr);
 *   - setKrausMap(m, (qcomp[5][16][16]) {{{...}}});
 *   - inline temporary VLA remains impossible even in C99, however
 * and C++ users gain overloads:
 *   - int n=8; std::vector vec(n); setSuperOp(m, vec);
 *   - setKrausMap(m, {{{...}}} );
 * An unintended but harmless side-effect is the exposure of functions setKrausMapFromArr(),
 * setSuperOpFromArr(), validate_setCompMatrFromArr() and validate_setSuperOpFromArr to the user.
 */


#ifdef __cplusplus

    // C++ overloads to accept vectors, which also enables vector initialiser literals

    void setKrausMap(KrausMap map, std::vector<std::vector<std::vector<qcomp>>> matrices);

    void setSuperOp(SuperOp op, std::vector<std::vector<qcomp>> matrix);
    

    // C++ cannot accept VLAs so does not define 2D array overloads

#else

    // C first defines bespoke functions which accept C99 VLAs, which we have to define here in
    // the header becauses the C++ source cannot use VLA, nor should we pass a 2D qcomp array
    // directly between C and C++ binaries (due to limited interoperability)


    // C must validate the struct fields before accessing passed 2D arrays to avoid seg-faults
    extern void _validateParamsToSetKrausMapFromArr(KrausMap map);
    extern void _validateParamsToSetSuperOpFromArr(SuperOp op);


    static inline void _setKrausMapFromArr(KrausMap map, qcomp matrices[map.numMatrices][map.numRows][map.numRows]) {
        _validateParamsToSetKrausMapFromArr(map);

        // create stack space for 2D collection of pointers, one to each input row
        qcomp* rows[map.numMatrices][map.numRows];
        qcomp** ptrs[map.numMatrices];

        // copy decayed array pointers into stack
        for (int n=0; n<map.numMatrices; n++) {
            for (qindex r=0; r<map.numRows; r++)
                rows[n][r] = matrices[n][r];
            ptrs[n] = rows[n];
        }

        setKrausMap(map, ptrs); // validation gauranteed to pass
    }

    static inline void _setSuperOpFromArr(SuperOp op, qcomp matrix[op.numRows][op.numRows]) {
        _validateParamsToSetSuperOpFromArr(op);

        // create stack space for pointers, one for each input row
        qcomp* ptrs[op.numRows];

        // copy decayed array pointers into stack
        for (qindex r=0; r<op.numRows; r++)
            ptrs[r] = matrix[r];

        setSuperOp(op, ptrs); // validation gauranteed to pass
    }


    // C then overloads setKrausMap() to call the above VLA when given arrays, using C11 Generics.
    // See the doc of getCompMatr1() in matrices.h for an explanation of Generic, and its nuances.

    #define setKrausMap(map, ...) \
        _Generic((__VA_ARGS__), \
            qcomp*** : setKrausMap, \
            default  : _setKrausMapFromArr \
        )((map), (__VA_ARGS__))

    #define setSuperOp(op, ...) \
        _Generic((__VA_ARGS__), \
            qcomp** : setSuperOp, \
            default : _setSuperOpFromArr \
        )((op), (__VA_ARGS__))

#endif



/*
 * LITERAL INITIALISERS
 *
 * which enable C users to give inline 2D and 3D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s vector overloads achieve the same thing.
 * 
 * These empower C and C++ users to call
 *   - setInlineSuperOp(m, 1, {{...}});
 *   - setInlineKrausMap(m, 2, 16, {{{...}}});
 */


#ifdef __cplusplus

    // C++ redirects to vector overloads, passing initialiser lists.  The args like 'numQb'
    // and 'numOps' are superfluous, but needed for consistency with the C API, so we additionally
    // validate that they match the struct dimensions (which requires validating the structs).

    void setInlineKrausMap(KrausMap map, int numQb, int numOps, std::vector<std::vector<std::vector<qcomp>>> matrices);

    void setInlineSuperOp(SuperOp op, int numQb, std::vector<std::vector<qcomp>> matrix);

#else

    // C defines macros which add compound literal syntax so that the user's passed lists
    // become compile-time-sized temporary arrays. C99 does not permit inline-initialised
    // VLAs, so we cannot have the macro expand to add (qcomp[matr.numRows][matr.numRows])
    // in order to preclude passing 'numQb'. We ergo accept and validate 'numQb' macro param.
    // We define private inner-functions of a macro, in lieu of writing multiline macros
    // using do-while, just to better emulate a function call for users - e.g. they
    // can wrap the macro invocations with another function call, etc.


    // C validators check 'numQb' and 'numOps' are consistent with the struct, but cannot check the user's passed literal sizes
    extern void _validateParamsToSetInlineKrausMap(KrausMap map, int numQb, int numOps);
    extern void _validateParamsToSetInlineSuperOp(SuperOp op, int numQb);


    static inline void _setInlineKrausMap(KrausMap map, int numQb, int numOps, qcomp elems[numOps][1<<numQb][1<<numQb]) {
        _validateParamsToSetInlineKrausMap(map, numQb, numOps);
        _setKrausMapFromArr(map, elems);
    }

    static inline void _setInlineSuperOp(SuperOp op, int numQb, qcomp elems[1<<(2*numQb)][1<<(2*numQb)] ) {
        _validateParamsToSetInlineSuperOp(op, numQb);
        _setSuperOpFromArr(op, elems);
    }


    #define setInlineKrausMap(map, numQb, numOps, ...) \
        _setInlineKrausMap((map), (numQb), (numOps), (qcomp[(numOps)][1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    #define setInlineSuperOp(matr, numQb, ...) \
        _setInlineSuperOp((matr), (numQb), (qcomp[1<<(2*(numQb))][1<<(2*(numQb))]) __VA_ARGS__)

#endif



/*
 * LITERAL CREATORS
 *
 * which combine creators and the inline initialisation functions, so that
 * both C and C++ users can call e.g.
 *   - SuperOp op = createInlineSuperOp(2, {{...}});
 */


#ifdef __cplusplus

    // C++ accepts vector initialiser lists

    KrausMap createInlineKrausMap(int numQubits, int numOperators, std::vector<std::vector<std::vector<qcomp>>> matrices);

    SuperOp createInlineSuperOp(int numQubits, std::vector<std::vector<qcomp>> matrix);

#else

    // C defines macros which add compound literal syntax so that the user's passed lists
    // become compile-time-sized temporary arrays. We use bespoke validation so that the
    // error messages reflect the name of the macro, rather than the inner called functions.
    // We define a private inner function per macro, in lieu of writing multiline macros
    // using do-while, just to better emulate a function call for users - e.g. they
    // can wrap the macro invocation with another function call.


    extern void _validateParamsToCreateInlineKrausMap(int numQb, int numOps);
    extern void _validateParamsToCreateInlineSuperOp(int numQb);


    static inline KrausMap _createInlineKrausMap(int numQb, int numOps, qcomp matrices[numOps][1<<numQb][1<<numQb]) {
        _validateParamsToCreateInlineKrausMap(numQb, numOps);
        KrausMap out = createKrausMap(numQb, numOps); // malloc failures will report 'createKrausMap', rather than 'inline' version. Alas!
        _setKrausMapFromArr(out, matrices);
        return out;
    }

    static inline SuperOp _createInlineSuperOp(int numQb, qcomp matrix[1<<numQb][1<<numQb]) {
        _validateParamsToCreateInlineSuperOp(numQb);
        SuperOp out = createSuperOp(numQb); // malloc failures will report 'createSuperOp', rather than 'inline' version. Alas!
        _setSuperOpFromArr(out, matrix);
        return out;
    }


    #define createInlineKrausMap(numQb, numOps, ...) \
        _createInlineKrausMap((numQb), (numOps), (qcomp[(numOps)][1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    #define createInlineSuperOp(numQb, ...) \
        _createInlineSuperOp((numQb), (qcomp[1<<(2*(numQb))][1<<(2*(numQb))]) __VA_ARGS__)

#endif



#endif // CHANNELS_H