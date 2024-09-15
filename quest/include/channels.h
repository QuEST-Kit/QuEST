/** @file
 * Data structures for representing arbitrary channels as
 * superoperators and Kraus maps, including their constructors, 
 * getters, setters and reporters. Note the functions to
 * actually simulate these channels are exposed in decoherence.h
 * 
 * Like matrices.h, this file makes extensive use of macros to
 * overload struct initialisers for user convenience. All macros
 * herein expand to single-line definitions for safety.
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

    // we must regrettably expose bespoke validation
    extern void validate_setKrausMapFromArr(KrausMap map);
    extern void validate_setSuperOpFromArr(SuperOp op);

    static inline void setKrausMapFromArr(KrausMap map, qcomp matrices[map.numMatrices][map.numRows][map.numRows]) {

        // we validate map's fields before the below stack allocation, since the fields could be 
        // invalid (like when the user declared struct without calling createKrausMap()), which
        // would otherwise trigger a seg-fault. It is fine that 'matrices' is already declared as 
        // a potentially invalid dimension (e.g. matrices[-1][-1][-1]) as long as we do not access 
        // any elements before validation
        validate_setKrausMapFromArr(map);

        // create stack space for 2D collection of pointers, one to each input row
        qcomp* rows[map.numMatrices][map.numRows];
        qcomp** ptrs[map.numMatrices];

        // copy decayed array pointers into stack
        for (int n=0; n<map.numMatrices; n++) {
            for (qindex r=0; r<map.numRows; r++)
                rows[n][r] = matrices[n][r];
            ptrs[n] = rows[n];
        }

        setKrausMap(map, ptrs);
    }

    static inline void setSuperOpFromArr(SuperOp op, qcomp matrix[op.numRows][op.numRows]) {

        // we validate op's fields before the below stack allocation, since the fields could be 
        // invalid (e.g. if op was declared but not initialised). The passed 'matrix' arrays
        // are invalid in that case, but it's fine if we do not access them pre-validation
        validate_setSuperOpFromArr(op);

        // create stack space for pointers, one for each input row
        qcomp* ptrs[op.numRows];

        // copy decayed array pointers into stack
        for (qindex r=0; r<op.numRows; r++)
            ptrs[r] = matrix[r];

        setSuperOp(op, ptrs);
    }


    // C then overloads setKrausMap() to call the above VLA when given arrays, using C11 Generics.
    // See the doc of getCompMatr1() in matrices.h for an explanation of Generic, and its nuances.

    #define setKrausMap(map, ...) \
        _Generic((__VA_ARGS__), \
            qcomp*** : setKrausMap, \
            default : setKrausMapFromArr \
        )((map), (__VA_ARGS__))

    #define setSuperOp(op, ...) \
        _Generic((__VA_ARGS__), \
            qcomp** : setSuperOp, \
            default : setSuperOpFromArr \
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

    // C++ uses the vector overloads, ignoring numOps and numQb parameters (blegh)

    #define setInlineKrausMap(map, numQb, numOps, ...) \
        setKrausMap((map), __VA_ARGS__);

    #define setInlineSuperOp(map, numQb, ...) \
        setSuperOp((map), __VA_ARGS__);

#else

    // C creates a compile-time-sized temporary array via a compound literal
    // (we sadly cannot use inline VLAs to preclude passing numQb and numOp)

    #define setInlineKrausMap(map, numQb, numOps, ...) \
        setKrausMapFromArr((map), (qcomp[numOps][1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    #define setInlineSuperOp(map, numQb, ...) \
        setSuperOpFromArr((map), (qcomp[1<<(2*(numQb))][1<<(2*(numQb))]) __VA_ARGS__)


    // TODO: currently, the macro parameters above are used only to create the 
    // temporary arrays given to setKrausMapFromArr() and are then discarded,
    // but in principle we could additionally validate that they match the
    // fields of argument map (a KrausMap). Mismatch is a reasonable error;
    // the user simply has to pass different numbers to create() and set(),
    // as they might do from copying code, and without validation, will
    // encounter an uninterpretable seg-fault. Admittedly, I did it myself!

#endif



#endif // CHANNELS_H