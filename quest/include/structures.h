/** @file
 * Signatures of API data structures like gate matrices. 
 * Note QuESTEnv and Qureg structs have their own signatures 
 * in environment.h and qureg.h respectively.
 */

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "quest/include/types.h"



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
 * MATRIX CONSTRUCTORS
 *
 * some of which are exposed only to C++ because they are incompatile with C binaries
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

    CompMatr1 getCompMatr1(qcomp in[2][2]);

    CompMatr2 getCompMatr2(qcomp in[4][4]);

#endif



/*
 * VARIABLE-SIZE MATRIX CONSTRUCTORS
 */

// de-mangle so below are directly callable by C binary
#ifdef __cplusplus
extern "C" {
#endif

    CompMatrN createCompMatrN(int numQubits);

    void destroyCompMatrN(CompMatrN matrix);

    void syncCompMatrN(CompMatrN matr);

    void setCompMatrN(CompMatrN matr, qcomp** vals);

#ifdef __cplusplus
}
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