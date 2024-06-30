/** @file
 * Definitions of API data structures like gate matrices. 
 * Note QuESTEnv and Qureg structs have their own definitions 
 * in environment.cpp and qureg.cpp respectively.
 * 
 * This file defines many "layers" of initialisation of complex
 * matrices, as explained in the header file.
 */

#include "quest/include/structures.h"
#include "quest/include/environment.h"
#include "quest/include/types.h"

#include "quest/src/comm/comm_config.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>


// for concise reporters
using namespace form_substrings;



/*
 * PRIVATE UTILITES 
 */

// A and B can both be qcomp** or qcomp[][] (mixed)
template<typename A, typename B> 
void populateCompMatrElems(A out, B in, qindex dim) {

    for (qindex r=0; r<dim; r++)
        for (qindex c=0; c<dim; c++)
            out[r][c] = in[r][c];
}

// T can be CompMatr1/2, B can be qcomp** or qcomp[][]
template<class T, typename B>
T getCompMatrFromElems(B in, int num) {

    T out;
    out.numQubits = num;
    out.numRows = powerOf2(num);
    populateCompMatrElems(out.elems, in, out.numRows);
    return out;
}



/*
 * EXPLICIT FIXED-SIZE MATRIX INITIALISERS
 * 
 * and their corresponding C-compatible alternatives which are wrapped by wrappers.h.
 * The wrappers are necessary because of the non-interchangeable C and C++ qcomp types.
 * See structures.h for an explanation. These definitions will always exist (even when
 * the user compiles their own code in C), but will be name-mangled.
 * 
 * The separate *FromArr an *FromPtr definitions are so that we can safe and simply
 * wrap these in a C-compatible param-agnostic macro, and C++-compatible overloads,
 * so that users can avoid having to use compound literals and pass {{...}} in-place.
 */

CompMatr1 getCompMatr1FromArr(qcomp in[2][2]) {
    return getCompMatrFromElems<CompMatr1>(in, 1);
}
CompMatr1 getCompMatr1FromPtr(qcomp** in) {
    return getCompMatrFromElems<CompMatr1>(in, 1);
}
CompMatr2 getCompMatr2FromArr(qcomp in[4][4]) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2FromPtr(qcomp** in) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}

extern "C" void wrap_getCompMatr1FromArr(CompMatr1* out, qcomp in[2][2]) {
    *out = getCompMatr1FromArr(in);
}
extern "C" void wrap_getCompMatr1FromPtr(CompMatr1* out, qcomp** in) {
    *out = getCompMatr1FromPtr(in);
}
extern "C" void wrap_getCompMatr2FromArr(CompMatr2* out, qcomp in[4][4]) {
    *out = getCompMatr2FromArr(in);
}
extern "C" void wrap_getCompMatr2FromPtr(CompMatr2* out, qcomp** in) {
    *out = getCompMatr2FromPtr(in);
}



/*
 * OVERLOADED FIXED-SIZE MATRIX INITIALISERS
 *
 * Only exposed to C++; equivalent C macros are defined in the header.
 * Only the C++ vector overloads can validate their literal dimensions.
 */

CompMatr1 getCompMatr1(qcomp in[2][2]) {
    return getCompMatrFromElems<CompMatr1>(in, 1);
}
CompMatr1 getCompMatr1(qcomp** in) {
    return getCompMatrFromElems<CompMatr1>(in, 1);
}
CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in) {
    validate_numMatrixElems(1, in, __func__);

    return getCompMatrFromElems<CompMatr1>(in, 1);
}

CompMatr2 getCompMatr2(qcomp in[4][4]) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2(qcomp** in) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in) {
    validate_numMatrixElems(2, in, __func__);

    return getCompMatrFromElems<CompMatr2>(in, 2);
}



/*
 * VARIABLE SIZE MATRIX CONSTRUCTORS
 *
 * all of which are de-mangled for both C++ and C compatibility
 */


extern "C" CompMatrN createCompMatrN(int numQubits) {
    validate_envInit(__func__);
    validate_newMatrixNumQubits(numQubits, __func__);

    // allocate GPU memory if the env is GPU-accelerated, regardless of (unknown) Qureg allocations
    bool isGpuAccel = getQuESTEnv().isGpuAccelerated;

    // validation ensures these (and below mem sizes) never overflow
    qindex numRows = powerOf2(numQubits);
    qindex numElems = numRows * numRows;

    CompMatrN out = {
        .numQubits = numQubits,
        .numRows = numRows,

        // 2D CPU memory (NULL if failed)
        .elems = (qcomp**) malloc(numRows * sizeof *out.elems),

        // 1D GPU memory (NULL if failed or not needed)
        .gpuElems = (isGpuAccel)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    // only if outer CPU allocation succeeded, attempt to allocate each row array
    if (out.elems != NULL)
        for (qindex r=0; r < numRows; r++)
            out.elems[r] = (qcomp*) calloc(numRows, sizeof **out.elems); // NULL if failed

    // check all CPU & GPU malloc and calloc's succeeded (no attempted freeing if not)
    bool isNewMatr = true;
    validate_newOrExistingMatrixAllocs(out, isNewMatr, __func__);

    return out;
}


extern "C" void destroyCompMatrN(CompMatrN matrix) {
    validate_matrixInit(matrix, __func__);

    // free each CPU row array
    for (qindex r=0; r < matrix.numRows; r++)
        free(matrix.elems[r]);

    // free CPU array of rows
    free(matrix.elems);

    // free flat GPU array if it exists
    if (matrix.gpuElems != NULL)
        gpu_deallocAmps(matrix.gpuElems);
}


extern "C" void syncCompMatrN(CompMatrN matr) {
    validate_matrixInit(matr, __func__);
    validate_matrixElemsDontContainUnsyncFlag(matr.elems[0][0], __func__);

    gpu_copyCpuToGpu(matr);
}



/*
 * EXPLICIT VARIABLE-SIZE MATRIX INITIALISERS
 *
 * all of which are demangled for C and C++ compatibility
 */


template <typename T> 
void validateAndSetCompMatrNElems(CompMatrN out, T elems, const char* caller) {
    validate_matrixInit(out, __func__);
    validate_matrixElemsDontContainUnsyncFlag(elems[0][0], caller);

    // serially copy values to CPU memory
    populateCompMatrElems(out.elems, elems, out.numRows);

    // overwrite GPU elements (including unsync flag)
    if (out.gpuElems != NULL)
        gpu_copyCpuToGpu(out);
}

extern "C" void setCompMatrNFromPtr(CompMatrN out, qcomp** elems) {

    validateAndSetCompMatrNElems(out, elems, __func__);
}

// the corresponding setCompMatrNFromArr() function must use VLAs and so
// is C++ incompatible, and is subsequently defined inline in the header file.
// Because it needs to create stack memory with size given by a CompMatrN field,
// we need to first validate that field via this exposed validation function. Blegh!

extern "C" void validate_setCompMatrNFromArr(CompMatrN out) {

    // the user likely invoked this function from the setLiteralCompMatrN()
    // macro, but we cannot know for sure so it's better to fall-back to
    // reporting the definitely-involved inner function, as we do elsewhere
    validate_matrixInit(out, "setCompMatrNFromArr");
}



/*
 * OVERLOADED VARIABLE-SIZE MATRIX INITIALISERS
 *
 * which are C++ only; equivalent C overloads are defined using
 * macros in the header file. Note the explicit overloads below
 * excludes a 2D qcomp[][] array because VLA is not supported by C++.
 */


void setCompMatrN(CompMatrN out, qcomp** in) {

    validateAndSetCompMatrNElems(out, in, __func__);
}

void setCompMatrN(CompMatrN out, std::vector<std::vector<qcomp>> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixInit(out, __func__);
    validate_numMatrixElems(1, in, __func__);

    validateAndSetCompMatrNElems(out, in, __func__);
}



/*
 * STRING PARSING VARIABLE-SIZE MATRIX INITIALISERS
 */


// de-mangle for both C and C++ invocation
extern "C" void setStringCompMatrN(CompMatrN matr, char* string) {
    validate_matrixInit(matr, __func__);

    // TODO: validate post-parsed string doesn't contain unsync flag

    // TODO: invoke parser
}



/*
 * C & C++ MATRIX REPORTERS
 *
 * and private (permittedly name-mangled) inner functions
 */


void printMatrixHeader(int numQubits) {

    // find memory used by matrix; equal to that of a non-distributed density matrix
    bool isMatr = true;
    int numNodes = 1;
    size_t mem = mem_getLocalMemoryRequired(numQubits, isMatr, numNodes);
    std::string memStr = form_str(mem) + by;

    // prepare dim substring
    qindex dim = powerOf2(numQubits);
    std::string dimStr = form_str(dim) + mu + form_str(dim);

    // e.g. CompMatr2 (4 x 4, 256 bytes):
    std::cout << "CompMatr" << numQubits << " (" << dimStr << ", " << memStr << "):" << std::endl;
}


template<class T> 
void rootPrintMatrix(T matrix) {
    
    if (comm_getRank() != 0)
        return;

    printMatrixHeader(matrix.numQubits);
    form_printMatrix(matrix);
}


extern "C" void reportCompMatr1(CompMatr1 matr) {

    rootPrintMatrix(matr);
}

extern "C" void reportCompMatr2(CompMatr2 matr) {

    rootPrintMatrix(matr);
}

extern "C" void reportCompMatrN(CompMatrN matr) {

    rootPrintMatrix(matr);
}

