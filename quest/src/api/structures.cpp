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
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>



/*
 * PRIVATE UTILITES 
 *
 * noting that some public matrix utilities (like isDenseMatrix) 
 * are defined in utilities.hpp
 */


// types A and B can both be qcomp** or qcomp[][] (mixed)
template<typename A, typename B> 
void populateCompMatrElems(A out, B in, qindex dim) {

    for (qindex r=0; r<dim; r++)
        for (qindex c=0; c<dim; c++)
            out[r][c] = in[r][c];
}


// type T can be CompMatr1/2, B can be qcomp** or qcomp[][]
template<class T, typename B>
T getCompMatrFromElems(B in, int num) {

    // must initialize the const fields inline
    T out = (T) {
        .numQubits = num,
        .numRows = powerOf2(num)
    };

    // elems is pre-allocated and non-const (elements can be modified)
    populateCompMatrElems(out.elems, in, out.numRows);
    return out;
}


// type T can be CompMatr or DiagMatr
template <class T>
bool didAnyAllocsFail(T matr) {

    // outer CPU memory should always be allocated
    if (matr.cpuElems == NULL)
        return true;

    // if memory is 2D, we must also check each inner array was allocated
    if constexpr (isDenseMatrixType<T>())
        for (qindex r=0; r<matr.numRows; r++)
            if (matr.cpuElems[r] == NULL)
                return true;

    // if env is GPU-accelerated, we should have allocated persistent GPU memory
    if (getQuESTEnv().isGpuAccelerated && matr.gpuElems == NULL)
        return true;

    // otherwise, all pointers were non-NULL and ergo all allocs were successful
    return false;
}


// type T can be CompMatr or DiagMatr
template <class T>
void freeAllMemoryIfAnyAllocsFailed(T matr) {

    // do nothing if everything allocated successfully
    if (!didAnyAllocsFail(matr))
        return;

    // otherwise, free all successfully allocated rows of 2D structures (if outer list allocated)
    if constexpr (isDenseMatrixType<T>())
        if (matr.cpuElems != NULL)
            for (qindex r=0; r<matr.numRows; r++)
                if (matr.cpuElems[r] != NULL)
                    free(matr.cpuElems[r]);

    // free the outer CPU array itself
    if (matr.cpuElems != NULL)
        free(matr.cpuElems);
    
    // and the GPU memory (gauranteed NULL in non-GPU mode)
    if (matr.gpuElems != NULL)
        gpu_deallocAmps(matr.gpuElems);
}



/*
 * OVERLOADED FIXED-SIZE COMPLEX MATRIX INITIALISERS
 *
 * enabling getCompMatr1/2 to receive qcomp**, qcomp(*)[] = qcomp[][]
 * (due to pointer decay) and vector<vector<qcomp>>. These signatures
 * are only exposed to C++; equivalent C macros are defined in the header.
 * The vector signatures permit inline initialisers, and can also
 * validate the dimensions of the passed lists.
 */


CompMatr1 getCompMatr1(qcomp in[2][2]) {
    return getCompMatrFromElems<CompMatr1>(in, 1);
}
CompMatr1 getCompMatr1(qcomp** in) {
    return getCompMatrFromElems<CompMatr1>(in, 1);
}
CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in) {
    validate_matrixNumNewElems(1, in, __func__);

    return getCompMatrFromElems<CompMatr1>(in, 1);
}

CompMatr2 getCompMatr2(qcomp in[4][4]) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2(qcomp** in) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in) {
    validate_matrixNumNewElems(2, in, __func__);

    return getCompMatrFromElems<CompMatr2>(in, 2);
}



/*
 * OVERLOADED FIXED-SIZE DIAGONAL MATRIX INITIALISERS
 *
 * enabling getDiagMatr1/2 to receive qcomp* = qcomp[], or vector<vector<qcomp>>. 
 * Only the C++-compatible vector definition is here; the C-compatible pointer 
 * overloads are defined in the header to circumvent C & C++ qcomp interopability 
 * issues. The vector signatures permit inline initialisers, and can also
 * validate the dimensions of the passed lists.
 */


DiagMatr1 getDiagMatr1(std::vector<qcomp> in) {
    validate_matrixNumNewElems(1, in, __func__);
    return getDiagMatr1(in.data());
}

DiagMatr2 getDiagMatr2(std::vector<qcomp> in) {
    validate_matrixNumNewElems(2, in, __func__);
    return getDiagMatr2(in.data());
}



/*
 * VARIABLE SIZE COMPLEX MATRIX CONSTRUCTORS
 *
 * all of which are de-mangled for both C++ and C compatibility
 */


extern "C" CompMatr createCompMatr(int numQubits) {
    validate_envIsInit(__func__);
    validate_newMatrixNumQubits(numQubits, __func__);

    // validation ensures these (and below mem sizes) never overflow
    qindex numRows = powerOf2(numQubits);
    qindex numElems = numRows * numRows;
    qindex numBytes = numElems * sizeof(qcomp); // used only by validation report

    // we will always allocate GPU memory if the env is GPU-accelerated
    bool isGpuAccel = getQuESTEnv().isGpuAccelerated;

    // initialise all CompMatr fields inline because struct is const
    CompMatr out = {
        .numQubits = numQubits,
        .numRows = numRows,

        // const 2D CPU memory (NULL if failed, or containing NULLs)
        .cpuElems = (qcomp**) malloc(numRows * sizeof *out.cpuElems), // NULL if failed,

        // const 1D GPU memory (NULL if failed or not needed)
        .gpuElems = (isGpuAccel)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    // only if outer CPU allocation succeeded, attempt to allocate each row array
    if (out.cpuElems != NULL)
        for (qindex r=0; r < numRows; r++)
            out.cpuElems[r] = (qcomp*) calloc(numRows, sizeof **out.cpuElems); // NULL if failed

    // if any of the above mallocs failed, below validation will memory leak; so free first (but don't set to NULL)
    freeAllMemoryIfAnyAllocsFailed(out);
    validate_newMatrixAllocs(out, numBytes, __func__);

    return out;
}


extern "C" void destroyCompMatr(CompMatr matrix) {
    validate_matrixFields(matrix, __func__);

    // free each CPU row array
    for (qindex r=0; r < matrix.numRows; r++)
        free(matrix.cpuElems[r]);

    // free CPU array of rows
    free(matrix.cpuElems);

    // free flat GPU array if it exists
    if (matrix.gpuElems != NULL)
        gpu_deallocAmps(matrix.gpuElems);
}


extern "C" void syncCompMatr(CompMatr matr) {
    validate_matrixFields(matr, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(matr.cpuElems[0][0], __func__);

    gpu_copyCpuToGpu(matr);
}



/*
 * VARIABLE SIZE DIAGONAL MATRIX CONSTRUCTORS
 *
 * all of which are de-mangled for both C++ and C compatibility
 */


extern "C" DiagMatr createDiagMatr(int numQubits) {
    validate_envIsInit(__func__);
    validate_newMatrixNumQubits(numQubits, __func__);

    // validation ensures these (and below mem sizes) never overflow
    qindex numElems = powerOf2(numQubits);
    qindex numBytes = numElems * sizeof(qcomp);

    // we will always allocate GPU memory if the env is GPU-accelerated
    bool isGpuAccel = getQuESTEnv().isGpuAccelerated;

    // initialise all CompMatr fields inline because struct is const
    DiagMatr out = {
        .numQubits = numQubits,
        .numElems = numElems,

        // const 2D CPU memory (NULL if failed, or containing NULLs)
        .cpuElems = (qcomp*) malloc(numElems * sizeof *out.cpuElems), // NULL if failed,

        // const 1D GPU memory (NULL if failed or not needed)
        .gpuElems = (isGpuAccel)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    // if either of these mallocs failed, below validation will memory leak; so free first (but don't set to NULL)
    freeAllMemoryIfAnyAllocsFailed(out);

    // check all CPU & GPU malloc and callocs succeeded (no attempted freeing if not)
    validate_newMatrixAllocs(out, numBytes, __func__);

    return out;
}


extern "C" void destroyDiagMatr(DiagMatr matrix) {
    validate_matrixFields(matrix, __func__);

    // free CPU array of rows
    free(matrix.cpuElems);

    // free flat GPU array if it exists
    if (matrix.gpuElems != NULL)
        gpu_deallocAmps(matrix.gpuElems);
}


extern "C" void syncDiagMatr(DiagMatr matr) {
    validate_matrixFields(matr, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(matr.cpuElems[0], __func__);

    gpu_copyCpuToGpu(matr);
}



/*
 * EXPLICIT VARIABLE-SIZE COMPLEX MATRIX INITIALISERS
 *
 * all of which are demangled for C and C++ compatibility.
 * Similar functions for DiagMatr are not necessary since
 * it contains a 1D pointer (which arrays decay to).
 */


// type T can be qcomp** or vector<vector<qcomp>>, but qcomp(*)[] is handled by header
template <typename T> 
void validateAndSetCompMatrElems(CompMatr out, T elems, const char* caller) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(elems[0][0], caller);

    // serially copy values to CPU memory
    populateCompMatrElems(out.cpuElems, elems, out.numRows);

    // overwrite GPU elements (including unsync flag)
    if (out.gpuElems != NULL)
        gpu_copyCpuToGpu(out);
}

extern "C" void setCompMatrFromPtr(CompMatr out, qcomp** elems) {

    validateAndSetCompMatrElems(out, elems, __func__);
}

// the corresponding setCompMatrFromArr() function must use VLAs and so
// is C++ incompatible, and is subsequently defined inline in the header file.
// Because it needs to create stack memory with size given by a CompMatr field,
// we need to first validate that field via this exposed validation function. Blegh!

extern "C" void validate_setCompMatrFromArr(CompMatr out) {

    // the user likely invoked this function from the setInlineCompMatr()
    // macro, but we cannot know for sure so it's better to fall-back to
    // reporting the definitely-involved inner function, as we do elsewhere
    validate_matrixFields(out, "setCompMatrFromArr");
}



/*
 * OVERLOADED VARIABLE-SIZE COMPLEX MATRIX INITIALISERS
 *
 * which are C++ only; equivalent C overloads are defined using
 * macros in the header file. Note the explicit overloads below
 * excludes a 2D qcomp[][] array because VLA is not supported by C++.
 */


void setCompMatr(CompMatr out, qcomp** in) {

    validateAndSetCompMatrElems(out, in, __func__);
}

void setCompMatr(CompMatr out, std::vector<std::vector<qcomp>> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    // then we unimportantly repeat some of this validation; alas!
    validateAndSetCompMatrElems(out, in, __func__);
}



/*
 * OVERLOADED VARIABLE-SIZE DIAGONAL MATRIX INITIALISERS
 *
 * visible to both C and C++, although C++ additionally gets a vector overload.
 */


extern "C" void setDiagMatr(DiagMatr out, qcomp* in) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(in[0], __func__);

    // overwrite CPU memory
    memcpy(out.cpuElems, in, out.numElems * sizeof(qcomp));

    // overwrite GPU elements (including unsync flag)
    if (out.gpuElems != NULL)
        gpu_copyCpuToGpu(out);
}

void setDiagMatr(DiagMatr out, std::vector<qcomp> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    // then we unimportantly repeat some of this validation; alas!
    setDiagMatr(out, in.data());
}



/*
 * C & C++ MATRIX REPORTERS
 *
 * and their private (permittedly name-mangled) inner functions
 */


// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2 or DiagMatr
template <class T>
void printMatrixHeader(T matr) {

    // produce exact type string, e.g. DiagMatr2
    std::string nameStr = (isDenseMatrixType<T>())? "CompMatr" : "DiagMatr";
    if (isFixedSizeMatrixType<T>())
        nameStr += form_str(matr.numQubits);

    // find memory (bytes) to store elements; equal to that of a non-distributed Qureg (statevec or dens-matr)
    size_t elemMem = mem_getLocalMemoryRequired(matr.numQubits, isDenseMatrixType<T>(), 1); // 1 node

    // find memory (bytes) of other struct fields; fixed-size sizeof includes arrays, var-size does not
    size_t otherMem = sizeof(matr);
    if (isFixedSizeMatrixType<T>())
        otherMem -= elemMem;

    // find dimension; illegal field access in wrong branch removed at compile-time
    qindex dim;
    if constexpr (isDenseMatrixType<T>())
        dim = matr.numRows;
    else
        dim = matr.numElems;

    form_printMatrixInfo(nameStr, matr.numQubits, dim, elemMem, otherMem);
}


// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2 or DiagMatr
template<class T> 
void rootPrintMatrix(T matrix) {
    
    // only root note prints, to avoid spam in distributed settings
    if (comm_getRank() != 0)
        return;

    printMatrixHeader(matrix);
    form_printMatrix(matrix);
}


// all reporters are C and C++ accessible, so are de-mangled
extern "C" {
    
    void reportCompMatr1(CompMatr1 matr) {
        rootPrintMatrix(matr);
    }
    void reportCompMatr2(CompMatr2 matr) {
        rootPrintMatrix(matr);
    }
    void reportCompMatr(CompMatr matr) {
        rootPrintMatrix(matr);
    }
    void reportDiagMatr1(DiagMatr1 matr) {
        rootPrintMatrix(matr);
    }
    void reportDiagMatr2(DiagMatr2 matr) {
        rootPrintMatrix(matr);
    }
    void reportDiagMatr(DiagMatr matr) {
        rootPrintMatrix(matr);
    }
}
