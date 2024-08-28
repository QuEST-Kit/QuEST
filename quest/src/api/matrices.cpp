/** @file
 * Definitions of API matrix data structures, and their getters
 * and setters, as well as reporting utilities.
 * 
 * This file defines many "layers" of initialisation of complex
 * matrices, as explained in the header file.
 */

#include "quest/include/matrices.h"
#include "quest/include/environment.h"
#include "quest/include/types.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>

using std::vector;



/*
 * PRIVATE MATRIX-AGNOSTIC UTILITES 
 *
 * noting that some public matrix utilities (like util_isDenseMatrix) 
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


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
bool didAnyLocalAllocsFail(T matr) {

    // god help us if this single-integer malloc failed
    if (matr.isUnitary == NULL)
        return true;

    // outer CPU memory should always be allocated
    if (matr.cpuElems == NULL)
        return true;

    // if memory is 2D, we must also check each inner array was allocated
    if constexpr (util_isDenseMatrixType<T>())
        for (qindex r=0; r<matr.numRows; r++)
            if (matr.cpuElems[r] == NULL)
                return true;

    // if env is GPU-accelerated, we should have allocated persistent GPU memory
    if (getQuESTEnv().isGpuAccelerated && matr.gpuElems == NULL)
        return true;

    // otherwise, all pointers were non-NULL and ergo all allocs were successful
    return false;
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
bool didAnyAllocsFailOnAnyNode(T matr) {

    bool anyFail = didAnyLocalAllocsFail(matr);
    if (comm_isInit())
        anyFail = comm_isTrueOnAllNodes(anyFail);

    return anyFail;
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void freeAllMemoryIfAnyAllocsFailed(T matr) {

    // do nothing if everything allocated successfully between all nodes
    if (!didAnyAllocsFailOnAnyNode(matr))
        return;

    // otherwise free all successfully allocated rows of 2D structures (if outer list allocated)
    if constexpr (util_isDenseMatrixType<T>())
        if (matr.cpuElems != NULL)
            for (qindex r=0; r<matr.numRows; r++)
                free(matr.cpuElems[r]);

    // freeing NULL is legal
    free(matr.cpuElems);
    free(matr.isUnitary);
    gpu_deallocAmps(matr.gpuElems);
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void validateMatrixAllocsAndSetUnitarity(T matr, size_t numBytes, const char* caller) {

    // free memory before throwing validation error to avoid memory leaks
    freeAllMemoryIfAnyAllocsFailed(matr);
    validate_newMatrixAllocs(matr, numBytes, caller);

    // set initial unitarity of the newly created matrix to unknown
    if (matr.isUnitary != NULL)
        *(matr.isUnitary) = validate_UNITARITY_UNKNOWN_FLAG;
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void validateAndSyncMatrix(T matr, const char* caller) {
    validate_matrixFields(matr, caller);
    validate_matrixNewElemsDontContainUnsyncFlag(util_getFirstLocalElem(matr), caller);

    // optionally overwrite GPU elements with user-modified CPU elements
    if (matr.gpuElems != NULL)
        gpu_copyCpuToGpu(matr);

    // optionally determine unitarity
    if (validate_isEnabled())
        *(matr.isUnitary) = util_isUnitary(matr);
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void validateAndDestroyMatrix(T matrix, const char* caller) {
    validate_matrixFields(matrix, caller);

    // free each CPU row array (if matrix is 2D)
    if constexpr (util_isDenseMatrixType<T>())
        for (qindex r=0; r < matrix.numRows; r++)
            free(matrix.cpuElems[r]);

    // free CPU array of rows (if 2D) or elems (if 1D)
    free(matrix.cpuElems);

    // free flat GPU array if it exists
    if (matrix.gpuElems != NULL)
        gpu_deallocAmps(matrix.gpuElems);
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
CompMatr1 getCompMatr1(vector<vector<qcomp>> in) {
    validate_matrixNumNewElems(1, in, __func__);

    return getCompMatrFromElems<CompMatr1>(in, 1);
}

CompMatr2 getCompMatr2(qcomp in[4][4]) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2(qcomp** in) {
    return getCompMatrFromElems<CompMatr2>(in, 2);
}
CompMatr2 getCompMatr2(vector<vector<qcomp>> in) {
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


DiagMatr1 getDiagMatr1(vector<qcomp> in) {
    validate_matrixNumNewElems(1, in, __func__);

    return getDiagMatr1(in.data());
}

DiagMatr2 getDiagMatr2(vector<qcomp> in) {
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
    validate_newCompMatrParams(numQubits, __func__);

    // validation ensures these (and below mem sizes) never overflow
    qindex numRows = powerOf2(numQubits);
    qindex numElems = numRows * numRows;
    qindex numBytes = numElems * sizeof(qcomp); // used only by validation report

    // initialise all CompMatr fields inline because struct is const
    CompMatr out = {
        .numQubits = numQubits,
        .numRows = numRows,

        // store unitary flag in the heap so that struct copies are mutable
        .isUnitary = (int*) malloc(sizeof *out.isUnitary), // NULL if failed

        // const 2D CPU memory (NULL if failed, or containing NULLs)
        .cpuElems = (qcomp**) malloc(numRows * sizeof *out.cpuElems), // NULL if failed,

        // const 1D GPU memory (NULL if failed or not needed)
        .gpuElems = (getQuESTEnv().isGpuAccelerated)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    // only if outer CPU allocation succeeded, attempt to allocate each row array
    if (out.cpuElems != NULL)
        for (qindex r=0; r < numRows; r++)
            out.cpuElems[r] = cpu_allocAmps(numRows); // NULL if failed

    validateMatrixAllocsAndSetUnitarity(out, numBytes, __func__);
    return out;
}


extern "C" void destroyCompMatr(CompMatr matr) {

    validateAndDestroyMatrix(matr, __func__);
}


extern "C" void syncCompMatr(CompMatr matr) {

    validateAndSyncMatrix(matr, __func__);
}



/*
 * VARIABLE SIZE DIAGONAL MATRIX CONSTRUCTORS
 *
 * all of which are de-mangled for both C++ and C compatibility
 */


extern "C" DiagMatr createDiagMatr(int numQubits) {
    validate_envIsInit(__func__);
    validate_newDiagMatrParams(numQubits, __func__);

    // validation ensures these (and below mem sizes) never overflow
    qindex numElems = powerOf2(numQubits);
    qindex numBytes = numElems * sizeof(qcomp);

    // initialise all CompMatr fields inline because struct is const
    DiagMatr out = {
        .numQubits = numQubits,
        .numElems = numElems,

        // store unitary flag in the heap so that struct copies are mutable
        .isUnitary = (int*) malloc(sizeof *out.isUnitary), // NULL if failed

        // const 1D CPU memory (NULL if failed, or containing NULLs)
        .cpuElems = cpu_allocAmps(numElems), // NULL if failed,

        // const 1D GPU memory (NULL if failed or not needed)
        .gpuElems = (getQuESTEnv().isGpuAccelerated)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    validateMatrixAllocsAndSetUnitarity(out, numBytes, __func__);
    return out;
}


extern "C" void destroyDiagMatr(DiagMatr matr) {

    validateAndDestroyMatrix(matr, __func__);
}


extern "C" void syncDiagMatr(DiagMatr matr) {

    validateAndSyncMatrix(matr, __func__);
}



/*
 * FULL-STATE MATRIX CONSTRUCTORS
 *
 * the public of which are de-mangled for both C++ and C compatibility
 */


FullStateDiagMatr validateAndCreateCustomFullStateDiagMatr(int numQubits, int useDistrib, const char* caller) {
    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    // must validate parameters before passing them to autodeployer
    validate_newFullStateDiagMatrParams(numQubits, useDistrib, caller);

    // overwrite useDistrib if it was left as AUTO_FLAG
    autodep_chooseFullStateDiagMatrDeployment(numQubits, useDistrib, env);

    qindex numElems = powerOf2(numQubits);
    qindex numElemsPerNode = numElems / (useDistrib? env.numNodes : 1); // divides evenly
    qindex numBytesPerNode = numElemsPerNode * sizeof(qcomp);

    FullStateDiagMatr out = {

        // data deployment configuration; disable distrib if deployed to 1 node
        .isDistributed = useDistrib && (env.numNodes > 1),

        .numQubits = numQubits,
        .numElems = numElems,
        .numElemsPerNode = numElemsPerNode,

        // store unitary flag in the heap so that struct copies are mutable
        .isUnitary = (int*) malloc(sizeof *out.isUnitary), // NULL if failed

        // 1D CPU memory (NULL if failed)
        .cpuElems = cpu_allocAmps(numElemsPerNode),

        // 1D GPU memory (NULL if failed or deliberately not allocated)
        .gpuElems = (env.isGpuAccelerated)? gpu_allocAmps(numElemsPerNode) : NULL, // first amp will be un-sync'd flag
    };

    validateMatrixAllocsAndSetUnitarity(out, numBytesPerNode, __func__);
    return out;
}


extern "C" FullStateDiagMatr createCustomFullStateDiagMatr(int numQubits, int useDistrib) {

    return validateAndCreateCustomFullStateDiagMatr(numQubits, useDistrib, __func__);
}


extern "C" FullStateDiagMatr createFullStateDiagMatr(int numQubits) {

    return validateAndCreateCustomFullStateDiagMatr(numQubits, modeflag::USE_AUTO, __func__);
}


extern "C" void destroyFullStateDiagMatr(FullStateDiagMatr matr) {

    validateAndDestroyMatrix(matr, __func__);
}


extern "C" void syncFullStateDiagMatr(FullStateDiagMatr matr) {

    validateAndSyncMatrix(matr, __func__);
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

    // overwrite GPU elements (including unsync flag); validation gauranteed to pass
    syncCompMatr(out); 
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

void setCompMatr(CompMatr out, vector<vector<qcomp>> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    // then we unimportantly repeat some of this validation; alas!
    validateAndSetCompMatrElems(out, in, __func__);
}



/*
 * OVERLOADED VARIABLE-SIZE DIAGONAL MATRIX INITIALISERS
 *
 * visible to both C and C++, although C++ additionally gets vector overloads.
 */


extern "C" void setDiagMatr(DiagMatr out, qcomp* in) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(in[0], __func__);

    // overwrite CPU memory
    memcpy(out.cpuElems, in, out.numElems * sizeof(qcomp));

    // overwrite GPU elements (including unsync flag); validation gauranteed to pass
    syncDiagMatr(out);
}

void setDiagMatr(DiagMatr out, vector<qcomp> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    // then we unimportantly repeat some of this validation; alas!
    setDiagMatr(out, in.data());
}


extern "C" void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(in[0], __func__);
    validate_fullStateDiagMatrNewElems(out, startInd, numElems, __func__);

    // if the matrix is non-distributed, we update every node's duplicated CPU amps
    if (!out.isDistributed)
        memcpy(&out.cpuElems[startInd], in, numElems * sizeof(qcomp));

    // only distributed nodes containing targeted elemsn need to do anything
    else if (util_areAnyElemsWithinThisNode(out.numElemsPerNode, startInd, numElems)) {
        util_IndexRange range = util_getLocalIndRangeOfElemsWithinThisNode(out.numElemsPerNode, startInd, numElems);
        memcpy(&out.cpuElems[range.localDistribStartInd], &in[range.localDuplicStartInd], range.numElems * sizeof(qcomp));
    }

    // all nodes overwrite GPU; validation gauranteed to succeed
    syncFullStateDiagMatr(out);
}

void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, vector<qcomp> in) {

    setFullStateDiagMatr(out, startInd, in.data(), in.size());
}



/*
 * MATRIX REPORTERS
 *
 * and their private (permittedly name-mangled) inner functions
 */


// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T>
void printMatrixHeader(T matr) {

    // determine how many nodes the matrix is distributed between (if it's a distributable type)
    int numNodes = (util_isDistributedMatrix(matr))? getQuESTEnv().numNodes : 1;

    // find memory (bytes) to store elements
    size_t elemMem = mem_getLocalMatrixMemoryRequired(matr.numQubits, util_isDenseMatrixType<T>(), numNodes);

    // find memory (bytes) of other struct fields; fixed-size sizeof includes arrays, var-size does not
    size_t otherMem = sizeof(matr);
    if (util_isFixedSizeMatrixType<T>())
        otherMem -= elemMem;

    print_matrixInfo(util_getMatrixTypeName<T>(), matr.numQubits, util_getMatrixDim(matr), elemMem, otherMem, numNodes);
}


// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template<class T> 
void validateAndPrintMatrix(T matr, const char* caller) {
    validate_matrixFields(matr, caller);

    // syncable matrices must be synced before reporting (though only CPU elems are printed)
    if constexpr (!util_isFixedSizeMatrixType<T>())
        validate_matrixIsSynced(matr, caller);

    // print_ functions will handle distributed coordination
    printMatrixHeader(matr);
    print_matrix(matr);
}


// all reporters are C and C++ accessible, so are de-mangled
extern "C" {
    
    void reportCompMatr1(CompMatr1 matr) {
        validateAndPrintMatrix(matr, __func__);
    }
    void reportCompMatr2(CompMatr2 matr) {
        validateAndPrintMatrix(matr, __func__);
    }
    void reportCompMatr(CompMatr matr) {
        validateAndPrintMatrix(matr, __func__);
    }
    void reportDiagMatr1(DiagMatr1 matr) {
        validateAndPrintMatrix(matr, __func__);
    }
    void reportDiagMatr2(DiagMatr2 matr) {
        validateAndPrintMatrix(matr, __func__);
    }
    void reportDiagMatr(DiagMatr matr) {
        validateAndPrintMatrix(matr, __func__);
    }
    void reportFullStateDiagMatr(FullStateDiagMatr matr) {
        validateAndPrintMatrix(matr, __func__);
    }
}
