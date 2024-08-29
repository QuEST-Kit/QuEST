/** @file
 * Definitions of API matrix data structures, and their getters
 * and setters, as well as reporting utilities.
 * 
 * This file defines several layers of initialisation of complex
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
 * FIZED-SIZE MATRIX VECTOR GETTERS
 *
 * enabling getCompMatr1 (etc) to receive vectors, in addition
 * to the existing pointer and array overloads in the header.
 * These are necessarily defined here in the source, unlike the
 * other header functions, so that they may use private valdiation.
 */


CompMatr1 getCompMatr1(vector<vector<qcomp>> in) {
    validate_matrixNumNewElems(1, in, __func__);

    qcomp* rowPtrs[] = {in[0].data(), in[1].data()};
    return getCompMatr1(rowPtrs);
}

CompMatr2 getCompMatr2(vector<vector<qcomp>> in) {
    validate_matrixNumNewElems(2, in, __func__);

    qcomp* rowPtrs[] = {in[0].data(), in[1].data(), in[2].data(), in[3].data()};
    return getCompMatr2(rowPtrs);
}


DiagMatr1 getDiagMatr1(vector<qcomp> in) {
    validate_matrixNumNewElems(1, in, __func__);

    return getDiagMatr1(in.data());
}

DiagMatr2 getDiagMatr2(vector<qcomp> in) {
    validate_matrixNumNewElems(2, in, __func__);

    return getDiagMatr2(in.data());
}



/*
 * PRIVATE FUNCTIONS OF VARIABLE-SIZE MATRIX CONSTRUCTORS
 */


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
    if (getQuESTEnv().isGpuAccelerated && util_getGpuMemPtr(matr) == NULL)
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
    gpu_deallocAmps(util_getGpuMemPtr(matr));
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void validateMatrixAllocsAndSetUnitarity(T matr, size_t numBytes, const char* caller) {

    // free memory before throwing validation error to avoid memory leaks
    freeAllMemoryIfAnyAllocsFailed(matr);
    validate_newMatrixAllocs(matr, numBytes, caller);

    // set initial unitarity of the newly created matrix to unknown
    *(matr.isUnitary) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
}



/*
 * VARIABLE-SIZE MATRIX CONSTRUCTORS
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
        .gpuElemsFlat = (getQuESTEnv().isGpuAccelerated)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    // only if outer CPU allocation succeeded, attempt to allocate each row array
    if (out.cpuElems != NULL)
        for (qindex r=0; r < numRows; r++)
            out.cpuElems[r] = cpu_allocAmps(numRows); // NULL if failed

    validateMatrixAllocsAndSetUnitarity(out, numBytes, __func__);
    return out;
}


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



/*
 * VARIABLE-SIZE MATRIX SYNC
 */


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void validateAndSyncMatrix(T matr, const char* caller) {
    validate_matrixFields(matr, caller);
    validate_matrixNewElemsDontContainUnsyncFlag(util_getFirstLocalElem(matr), caller);

    // optionally overwrite GPU elements with user-modified CPU elements
    if (util_getGpuMemPtr(matr) != NULL)
        gpu_copyCpuToGpu(matr);

    // indicate that we do not know whether the revised matrix is
    // is unitarity; we defer establishing that until a unitarity check
    *(matr.isUnitary) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
}


// de-mangled for both C++ and C compatibility
extern "C" {

    void syncCompMatr(CompMatr matr) { validateAndSyncMatrix(matr, __func__); }
    void syncDiagMatr(DiagMatr matr) { validateAndSyncMatrix(matr, __func__); }
    void syncFullStateDiagMatr(FullStateDiagMatr matr) { validateAndSyncMatrix(matr, __func__); }

}



/*
 * VARIABLE-SIZE MATRIX DESTRUCTION
 */


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
    if (util_getGpuMemPtr(matrix) != NULL)
        gpu_deallocAmps(util_getGpuMemPtr(matrix));
}


// de-mangled for C++ and C compatibility
extern "C" {
    
    void destroyCompMatr(CompMatr matr) { validateAndDestroyMatrix(matr, __func__); }
    void destroyDiagMatr(DiagMatr matr) { validateAndDestroyMatrix(matr, __func__); }
    void destroyFullStateDiagMatr(FullStateDiagMatr matr) { validateAndDestroyMatrix(matr, __func__); }

}



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA POINTERS
 */


// type T can be qcomp** or vector<vector<qcomp>>, but qcomp(*)[] is handled by header
template <typename T> 
void validateAndSetDenseMatrElems(CompMatr out, T elems, const char* caller) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(elems[0][0], caller);

    // copy elems into matrice's CPU memory
    for (qindex r=0; r<out.numRows; r++)
        for (qindex c=0; c<out.numRows; c++)
            out.cpuElems[r][c] = elems[r][c];

    // overwrite GPU elements (including unsync flag); validation gauranteed to pass
    syncCompMatr(out); 
}


extern "C" void setCompMatr(CompMatr out, qcomp** in) {

    // use the above template, which we will reuse for the vector overloads
    validateAndSetDenseMatrElems(out, in, __func__);
}


extern "C" void setDiagMatr(DiagMatr out, qcomp* in) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(in[0], __func__);

    // overwrite CPU memory
    memcpy(out.cpuElems, in, out.numElems * sizeof(qcomp));

    // overwrite GPU elements (including unsync flag); validation gauranteed to pass
    syncDiagMatr(out);
}


extern "C" void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems) {
    validate_matrixFields(out, __func__);
    validate_matrixNewElemsDontContainUnsyncFlag(in[0], __func__);
    validate_fullStateDiagMatrNewElems(out, startInd, numElems, __func__);

    // if the matrix is non-distributed, we update every node's duplicated CPU amps
    if (!out.isDistributed)
        memcpy(&out.cpuElems[startInd], in, numElems * sizeof(qcomp));

    // only distributed nodes containing targeted elements need to do anything
    else if (util_areAnyElemsWithinThisNode(out.numElemsPerNode, startInd, numElems)) {
        util_IndexRange range = util_getLocalIndRangeOfElemsWithinThisNode(out.numElemsPerNode, startInd, numElems);
        memcpy(&out.cpuElems[range.localDistribStartInd], &in[range.localDuplicStartInd], range.numElems * sizeof(qcomp));
    }

    // all nodes overwrite GPU; validation gauranteed to succeed
    syncFullStateDiagMatr(out);
}



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA VECTORS
 */


void setCompMatr(CompMatr out, vector<vector<qcomp>> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    // then we unimportantly repeat some of this validation; alas!
    validateAndSetDenseMatrElems(out, in, __func__);
}


void setDiagMatr(DiagMatr out, vector<qcomp> in) {

    // we validate dimension of 'in', which first requires validating 'out' fields
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    // then we unimportantly repeat some of this validation; alas!
    setDiagMatr(out, in.data());
}


void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, vector<qcomp> in) {

    setFullStateDiagMatr(out, startInd, in.data(), in.size());
}


// no bespoke array functions are necessary for diagonal matrices initialisation, 
// since passed arrays automatically decay to pointers



/*
 * MATRIX REPORTERS
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
    
    void reportCompMatr1(CompMatr1 matr) { validateAndPrintMatrix(matr, __func__); }
    void reportCompMatr2(CompMatr2 matr) { validateAndPrintMatrix(matr, __func__); }
    void reportCompMatr (CompMatr  matr) { validateAndPrintMatrix(matr, __func__); }
    void reportDiagMatr1(DiagMatr1 matr) { validateAndPrintMatrix(matr, __func__); }
    void reportDiagMatr2(DiagMatr2 matr) { validateAndPrintMatrix(matr, __func__); }
    void reportDiagMatr (DiagMatr  matr) { validateAndPrintMatrix(matr, __func__); }
    void reportFullStateDiagMatr(FullStateDiagMatr matr) { validateAndPrintMatrix(matr, __func__); }

}
