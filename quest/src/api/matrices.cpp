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
 * PRIVATE VARIABLE-SIZE MATRIX UTILITIES
 */


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void freeHeapMatrix(T matr) {

    // free the 1D or 2D matrix - safe even if nullptr
    if constexpr (util_isDenseMatrixType<T>())
        cpu_deallocMatrix(matr.cpuElems, matr.numRows);
    else
        cpu_deallocArray(matr.cpuElems);

    // we avoid invoking a GPU function in non-GPU mode
    auto gpuPtr = util_getGpuMemPtr(matr);
    if (mem_isAllocated(gpuPtr))
        gpu_deallocArray(gpuPtr);

    // free the teeny tiny heap flags
    cpu_deallocHeapFlag(matr.isUnitary);
    cpu_deallocHeapFlag(matr.wasGpuSynced);
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
bool didAnyLocalAllocsFail(T matr) {

    // god help us if these single-integer malloc failed
    if (!mem_isAllocated(matr.isUnitary))
        return true;
    if (!mem_isAllocated(matr.wasGpuSynced))
        return true;

    // outer CPU memory should always be allocated
    if constexpr (util_isDenseMatrixType<T>()) {
        if (!mem_isAllocated(matr.cpuElems, matr.numRows))
            return true;
    } else {
        if (!mem_isAllocated(matr.cpuElems))
            return true;
    }

    // if memory is 2D, we must also check each inner array was allocated
    if constexpr (util_isDenseMatrixType<T>()) {
        if (!mem_isAllocated(matr.cpuElems, matr.numRows))
            return true;
    } else {
        if (!mem_isAllocated(matr.cpuElems))
            return true;
    }

    // if env is GPU-accelerated, we should have allocated persistent GPU memory
    if (getQuESTEnv().isGpuAccelerated && !mem_isAllocated(util_getGpuMemPtr(matr)))
        return true;

    // otherwise, all pointers were non-NULL and ergo all allocs were successful
    return false;
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void freeAllMemoryIfAnyAllocsFailed(T matr) {

    // ascertain whether any allocs failed on any node
    bool anyFail = didAnyLocalAllocsFail(matr);
    if (comm_isInit())
        anyFail = comm_isTrueOnAllNodes(anyFail);

    // if so, free all heap fields
    if (anyFail)
        freeHeapMatrix(matr);
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void validateMatrixAllocs(T matr, const char* caller) {

    // free memory before throwing validation error to avoid memory leaks
    freeAllMemoryIfAnyAllocsFailed(matr);
    validate_newMatrixAllocs(matr, caller);
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void setInitialHeapFlags(T matr) {

    // set initial unitarity of the newly created matrix to unknown
    *(matr.isUnitary) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;

    // indicate that GPU memory has not yet been synchronised
    *(matr.wasGpuSynced) = 0;
}



/*
 * VARIABLE-SIZE MATRIX CONSTRUCTORS
 */


extern "C" CompMatr createCompMatr(int numQubits) {
    validate_envIsInit(__func__);
    validate_newCompMatrParams(numQubits, __func__);

    // validation ensures these never overflow
    qindex numRows = powerOf2(numQubits);
    qindex numElems = numRows * numRows;

    // initialise all CompMatr fields inline because most are const
    CompMatr out = {
        .numQubits = numQubits,
        .numRows = numRows,

        // allocate flags in the heap so that struct copies are mutable
        .isUnitary = cpu_allocHeapFlag(), // nullptr if failed
        .wasGpuSynced = cpu_allocHeapFlag(), // nullptr if failed

        // 2D CPU memory
        .cpuElems = cpu_allocMatrix(numRows), // nullptr if failed, or may contain nullptr

        // 1D GPU memory
        .gpuElemsFlat = (getQuESTEnv().isGpuAccelerated)? gpu_allocArray(numElems) : nullptr // nullptr if failed or not needed
    };

    validateMatrixAllocs(out, __func__);
    setInitialHeapFlags(out);
    return out;
}


extern "C" DiagMatr createDiagMatr(int numQubits) {
    validate_envIsInit(__func__);
    validate_newDiagMatrParams(numQubits, __func__);

    // validation ensures this never overflows
    qindex numElems = powerOf2(numQubits);

    // initialise all CompMatr fields inline because most are const
    DiagMatr out = {
        .numQubits = numQubits,
        .numElems = numElems,

        // allocate flags in the heap so that struct copies are mutable
        .isUnitary = cpu_allocHeapFlag(), // nullptr if failed
        .wasGpuSynced = cpu_allocHeapFlag(), // nullptr if failed

        // 1D CPU memory
        .cpuElems = cpu_allocArray(numElems), // nullptr if failed

        // 1D GPU memory
        .gpuElems = (getQuESTEnv().isGpuAccelerated)? gpu_allocArray(numElems) : nullptr // nullptr if failed or not needed
    };

    validateMatrixAllocs(out, __func__);
    setInitialHeapFlags(out);
    return out;
}


FullStateDiagMatr validateAndCreateCustomFullStateDiagMatr(int numQubits, int useDistrib, const char* caller) {
    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    // must validate parameters before passing them to autodeployer
    validate_newFullStateDiagMatrParams(numQubits, useDistrib, caller);

    // overwrite useDistrib if it was left as AUTO_FLAG
    autodep_chooseFullStateDiagMatrDeployment(numQubits, useDistrib, env);

    // validation ensures this never overflows
    qindex numElems = powerOf2(numQubits);
    qindex numElemsPerNode = numElems / (useDistrib? env.numNodes : 1); // divides evenly

    FullStateDiagMatr out = {

        // data deployment configuration; disable distrib if deployed to 1 node
        .isDistributed = useDistrib && (env.numNodes > 1),

        .numQubits = numQubits,
        .numElems = numElems,
        .numElemsPerNode = numElemsPerNode,

        // allocate flags in the heap so that struct copies are mutable
        .isUnitary = cpu_allocHeapFlag(), // nullptr if failed
        .wasGpuSynced = cpu_allocHeapFlag(), // nullptr if failed

        // 1D CPU memory
        .cpuElems = cpu_allocArray(numElemsPerNode), // nullptr if failed

        // 1D GPU memory
        .gpuElems = (env.isGpuAccelerated)? gpu_allocArray(numElemsPerNode) : nullptr, // nullptr if failed or not needed
    };

    validateMatrixAllocs(out, __func__);
    setInitialHeapFlags(out);
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

    // optionally overwrite GPU elements with user-modified CPU elements
    if (mem_isAllocated(util_getGpuMemPtr(matr)))
        gpu_copyCpuToGpu(matr);

    // indicate that we do not know whether the revised matrix is
    // is unitarity; we defer establishing that until a unitarity check
    *(matr.isUnitary) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;

    // indicate that the matrix is now permanently GPU synchronised, even
    // if we are not in GPU-accelerated mode (in which case it's never consulted)
    *(matr.wasGpuSynced) = 1;
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

    freeHeapMatrix(matrix);
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

    // copy elems into matrix's CPU memory
    cpu_copyMatrix(out.cpuElems, elems, out.numRows);

    // overwrite GPU elements; validation gauranteed to pass
    syncCompMatr(out); 
}


extern "C" void setCompMatr(CompMatr out, qcomp** in) {

    // use the above template, which we will reuse for the vector overloads
    validateAndSetDenseMatrElems(out, in, __func__);
}


extern "C" void setDiagMatr(DiagMatr out, qcomp* in) {
    validate_matrixFields(out, __func__);

    // overwrite CPU memory
    cpu_copyArray(out.cpuElems, in, out.numElems);

    // overwrite GPU elements; validation gauranteed to pass
    syncDiagMatr(out);
}


extern "C" void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems) {
    validate_matrixFields(out, __func__);
    validate_fullStateDiagMatrNewElems(out, startInd, numElems, __func__);

    // if the matrix is non-distributed, we update every node's duplicated CPU amps
    if (!out.isDistributed)
        cpu_copyArray(&out.cpuElems[startInd], in, numElems);

    // only distributed nodes containing targeted elements need to do anything
    else if (util_areAnyElemsWithinThisNode(out.numElemsPerNode, startInd, numElems)) {
        util_IndexRange range = util_getLocalIndRangeOfElemsWithinThisNode(out.numElemsPerNode, startInd, numElems);
        cpu_copyArray(&out.cpuElems[range.localDistribStartInd], &in[range.localDuplicStartInd], range.numElems);
    }

    // all nodes overwrite GPU; validation gauranteed to succeed
    syncFullStateDiagMatr(out);
}



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA VECTORS
 */


extern "C" void validate_setCompMatrFromArr(CompMatr matr) {

    // we define bespoke validation used by setCompMatrFromArr(), which
    // is necessarily exposed in the header. We report the caller as setCompMatr
    // although the user may have actually called setInlineCompMatr()
    validate_matrixFields(matr, "setCompMatr");
}


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
