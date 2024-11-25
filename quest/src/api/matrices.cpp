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
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions

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
    cpu_deallocHeapFlag(matr.isHermitian);
    cpu_deallocHeapFlag(matr.wasGpuSynced);
}


// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
bool didAnyLocalAllocsFail(T matr) {

    // god help us if these single-integer malloc failed
    if (!mem_isAllocated(matr.isUnitary))    return true;
    if (!mem_isAllocated(matr.isHermitian))  return true;
    if (!mem_isAllocated(matr.wasGpuSynced)) return true;

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

    // if GPU memory is not allocated in a GPU environment...
    bool isGpuAlloc = mem_isAllocated(util_getGpuMemPtr(matr));
    if (getQuESTEnv().isGpuAccelerated && !isGpuAlloc) {

        // then FullStateDiagMatr GPU alloc failed only if it tried...
        if constexpr (util_isFullStateDiagMatr<T>()) {
            if (matr.isGpuAccelerated)
                return true;
            
        // but all other matrices always try to alloc, so must have failed
        } else
            return true;
    }

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

    // set initial propreties of the newly created matrix to unknown
    *(matr.isUnitary)   = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
    *(matr.isHermitian) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;

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
        .isUnitary    = cpu_allocHeapFlag(), // nullptr if failed
        .isHermitian  = cpu_allocHeapFlag(), // nullptr if failed
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
        .isUnitary    = cpu_allocHeapFlag(), // nullptr if failed
        .isHermitian  = cpu_allocHeapFlag(), // nullptr if failed
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


FullStateDiagMatr validateAndCreateCustomFullStateDiagMatr(int numQubits, int useDistrib, int useGpuAccel, const char* caller) {
    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    // validate parameters before passing them to autodeployer
    validate_newFullStateDiagMatrParams(numQubits, useDistrib, useGpuAccel, caller);

    // overwrite useDistrib and useGpuAccel if they were left as AUTO_FLAG
    autodep_chooseFullStateDiagMatrDeployment(numQubits, useDistrib, useGpuAccel, env);

    // validation ensures this never overflows
    qindex numElems = powerOf2(numQubits);
    qindex numElemsPerNode = numElems / (useDistrib? env.numNodes : 1); // divides evenly

    FullStateDiagMatr out = {

        .numQubits = numQubits,
        .numElems = numElems,

        // data deployment configuration; disable distrib if deployed to 1 node
        .isGpuAccelerated = useGpuAccel,
        .isDistributed = useDistrib && (env.numNodes > 1),
        .numElemsPerNode = numElemsPerNode,

        // allocate flags in the heap so that struct copies are mutable
        .isUnitary    = cpu_allocHeapFlag(), // nullptr if failed
        .isHermitian  = cpu_allocHeapFlag(), // nullptr if failed
        .wasGpuSynced = cpu_allocHeapFlag(), // nullptr if failed

        // 1D CPU memory
        .cpuElems = cpu_allocArray(numElemsPerNode), // nullptr if failed

        // 1D GPU memory
        .gpuElems = (useGpuAccel)? gpu_allocArray(numElemsPerNode) : nullptr, // nullptr if failed or not needed
    };

    validateMatrixAllocs(out, __func__);
    setInitialHeapFlags(out);
    return out;
}

extern "C" FullStateDiagMatr createCustomFullStateDiagMatr(int numQubits, int useDistrib, int useGpuAccel) {

    return validateAndCreateCustomFullStateDiagMatr(numQubits, useDistrib, useGpuAccel, __func__);
}

extern "C" FullStateDiagMatr createFullStateDiagMatr(int numQubits) {

    return validateAndCreateCustomFullStateDiagMatr(numQubits, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);
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

    // indicate that we do not know the revised matrix properties;
    // we defer establishing that until validation needs to check them
    *(matr.isUnitary)   = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
    *(matr.isHermitian) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;

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
void setAndSyncDenseMatrElems(CompMatr out, T elems) {
    
    // copy elems into matrix's CPU memory
    cpu_copyMatrix(out.cpuElems, elems, out.numRows);

    // overwrite GPU elements; validation gauranteed to pass
    syncCompMatr(out); 
}


extern "C" void setCompMatr(CompMatr out, qcomp** in) {
    validate_matrixFields(out, __func__);

    setAndSyncDenseMatrElems(out, in);
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

    // overwrites both the CPU and GPU memory (if it exists), maintaining consistency
    localiser_fullstatediagmatr_setElems(out, startInd, in, numElems);
}



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA VECTORS
 */


void setCompMatr(CompMatr out, vector<vector<qcomp>> in) {
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__);

    setAndSyncDenseMatrElems(out, in);
}


void setDiagMatr(DiagMatr out, vector<qcomp> in) {
    validate_matrixFields(out, __func__);
    validate_matrixNumNewElems(out.numQubits, in, __func__); // validates 'in' dim

    setDiagMatr(out, in.data()); // harmessly re-validates
}


void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, vector<qcomp> in) {

    setFullStateDiagMatr(out, startInd, in.data(), in.size());
}


// no bespoke array functions are necessary for diagonal matrices initialisation, 
// since passed arrays automatically decay to pointers



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA LITERALS
 *
 * Only the C++ versions are defined here, while the C versions are macros
 * defined in the header. Note the C++ versions themselves are entirely
 * superfluous and merely call the above vector setters, but we still define
 * them for API consistency, and we additionally validate the superfluous
 * additional parameters they pass.
 */


void setInlineCompMatr(CompMatr matr, int numQb, vector<vector<qcomp>> in) {
    validate_matrixFields(matr, __func__);
    validate_matrixNumQubitsMatchesParam(matr.numQubits, numQb, __func__);
    validate_matrixNumNewElems(matr.numQubits, in, __func__);

    setAndSyncDenseMatrElems(matr, in);
}

void setInlineDiagMatr(DiagMatr matr, int numQb, vector<qcomp> in) {
    validate_matrixFields(matr, __func__);
    validate_matrixNumQubitsMatchesParam(matr.numQubits, numQb, __func__);
    validate_matrixNumNewElems(matr.numQubits, in, __func__);

    setDiagMatr(matr, in.data()); // validation gauranteed to pass
}

void setInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems, vector<qcomp> in) {
    validate_matrixFields(matr, __func__);
    validate_declaredNumElemsMatchesVectorLength(numElems, in.size(), __func__);
    validate_fullStateDiagMatrNewElems(matr, startInd, numElems, __func__);

    setFullStateDiagMatr(matr, startInd, in); // validation gauranteed to pass
}



/*
 * VARIABLE-SIZE MATRIX INLINE-SETTER CONSTRUCTORS 
 *
 * Only the C++ versions are defined here; the C versions are header macros
 */


CompMatr createInlineCompMatr(int numQb, vector<vector<qcomp>> elems) {
    validate_envIsInit(__func__);
    validate_newCompMatrParams(numQb, __func__);
    validate_matrixNumNewElems(numQb, elems, __func__);

    // pre-validation gauranteed to pass, but malloc failures will trigger an error 
    // message specific to 'createCompMatr', rather than this 'inline' version. Alas!
    CompMatr matr = createCompMatr(numQb);
    setAndSyncDenseMatrElems(matr, elems);
    return matr;
}


DiagMatr createInlineDiagMatr(int numQb, vector<qcomp> elems) {
    validate_envIsInit(__func__);
    validate_newDiagMatrParams(numQb, __func__);
    validate_matrixNumNewElems(numQb, elems, __func__);

    // pre-validation gauranteed to pass, but malloc failures will trigger an error 
    // message specific to 'createCompMatr', rather than this 'inline' version. Alas!
    DiagMatr matr = createDiagMatr(numQb);
    setDiagMatr(matr, elems.data()); // validation gauranteed to pass
    return matr;
}



/*
 * EXPOSING SOME SETTER VALIDATION TO HEADER
 *
 * Some setters are necessarily defined in the header, because they accept 
 * C-only VLAs which need to be cast into pointers before being passed to 
 * this C++ backend (which does not support VLA). These setters need their
 * validators exposed, though we cannot expose the entirety of validation.hpp 
 * because it cannot be parsed by C; so we here wrap specific functions.
 */


extern "C" {

    void _validateParamsToSetCompMatrFromArr(CompMatr matr) { 

        validate_matrixFields(matr, "setCompMatr");
    }

    void _validateParamsToSetInlineCompMatr(CompMatr matr, int numQb) {

        const char* caller = "setInlineCompMatr";
        validate_matrixFields(matr, caller);
        validate_matrixNumQubitsMatchesParam(matr.numQubits, numQb, caller);
    }

    void _validateParamsToSetInlineDiagMatr(DiagMatr matr, int numQb) {

        const char* caller = "setInlineDiagMatr";
        validate_matrixFields(matr, caller);
        validate_matrixNumQubitsMatchesParam(matr.numQubits, numQb, caller);
    }

    void _validateParamsToSetInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems) {

        const char* caller = "setInlineFullStateDiagMatr";
        validate_matrixFields(matr, caller);
        validate_fullStateDiagMatrNewElems(matr, startInd, numElems, caller);
    }

    void _validateParamsToCreateInlineCompMatr(int numQb) {

        const char* caller = "createInlineCompMatr";
        validate_envIsInit(caller);
        validate_newCompMatrParams(numQb, caller);
    }

    void _validateParamsToCreateInlineDiagMatr(int numQb) {

        const char* caller = "createInlineDiagMatr";
        validate_envIsInit(caller);
        validate_newDiagMatrParams(numQb, caller);
    }

}



/*
 * SPECIAL CREATORS AND SETTERS
 */


extern "C" void setFullStateDiagMatrFromPauliStrSum(FullStateDiagMatr out, PauliStrSum in) {

    // TODO
    error_functionNotImplemented(__func__);
}


extern "C" FullStateDiagMatr createFullStateDiagMatrFromPauliStrSumFile(char* fn) {

    // TODO
    error_functionNotImplemented(__func__);
    return createFullStateDiagMatr(-1);
}


vector<qindex> getMultiVarValuesFromIndex(qindex index, int* numQubitsPerVar, int numVars, bool areSigned) {

    vector<qindex> varValues(numVars);
    qindex remainingValue = index;

    // find unsigned integer value of each var by partitioning index bits between vars
    for (int v=0; v<numVars; v++) {
        varValues[v] =  getBitsRightOfIndex(remainingValue, numQubitsPerVar[v]);
        remainingValue = getBitsLeftOfIndex(remainingValue, numQubitsPerVar[v]-1);
    }

    // two's-complement signed integers are negated if the leftmost variable sign bit is 1 
    if (areSigned)
        for (int v=0; v<numVars; v++)
            if (getBit(varValues[v], numQubitsPerVar[v]-1) == 1)
                varValues[v] -= powerOf2(numQubitsPerVar[v] - 1);

    return varValues;
}


extern "C" void setDiagMatrFromMultiVarFunc(DiagMatr out, qcomp (*callbackFunc)(qindex*), int* numQubitsPerVar, int numVars, int areSigned) {
    validate_matrixFields(out, __func__);
    validate_multiVarFuncQubits(out.numQubits, numQubitsPerVar, numVars, __func__);
    validate_funcVarSignedFlag(areSigned, __func__);

    // set each element of the diagonal in-turn; user's callback might not be thread-safe
    for (qindex elemInd=0; elemInd<out.numElems; elemInd++) {
        vector<qindex> varValues = getMultiVarValuesFromIndex(elemInd, numQubitsPerVar, numVars, areSigned);

        // call user function, and update only the CPU elems
        out.cpuElems[elemInd] = callbackFunc(varValues.data());
    }

    // overwrite all GPU elems
    syncDiagMatr(out);
}


extern "C" void setFullStateDiagMatrFromMultiVarFunc(FullStateDiagMatr out, qcomp (*callbackFunc)(qindex*), int* numQubitsPerVar, int numVars, int areSigned) {
    validate_matrixFields(out, __func__);
    validate_multiVarFuncQubits(out.numQubits, numQubitsPerVar, numVars, __func__);
    validate_funcVarSignedFlag(areSigned, __func__);

    int rank = getQuESTEnv().rank;
    int numLocalIndBits = logBase2(out.numElemsPerNode);

    // every node updates their local elems embarrassingly parallel, but user's callback might not be thread-safe
    for (qindex localInd=0; localInd<out.numElemsPerNode; localInd++) {

        // each local index corresponds to a unique global index which informs variable values
        qindex globalInd = concatenateBits(rank, localInd, numLocalIndBits);
        vector<qindex> varValues = getMultiVarValuesFromIndex(globalInd, numQubitsPerVar, numVars, areSigned);

        // call user function, and update only the CPU elems
        out.cpuElems[localInd] = callbackFunc(varValues.data());
    }

    // overwrite all GPU elems
    syncFullStateDiagMatr(out);
}


qcomp getElemFromNestedPtrs(void* in, qindex* inds, int numInds) {
    qindex ind = inds[0];

    if (numInds == 1)
        return ((qcomp*) in)[ind];

    qcomp* ptr = ((qcomp**) in)[ind];
    return getElemFromNestedPtrs(ptr, &inds[1], numInds-1); // compiler may optimise tail-recursion
}


extern "C" void setDiagMatrFromMultiDimLists(DiagMatr out, void* lists, int* numQubitsPerDim, int numDims) {
    validate_matrixFields(out, __func__);
    validate_multiVarFuncQubits(out.numQubits, numQubitsPerDim, numDims, __func__);

    // set each element of the diagonal in-turn, which is embarrassingly parallel, although we do not parallelise
    for (qindex elemInd=0; elemInd<out.numElems; elemInd++) {

        // nested list indices = unsigned integer values of variables
        vector<qindex> listInds = getMultiVarValuesFromIndex(elemInd, numQubitsPerDim, numDims, false);

        // update only the CPU elems
        out.cpuElems[elemInd] = getElemFromNestedPtrs(lists, listInds.data(), numDims);
    }

    // overwrite all GPU elems
    syncDiagMatr(out);
}


extern "C" void setFullStateDiagMatrFromMultiDimLists(FullStateDiagMatr out, void* lists, int* numQubitsPerDim, int numDims) {
    validate_matrixFields(out, __func__);
    validate_multiVarFuncQubits(out.numQubits, numQubitsPerDim, numDims, __func__);

    int rank = getQuESTEnv().rank;
    int numLocalIndBits = logBase2(out.numElemsPerNode);

    // every node updates their local elems embarrassingly parallel, but user's callback might not be thread-safe
    for (qindex localInd=0; localInd<out.numElemsPerNode; localInd++) {

        // each local diag index corresponds to a unique global index which informs list indices
        qindex globalInd = concatenateBits(rank, localInd, numLocalIndBits);
        vector<qindex> listInds = getMultiVarValuesFromIndex(globalInd, numQubitsPerDim, numDims, false);

        // update only the CPU elems using lists which are duplicated on every node
        out.cpuElems[localInd] = getElemFromNestedPtrs(lists, listInds.data(), numDims);
    }

    // overwrite all GPU elems
    syncFullStateDiagMatr(out);
}



/*
 * MATRIX REPORTERS
 */


// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template<class T> 
void validateAndPrintMatrix(T matr, const char* caller) {
    validate_matrixFields(matr, caller);

    // syncable matrices must be synced before reporting (though only CPU elems are printed)
    if constexpr (util_isHeapMatrixType<T>())
        validate_matrixIsSynced(matr, caller);

    // calculate the total memory (in bytes) consumed by the matrix on each
    // node, which will depend on whether the matrix is distributed, and
    // includes the size of the matrix struct itself. Note that GPU memory
    // is not included since it is less than or equal to the CPU memory, and
    // occupies different memory spaces, confusing capacity accounting
    int numNodes = (util_isDistributedMatrix(matr))? comm_getNumNodes() : 1;
    size_t elemMem = mem_getLocalMatrixMemoryRequired(matr.numQubits, util_isDenseMatrixType<T>(), numNodes);
    size_t structMem = sizeof(matr);

    // struct memory includes fixed-size arrays (qcomp[][]), so we undo double-counting
    if (util_isFixedSizeMatrixType<T>())
        structMem -= elemMem;

    size_t numBytesPerNode = elemMem + structMem;
    print_header(matr, numBytesPerNode);
    print_elems(matr);
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
