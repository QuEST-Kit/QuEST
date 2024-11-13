/** @file
 * Data structures for representing arbitrary channels as
 * superoperators and Kraus maps, including their constructors, 
 * getters, setters and reporters. Note the functions to
 * actually simulate these channels are exposed in decoherence.h
 */

#include "quest/include/channels.h"
#include "quest/include/types.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <vector>
#include <cstdlib>

using std::vector;



/*
 * PRIVATE UTILITIES
 */


void freeSuperOp(SuperOp op) {

    // free CPU memory, even if it is NULL
    cpu_deallocMatrix(op.cpuElems, op.numRows);

    // free teeniy-tiny heap flag
    cpu_deallocHeapFlag(op.wasGpuSynced);

    // we avoid invoking a GPU function in non-GPU mode
    auto gpuPtr = util_getGpuMemPtr(op);
    if (mem_isAllocated(gpuPtr))
        gpu_deallocArray(gpuPtr);
}


void freeKrausMap(KrausMap map) {

    // free CPU matrices of Kraus operators
    cpu_deallocMatrixList(map.matrices, map.numMatrices, map.numRows);

    // free teeny-tiny heap flag
    cpu_deallocHeapFlag(map.isCPTP);

    // free superoperator (which may include freeing GPU memory)
    freeSuperOp(map.superop);
}


void freeObj(SuperOp op) {
    freeSuperOp(op);
}
void freeObj(KrausMap map) {
    freeKrausMap(map);
}


bool didAnyLocalAllocsFail(SuperOp op) {

    // god help us if this single-integer malloc failed
    if (!mem_isAllocated(op.wasGpuSynced))
        return true;

    if (!mem_isAllocated(op.cpuElems, op.numRows))
        return true;

    if (getQuESTEnv().isGpuAccelerated && !mem_isAllocated(op.gpuElemsFlat))
        return true;
    
    return false;
}


bool didAnyLocalAllocsFail(KrausMap map) {

    // god help us if this single-integer malloc failed
    if (!mem_isAllocated(map.isCPTP))
        return true;

    // list of CPU matrices and all matrices/rows therein shoul dbe non-NULL
    if (!mem_isAllocated(map.matrices, map.numMatrices, map.numRows))
        return true;

    // check if anything in the superoperator failed to allocate
    if (didAnyLocalAllocsFail(map.superop))
        return true;

    // otherwise, all pointers were non-NULL and ergo all allocs were successful
    return false;
}


// T will be SuperOp or KrausMap
template <typename T>
void freeAllMemoryIfAnyAllocsFailed(T obj) {

    // determine whether any node experienced a failure
    bool anyFail = didAnyLocalAllocsFail(obj);
    if (comm_isInit())
        anyFail = comm_isTrueOnAllNodes(anyFail);

    // if so, free all memory before subsequent validation
    if (anyFail)
        freeObj(obj);
}



/*
 * CONSTRUCTORS
 */


SuperOp allocSuperOp(int numQubits) {

    // this inner function will NOT validate, verify whether 
    // allocations succeeded, nor attempt to free memory if not

    // prior validation ensures these never overflow
    qindex numRows = powerOf2(2 * numQubits);
    qindex numElems = numRows * numRows;

    SuperOp out = {
        .numQubits = numQubits,
        .numRows = numRows,

        .cpuElems = cpu_allocMatrix(numRows), // nullptr if failed
        .gpuElemsFlat = (getQuESTEnv().isGpuAccelerated)? gpu_allocArray(numElems) : nullptr, // nullptr if failed or not needed

        .wasGpuSynced = cpu_allocHeapFlag() // nullptr if failed
    };

    // if heap flag allocated, mark it as unsynced (caller will handle if allocation failed)
    if (mem_isAllocated(out.wasGpuSynced))
        *(out.wasGpuSynced) = 0;

    return out;
}


extern "C" SuperOp createSuperOp(int numQubits) {
    validate_envIsInit(__func__);
    validate_newSuperOpParams(numQubits, __func__);

    SuperOp out = allocSuperOp(numQubits); // fields may be or contain nullptr if failed

    // free all memory before we throw validation errors to avoid leaks
    freeAllMemoryIfAnyAllocsFailed(out);
    validate_newSuperOpAllocs(out, __func__);

    return out;
}


extern "C" KrausMap createKrausMap(int numQubits, int numOperators) {
    validate_envIsInit(__func__);
    validate_newKrausMapParams(numQubits, numOperators, __func__);

    // validation ensures this never overflows
    qindex numRows = powerOf2(numQubits);

    KrausMap out = {
        .numQubits = numQubits,
        .numMatrices = numOperators,
        .numRows = numRows,

        .matrices = cpu_allocMatrixList(numOperators, numRows), // is or contains nullptr if failed
        .superop = allocSuperOp(numQubits), // heap fields are or contain nullptr if failed

        .isCPTP = cpu_allocHeapFlag(), // nullptr if failed
    };

    // free memory before throwing validation error to avoid memory leaks
    freeAllMemoryIfAnyAllocsFailed(out);
    validate_newKrausMapAllocs(out, __func__);

    // mark CPTP as unknown; it will be lazily evaluated whene a function asserts CPTP-ness
    *(out.isCPTP) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;

    return out;
}


extern "C" void destroySuperOp(SuperOp op) {
    validate_superOpFields(op, __func__);

    freeSuperOp(op);
}


extern "C" void destroyKrausMap(KrausMap map) {
    validate_krausMapFields(map, __func__);

    freeKrausMap(map);
}



/*
 * SYNC
 */


extern "C" void syncSuperOp(SuperOp op) {
    validate_superOpFields(op, __func__);

    // optionally overwrite GPU elements with user-modified CPU elements
    if (mem_isAllocated(util_getGpuMemPtr(op)))
        gpu_copyCpuToGpu(op);

    // indicate that the matrix is now permanently GPU synchronised, even
    // if we are not in GPU-accelerated mode - it hardly matters
    *(op.wasGpuSynced) = 1;
}


extern "C" void syncKrausMap(KrausMap map) {
    validate_krausMapFields(map, __func__);

    // re-compute and GPU-sync the superoperator from the modified matrices
    util_setSuperoperator(map.superop.cpuElems, map.matrices, map.numMatrices, map.numQubits);
    syncSuperOp(map.superop);

    // indicate that we do not know whether the revised map is
    // is CPTP; we defer establishing that until a CPTP check
    *(map.isCPTP) = validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
}



/*
 * SUPEROPERATOR SETTERS
 *
 * defining the C and C++ compatible pointer version,
 * and the C++-only vector overloads. Note the C-only
 * VLA overloads are defined in the header
 */


// T can be qcomp** or vector<vector<qcomp>>
template <typename T> 
void setAndSyncSuperOpElems(SuperOp op, T matrix) {
    
    // overwrite the CPU matrix
    cpu_copyMatrix(op.cpuElems, matrix, op.numRows);

    // sync CPU matrix to flat GPU array
    syncSuperOp(op);
}


extern "C" void setSuperOp(SuperOp op, qcomp** matrix) {
    validate_superOpFields(op, __func__);

    setAndSyncSuperOpElems(op, matrix);
}

void setSuperOp(SuperOp op, vector<vector<qcomp>> matrix) {
    validate_superOpFields(op, __func__);
    validate_superOpNewMatrixDims(op, matrix, __func__);

    setAndSyncSuperOpElems(op, matrix);
}



/*
 * KRAUS MAP SETTERS
 *
 * defining the C and C++ compatible pointer version,
 * and the C++-only vector overloads. Note the C-only
 * VLA overloads are defined in the header
 */


// type T can be qcomp*** or vector<vector<vector<qcomp>>>
template <typename T> 
void setAndSyncKrausMapElems(KrausMap map, T matrices) {

    // copy over the given matrices into the map's CPU mmemory
    for (int n=0; n<map.numMatrices; n++)
        cpu_copyMatrix(map.matrices[n], matrices[n], map.numRows);

    // update the superoperator, including its GPU memory, and its isCPTP flag
    syncKrausMap(map);
}


extern "C" void setKrausMap(KrausMap map, qcomp*** matrices) {
    validate_krausMapFields(map, __func__);

    setAndSyncKrausMapElems(map, matrices);
}


void setKrausMap(KrausMap map, vector<vector<vector<qcomp>>> matrices) {
    validate_krausMapFields(map, __func__);
    validate_krausMapNewMatrixDims(map, matrices, __func__);

    setAndSyncKrausMapElems(map, matrices);
}



/*
 * LITERAL SETTERS
 *
 * Only the C++ versions are defined here, because the C versions are macros
 * defined in the header. Note the C++ versions themselves are entirely
 * superfluous and merely call the above vector setters, but we still define
 * them for API consistency, and we additionally validate the superfluous
 * additional parameters they pass.
 */


void setInlineKrausMap(KrausMap map, int numQb, int numOps, vector<vector<vector<qcomp>>> matrices) {
    validate_krausMapFields(map, __func__);
    validate_krausMapFieldsMatchPassedParams(map, numQb, numOps, __func__);
    validate_krausMapNewMatrixDims(map, matrices, __func__);

    setAndSyncKrausMapElems(map, matrices);
}


void setInlineSuperOp(SuperOp op, int numQb, vector<vector<qcomp>> matrix) {
    validate_superOpFields(op, __func__);
    validate_superOpFieldsMatchPassedParams(op, numQb, __func__);
    validate_superOpNewMatrixDims(op, matrix, __func__);

    setAndSyncSuperOpElems(op, matrix);
}



/*
 * LITERAL CREATORS
 *
 * Only the C++ versions are defined here; the C versions are in-header macros.
 */


KrausMap createInlineKrausMap(int numQubits, int numOperators, vector<vector<vector<qcomp>>> matrices) {
    validate_envIsInit(__func__);
    validate_newKrausMapParams(numQubits, numOperators, __func__);
    validate_newInlineKrausMapDimMatchesVectors(numQubits, numOperators, matrices, __func__);

    // pre-validation gauranteed to pass, but malloc failures will trigger an error 
    // message specific to 'createKrausMap', rather than this 'inline' version. Alas!
    KrausMap map = createKrausMap(numQubits, numOperators);
    setAndSyncKrausMapElems(map, matrices);
    return map;
}


SuperOp createInlineSuperOp(int numQubits, vector<vector<qcomp>> matrix) {
    validate_envIsInit(__func__);
    validate_newSuperOpParams(numQubits, __func__);
    validate_newInlineSuperOpDimMatchesVectors(numQubits, matrix, __func__);

    // pre-validation gauranteed to pass, but malloc failures will trigger an error 
    // message specific to 'createSuperOp', rather than this 'inline' version. Alas!
    SuperOp op = createSuperOp(numQubits);
    setAndSyncSuperOpElems(op, matrix);
    return op;
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

    void _validateParamsToSetKrausMapFromArr(KrausMap map) {
        validate_krausMapFields(map, "setKrausMap");
    }

    void _validateParamsToSetSuperOpFromArr(SuperOp op) {
        validate_superOpFields(op, "setSuperOp");
    }

    void _validateParamsToSetInlineKrausMap(KrausMap map, int numQb, int numOps) {

        const char* caller = "setInlineKrausMap";
        validate_krausMapFields(map, caller);
        validate_krausMapFieldsMatchPassedParams(map, numQb, numOps, caller);
    }

    void _validateParamsToSetInlineSuperOp(SuperOp op, int numQb) {

        const char* caller = "setInlineSuperOp";
        validate_superOpFields(op, caller);
        validate_superOpFieldsMatchPassedParams(op, numQb, caller);
    }

    void _validateParamsToCreateInlineKrausMap(int numQb, int numOps) {

        const char* caller = "createInlineKrausMap";
        validate_envIsInit(caller);
        validate_newKrausMapParams(numQb, numOps, caller);
    }

    void _validateParamsToCreateInlineSuperOp(int numQb) {

        const char* caller = "createInlineSuperOp";
        validate_envIsInit(caller);
        validate_newSuperOpParams(numQb, caller);
    }

}



/*
 * REPORTING
 */


extern "C" void reportSuperOp(SuperOp op) {
    validate_superOpFields(op, __func__);

    // demand that SuperOp is GPU-synced, since the
    // to-be-printed GPU elems will overwrite CPU
    validate_superOpIsSynced(op, __func__);

    // determine memory costs; heap memory and struct fields
    size_t elemMem = mem_getLocalSuperOpMemoryRequired(op.numQubits);
    size_t structMem = sizeof(op);

    print_header(op, elemMem + structMem);
    print_elems(op);
}


extern "C" void reportKrausMap(KrausMap map) {
    validate_krausMapFields(map, __func__);

    // pedantically demand that the map's SuperOP is GPU-synced,
    // even though only the CPU Kraus operators are printed
    validate_krausMapIsSynced(map, __func__);

    // determine memory costs (gauranteed not to overflow)
    bool isDense = true;
    int numNodes = 1;
    size_t krausMem = mem_getLocalMatrixMemoryRequired(map.numQubits, isDense, numNodes) * map.numMatrices;
    size_t superMem = mem_getLocalSuperOpMemoryRequired(map.superop.numQubits);
    size_t strucMem = sizeof(map);

    // gauranteed not to overflow
    size_t totalMem = krausMem + superMem + strucMem;
    print_header(map, totalMem);
    print_elems(map);
}
