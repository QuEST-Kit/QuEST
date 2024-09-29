/** @file
 * API definitions for creating and managing Quregs, and automatically
 * choosing their deployment modes.
 */

#include "quest/include/qureg.h"
#include "quest/include/environment.h"
#include "quest/include/initialisations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include "quest/src/core/errors.hpp" // only needed for function-not-implemented error

#include <string>

// provides substrings (by, na, pm, etc) used by reportQureg
using namespace printer_substrings;

using std::string;



/*
 * PRIVATE INNER FUNCTIONS (C++)
 */


bool didAnyLocalAllocsFail(Qureg qureg) {

    // CPU memory should always be allocated
    if (! mem_isAllocated(qureg.cpuAmps))
        return true;

    // when distributed, the CPU communication buffer must be allocated
    if (qureg.isDistributed && ! mem_isAllocated(qureg.cpuCommBuffer))
        return true;

    // when GPU-accelerated, the GPU memory should be allocated
    if (qureg.isGpuAccelerated && ! mem_isAllocated(qureg.gpuAmps))
        return true;

    // when both distributed and GPU-accelerated, the GPU communication buffer must be allocated
    if (qureg.isDistributed && qureg.isGpuAccelerated && ! mem_isAllocated(qureg.gpuCommBuffer))
        return true;

    // otherwise all pointers were non-NULL so no allocations failed
    return false;
}


bool didAnyAllocsFailOnAnyNode(Qureg qureg) {

    bool anyFail = didAnyLocalAllocsFail(qureg);
    if (comm_isInit())
        anyFail = comm_isTrueOnAllNodes(anyFail);

    return anyFail;
}


void freeAllMemoryIfAnyAllocsFailed(Qureg qureg) {

    // do nothing if everything allocated successfully between all nodes
    if (!didAnyAllocsFailOnAnyNode(qureg))
        return;

    // otherwise, free everything that was successfully allocated (freeing nullptr is legal)
    cpu_deallocArray(qureg.cpuAmps);
    cpu_deallocArray(qureg.cpuCommBuffer);

    // although we avoid calling GPU deallocation in non-GPU mode
    if (qureg.isGpuAccelerated) {
        gpu_deallocArray(qureg.gpuAmps);
        gpu_deallocArray(qureg.gpuCommBuffer);
    }
}


Qureg validateAndCreateCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread, const char* caller) {

    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    // ensure deployment is compatible with environment, considering available hardware and their memory capacities
    validate_newQuregParams(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env, caller);

    // automatically overwrite distrib, GPU, and multithread fields which were left as modeflag::USE_AUTO
    autodep_chooseQuregDeployment(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, env);

    // throw error if the user had forced multithreading but GPU accel was auto-chosen
    validate_newQuregNotBothMultithreadedAndGpuAccel(useGpuAccel, useMultithread, caller);

    // pre-prepare some struct fields (to avoid circular initialisation)
    int logNumNodes = (useDistrib)? 
        logBase2(env.numNodes) : 0;
    qindex logNumAmpsPerNode = (isDensMatr)? 
        (2*numQubits - logNumNodes) :
        (  numQubits - logNumNodes);

    // const struct fields must be initialised inline
    Qureg qureg = {

        // bind deployment info
        .isMultithreaded  = useMultithread,
        .isGpuAccelerated = useGpuAccel,
        .isDistributed    = useDistrib,

        // optionally bind distributed info, but always etain the env's rank because non-distributed
        // quregs are still duplicated between every node, and have duplicate processes
        .rank        = env.rank,
        .numNodes    = (useDistrib)? env.numNodes : 1,
        .logNumNodes = (useDistrib)? logBase2(env.numNodes) : 0, // duplicated for clarity

        // set dimensions
        .isDensityMatrix = isDensMatr,
        .numQubits  = numQubits,
        .numAmps    = (isDensMatr)? powerOf2(2*numQubits) : powerOf2(numQubits),
        .logNumAmps = (isDensMatr)?          2*numQubits  :          numQubits,

        // set dimensions per node (even if not distributed)
        .numAmpsPerNode = powerOf2(logNumAmpsPerNode),
        .logNumAmpsPerNode = logNumAmpsPerNode,
        .logNumColsPerNode = (isDensMatr)? numQubits - logNumNodes : 0, // used only by density matrices

        // always allocate CPU memory
        .cpuAmps = cpu_allocArray(qureg.numAmpsPerNode), // nullptr if failed

        // conditionally allocate GPU memory and communication buffers (even if numNodes == 1).
        // note that in distributed settings but where useDistrib=false, each node will have a
        // full copy of the amplitudes, but will NOT have the communication buffers allocated.
        .gpuAmps       = (useGpuAccel)?               gpu_allocArray(qureg.numAmpsPerNode) : nullptr,
        .cpuCommBuffer = (useDistrib)?                cpu_allocArray(qureg.numAmpsPerNode) : nullptr,
        .gpuCommBuffer = (useGpuAccel && useDistrib)? gpu_allocArray(qureg.numAmpsPerNode) : nullptr,
    };

    // if any of the above mallocs failed, below validation will memory leak; so free first (but don't set to nullptr)
    freeAllMemoryIfAnyAllocsFailed(qureg);
    validate_newQuregAllocs(qureg, __func__);

    // initialise state to |0> or |0><0|
    initZeroState(qureg); 

    return qureg;
}



/*
 * PRIVATE QUREG REPORTING INNER FUNCTIONS
 */


void printDeploymentInfo(Qureg qureg) {

    print_table(
        "deployment", {
        {"isMpiEnabled", qureg.isDistributed},
        {"isGpuEnabled", qureg.isGpuAccelerated},
        {"isOmpEnabled", qureg.isMultithreaded},
    });
}

void printDimensionInfo(Qureg qureg) {

    // 2^N = M or 2^N x 2^N = M
    string str;
    str  = bt + printer_toStr(qureg.numQubits);
    str += (qureg.isDensityMatrix)? mu + str : "";
    str += eq + printer_toStr(qureg.numAmps);

    print_table(
        "dimension", {
        {"isDensMatr", printer_toStr(qureg.isDensityMatrix)},
        {"numQubits",  printer_toStr(qureg.numQubits)},
        {"numAmps",    str},
    });
}


void printDistributionInfo(Qureg qureg) {

    // not applicable when not distributed
    string nodesStr = na;
    string ampsStr  = na;

    // 2^N = M per node
    if (qureg.isDistributed) {
        nodesStr = bt + printer_toStr(qureg.logNumNodes)       + eq + printer_toStr(qureg.numNodes);
        ampsStr  = bt + printer_toStr(qureg.logNumAmpsPerNode) + eq + printer_toStr(qureg.numAmpsPerNode);
        ampsStr += pn;
    }

    print_table(
        "distribution", {
        {"numNodes", nodesStr},
        {"numAmps",  ampsStr},
    });
}


void printMemoryInfo(Qureg qureg) {

    size_t localArrayMem = mem_getLocalQuregMemoryRequired(qureg.numAmpsPerNode);
    string localMemStr = printer_toStr(localArrayMem) + by + ((qureg.isDistributed)? pn : "");

    // precondition: no reportable fields are at risk of overflow as a qindex
    // type, EXCEPT aggregate total memory between distributed nodes (in bytes)
    qindex globalTotalMem = mem_getTotalGlobalMemoryUsed(qureg);
    string globalMemStr = (globalTotalMem == 0)? "overflowed" : (printer_toStr(globalTotalMem) + by);

    print_table(
        "memory", {
        {"cpuAmps",       mem_isAllocated(qureg.cpuAmps)?       localMemStr : na},
        {"gpuAmps",       mem_isAllocated(qureg.gpuAmps)?       localMemStr : na},
        {"cpuCommBuffer", mem_isAllocated(qureg.cpuCommBuffer)? localMemStr : na},
        {"gpuCommBuffer", mem_isAllocated(qureg.gpuCommBuffer)? localMemStr : na},
        {"globalTotal",   globalMemStr},
    });
}



/*
 * PUBLIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {


Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread) {

    return validateAndCreateCustomQureg(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread, __func__);
}


Qureg createQureg(int numQubits) {

    int isDensMatr = 0;
    int autoMode = modeflag::USE_AUTO;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, autoMode, autoMode, autoMode, __func__);
}


Qureg createDensityQureg(int numQubits) {

    int isDensMatr = 1;
    int autoMode = modeflag::USE_AUTO;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, autoMode, autoMode, autoMode, __func__);
}


Qureg createCloneQureg(Qureg qureg) {

    // create a new Qureg with identical fields, but zero'd memory
    Qureg clone = validateAndCreateCustomQureg(
        qureg.numQubits, qureg.isDensityMatrix, qureg.isDistributed, 
        qureg.isGpuAccelerated, qureg.isMultithreaded, __func__);

    // copy CPU memory from qureg to clone
    cpu_copyArray(clone.cpuAmps, qureg.cpuAmps, qureg.numAmpsPerNode);

    // update clone's GPU memory
    if (clone.isGpuAccelerated)
        gpu_copyCpuToGpu(clone);

    return clone;
}


void destroyQureg(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // free CPU memory
    cpu_deallocArray(qureg.cpuAmps);

    // free CPU communication buffer
    if (qureg.isDistributed)
        cpu_deallocArray(qureg.cpuCommBuffer);

    // free GPU memory
    if (qureg.isGpuAccelerated)
        gpu_deallocArray(qureg.gpuAmps);

    // free GPU communication buffer
    if (qureg.isGpuAccelerated && qureg.isDistributed)
        gpu_deallocArray(qureg.gpuCommBuffer);

    // cannot set free'd fields to nullptr because qureg
    // wasn't passed-by-reference, and isn't returned.
}


void reportQuregParams(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // TODO: add function to write this output to file (useful for HPC debugging)

    print("Qureg:");
    printDeploymentInfo(qureg);
    printDimensionInfo(qureg);
    printDistributionInfo(qureg);
    printMemoryInfo(qureg);
}


void reportQureg(Qureg qureg) {

    // TODO
    error_functionNotImplemented(__func__);
}

void reportQuregToFile(Qureg qureg, char* fn) {

    // TODO
    error_functionNotImplemented(__func__);
}


void syncQuregToGpu(Qureg qureg) {

    // TODO
    error_functionNotImplemented(__func__);
}

void syncQuregFromGpu(Qureg qureg) {

    // TODO
    error_functionNotImplemented(__func__);
}


void syncSubQuregToGpu  (Qureg qureg, qindex startInd, qindex numAmps) {

    // TODO
    error_functionNotImplemented(__func__);
}
void syncSubQuregFromGpu(Qureg qureg, qindex startInd, qindex numAmps) {

    // TODO
    error_functionNotImplemented(__func__);
}

void syncSubDensityQuregToGpu  (Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols) {

    // TODO
    error_functionNotImplemented(__func__);
}
void syncSubDensityQuregFromGpu(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols) {

    // TODO
    error_functionNotImplemented(__func__);
}





// end de-mangler
}




/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because of limited
 * interoperability of the qcomp type. See calculations.h 
 * for more info. We here define a C++-only signature (with
 * name-mangling), and a C-friendly wrapper which passes by
 * pointer; the C-friendly interface in wrappers.h which itself
 * wrap this.
 */


qcomp getQuregAmp(Qureg qureg, qindex index) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_getQuregAmp(qcomp* out, Qureg qureg, qindex index) {

    *out = getQuregAmp(qureg, index);
}


qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_getDensityQuregAmp(qcomp* out, Qureg qureg, qindex row, qindex column) {

    *out = getDensityQuregAmp(qureg, row, column);
}
