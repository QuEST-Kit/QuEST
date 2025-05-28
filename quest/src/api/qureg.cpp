/** @file
 * API definitions for creating and managing Quregs, 
 * and automatically choosing their deployment modes.
 * 
 * @author Tyson Jones
 */

#include "quest/include/qureg.h"
#include "quest/include/environment.h"
#include "quest/include/initialisations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <string>
#include <vector>

using std::string;
using std::vector;



/*
 * INTERNALLY EXPOSED FUNCTION
 */


Qureg qureg_populateNonHeapFields(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread) {

    QuESTEnv env = getQuESTEnv();

    // pre-prepare some struct fields (to avoid circular initialisation)
    int logNumNodes = (useDistrib)? 
        logBase2(env.numNodes) : 0;
    qindex logNumAmpsPerNode = (isDensMatr)? 
        (2*numQubits - logNumNodes) :
        (  numQubits - logNumNodes);

    // prepare output Qureg (avoiding C++20 designated initialiser)
    Qureg out;

    // bind deployment info
    out.isMultithreaded  = useMultithread;
    out.isGpuAccelerated = useGpuAccel;
    out.isDistributed    = useDistrib;

    // optionally bind distributed info, noting that in distributed environments,
    // the non-distributed quregs are duplicated on each node and each believe
    // they are the root node, with no other nodes existing; this is essential so
    // that these quregs can agnostically use distributed routines which consult
    // the rank, but it will interfere with naive root-only printing logic
    out.rank        = (useDistrib)? env.rank : 0;
    out.numNodes    = (useDistrib)? env.numNodes : 1;
    out.logNumNodes = (useDistrib)? logBase2(env.numNodes) : 0; // duplicated for clarity

    // bind dimensions
    out.isDensityMatrix = isDensMatr;
    out.numQubits  = numQubits;
    out.numAmps    = (isDensMatr)? powerOf2(2*numQubits) : powerOf2(numQubits);
    out.logNumAmps = (isDensMatr)?          2*numQubits  :          numQubits;

    // bind dimensions per node (even if not distributed)
    out.numAmpsPerNode = powerOf2(logNumAmpsPerNode);
    out.logNumAmpsPerNode = logNumAmpsPerNode;
    out.logNumColsPerNode = (isDensMatr)? numQubits - logNumNodes : 0; // used only by density matrices

    // caller will allocate heap memory as necessary
    out.cpuAmps       = nullptr;
    out.gpuAmps       = nullptr;
    out.cpuCommBuffer = nullptr;
    out.gpuCommBuffer = nullptr;

    return out;
}



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

    Qureg qureg = qureg_populateNonHeapFields(numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread);

    // always allocate CPU memory
    qureg.cpuAmps = cpu_allocArray(qureg.numAmpsPerNode); // nullptr if failed

    // conditionally allocate GPU memory and communication buffers (even if numNodes == 1).
    // note that in distributed settings but where useDistrib=false, each node will have a
    // full copy of the amplitudes, but will NOT have the communication buffers allocated.
    qureg.gpuAmps       = (useGpuAccel)?               gpu_allocArray(qureg.numAmpsPerNode) : nullptr;
    qureg.cpuCommBuffer = (useDistrib)?                cpu_allocArray(qureg.numAmpsPerNode) : nullptr;
    qureg.gpuCommBuffer = (useGpuAccel && useDistrib)? gpu_allocArray(qureg.numAmpsPerNode) : nullptr;

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

    using namespace printer_substrings;

    // 2^N = M
    string ampsStr;
    ampsStr  = bt + printer_toStr(qureg.numQubits * (qureg.isDensityMatrix? 2 : 1));
    ampsStr += eq + printer_toStr(qureg.numAmps);
    
    string colsStr = na;
    if (qureg.isDensityMatrix)
        colsStr = (
            bt + printer_toStr(qureg.numQubits) + 
            eq + printer_toStr(powerOf2(qureg.numQubits)));

    print_table(
        "dimension", {
        {"isDensMatr", printer_toStr(qureg.isDensityMatrix)},
        {"numQubits",  printer_toStr(qureg.numQubits)},
        {"numCols",    colsStr},
        {"numAmps",    ampsStr},
    });
}


void printDistributionInfo(Qureg qureg) {

    using namespace printer_substrings;

    // not applicable when not distributed
    string nodesStr = na;
    string ampsStr  = na;
    string colsStr  = na;

    // 2^N = M per node
    if (qureg.isDistributed) {
        nodesStr = bt + printer_toStr(qureg.logNumNodes)       + eq + printer_toStr(qureg.numNodes);
        ampsStr  = bt + printer_toStr(qureg.logNumAmpsPerNode) + eq + printer_toStr(qureg.numAmpsPerNode) + pn;
        if (qureg.isDensityMatrix)
            colsStr = bt + printer_toStr(qureg.logNumColsPerNode) + eq + printer_toStr(powerOf2(qureg.logNumColsPerNode)) + pn;
    }

    print_table(
        "distribution", {
        {"numNodes", nodesStr},
        {"numCols",  colsStr},
        {"numAmps",  ampsStr},
    });
}


void printMemoryInfo(Qureg qureg) {

    using namespace printer_substrings;

    size_t localArrayMem = mem_getLocalQuregMemoryRequired(qureg.numAmpsPerNode);
    string localMemStr = printer_getMemoryWithUnitStr(localArrayMem) + (qureg.isDistributed? pn : "");

    // precondition: no reportable fields are at risk of overflow as a qindex
    // type, EXCEPT aggregate total memory between distributed nodes (in bytes)
    qindex globalTotalMem = mem_getTotalGlobalMemoryUsed(qureg);
    string globalMemStr = (globalTotalMem == 0)? "overflowed" : printer_getMemoryWithUnitStr(globalTotalMem);

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
 * PUBLIC C & C++ AGNOSTIC FUNCTIONS
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


Qureg createForcedQureg(int numQubits) {
    validate_envIsInit(__func__);

    QuESTEnv env = getQuESTEnv();

    int isDensMatr = 0;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded, __func__);
}


Qureg createForcedDensityQureg(int numQubits) {
    validate_envIsInit(__func__);

    QuESTEnv env = getQuESTEnv();

    int isDensMatr = 1;
    return validateAndCreateCustomQureg(numQubits, isDensMatr, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded, __func__);
}


Qureg createCloneQureg(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // create a new Qureg with identical fields, but zero'd memory
    Qureg clone = validateAndCreateCustomQureg(
        qureg.numQubits, qureg.isDensityMatrix, qureg.isDistributed, 
        qureg.isGpuAccelerated, qureg.isMultithreaded, __func__);

    setQuregToClone(clone, qureg); // harmlessly re-validates

    // if GPU-accelerated, clone's CPU amps are NOT updated
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
    validate_numReportedNewlinesAboveZero(__func__); // because trailing newline mandatory

    /// @todo add function to write this output to file (useful for HPC debugging)

    // printer routines will consult env rank to avoid duplicate printing
    print_label("Qureg");
    printDeploymentInfo(qureg);
    printDimensionInfo(qureg);
    printDistributionInfo(qureg);
    printMemoryInfo(qureg);

    // exclude mandatory newline above
    print_oneFewerNewlines();
}


void reportQureg(Qureg qureg) {
    validate_quregFields(qureg, __func__);
    validate_numReportedNewlinesAboveZero(__func__); // because trailing newline mandatory

    // account all local CPU memory (including buffer), neglecting GPU memory
    // because it occupies distinct memory spaces, confusing accounting
    size_t localMem = mem_getLocalQuregMemoryRequired(qureg.numAmpsPerNode);
    if (qureg.isDistributed)
        localMem *= 2; // include buffer. @todo will this ever overflow?!?!
    
    // include struct size (expected negligibly tiny)
    localMem += sizeof(qureg);

    print_header(qureg, localMem);
    print_elems(qureg);

    // exclude mandatory newline above
    print_oneFewerNewlines();
}


void syncQuregToGpu(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // permit this to be called even in non-GPU mode
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg); // syncs then overwrites all local GPU amps
}
void syncQuregFromGpu(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // permit this to be called even in non-GPU mode
    if (qureg.isGpuAccelerated)
        gpu_copyGpuToCpu(qureg); // syncs then overwrites all local CPU amps
}


void syncSubQuregToGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps) {
    validate_quregFields(qureg, __func__);
    validate_localAmpIndices(qureg, localStartInd, numLocalAmps, __func__); 
    
    // the above validation communicates for node consensus in
    // distributed settings, because params can differ per-node.
    // note also this function accepts statevectors AND density
    // matrices, because the latter does not need a bespoke
    // (row,col) interface, because the user can only access/modify
    // local density matrix amps via a flat index anyway!

    // we permit this function to do nothing when not GPU-accelerated
    if (!qureg.isGpuAccelerated)
        return;

    // otherwise, every node merely copies its local subset, which
    // may differ per-node, in an embarrassingly parallel manner
    gpu_copyCpuToGpu(&qureg.cpuAmps[localStartInd], &qureg.gpuAmps[localStartInd], numLocalAmps);
}
void syncSubQuregFromGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps) {
    validate_quregFields(qureg, __func__);
    validate_localAmpIndices(qureg, localStartInd, numLocalAmps, __func__);

    // the above validation communicates for node consensus in
    // distributed settings, because params can differ per-node.
    // note also this function accepts statevectors AND density
    // matrices, because the latter does not need a bespoke
    // (row,col) interface, because the user can only access/modify
    // local density matrix amps via a flat index anyway!

    // we permit this function to do nothing when not GPU-accelerated
    if (!qureg.isGpuAccelerated)
        return;

    // otherwise, every node merely copies its local subset, which
    // may differ per-node, in an embarrassingly parallel manner
    gpu_copyGpuToCpu(&qureg.gpuAmps[localStartInd], &qureg.cpuAmps[localStartInd], numLocalAmps);
}


void getQuregAmps(qcomp* outAmps, Qureg qureg, qindex startInd, qindex numAmps) {
    validate_quregFields(qureg, __func__);
    validate_quregIsStateVector(qureg, __func__);
    validate_basisStateIndices(qureg, startInd, numAmps, __func__);

    localiser_statevec_getAmps(outAmps, qureg, startInd, numAmps);
}


void getDensityQuregAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_basisStateRowCols(qureg, startRow, startCol, numRows, numCols, __func__);

    localiser_densmatr_getAmps(outAmps, qureg, startRow, startCol, numRows, numCols);
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
    validate_quregFields(qureg, __func__);
    validate_quregIsStateVector(qureg, __func__);
    validate_basisStateIndex(qureg, index, __func__);

    return localiser_statevec_getAmp(qureg, index);
}
extern "C" void _wrap_getQuregAmp(qcomp* out, Qureg qureg, qindex index) {

    *out = getQuregAmp(qureg, index);
}


qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_basisStateRowCol(qureg, row, column, __func__);

    qindex ind = util_getGlobalFlatIndex(qureg, row, column);
    qcomp amp = localiser_statevec_getAmp(qureg, ind);
    return amp;
}
extern "C" void _wrap_getDensityQuregAmp(qcomp* out, Qureg qureg, qindex row, qindex column) {

    *out = getDensityQuregAmp(qureg, row, column);
}



/*
 * C++ OVERLOADS
 */


vector<qcomp> getQuregAmps(Qureg qureg, qindex startInd, qindex numAmps) {

    // allocate the output vector, and validate successful
    vector<qcomp> out;
    auto callback = [&]() { validate_tempAllocSucceeded(false, numAmps, sizeof(qcomp), __func__); };
    util_tryAllocVector(out, numAmps, callback);

    // performs main validation
    getQuregAmps(out.data(), qureg, startInd, numAmps);
    return out;
}


vector<vector<qcomp>> getDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols) {

    // allocate the output matrix, and validate successful
    vector<vector<qcomp>> out;
    qindex numElems = numRows * numCols; // never overflows (else Qureg alloc would fail)
    auto callback1 = [&]() { validate_tempAllocSucceeded(false, numElems, sizeof(qcomp), __func__); };
    util_tryAllocMatrix(out, numRows, numCols, callback1);

    // we must pass nested pointers to core C function, requiring another temp array, also validated
    vector<qcomp*> ptrs;
    auto callback2 = [&]() { validate_tempAllocSucceeded(false, numRows, sizeof(qcomp*), __func__); };
    util_tryAllocVector(ptrs, numRows, callback2);

    // embed out pointers
    for (qindex i=0; i<numRows; i++)
        ptrs[i] = out[i].data();

    // modify out through its ptrs
    getDensityQuregAmps(ptrs.data(), qureg, startRow, startCol, numRows, numCols);
    return out;
}
