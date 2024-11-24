/** @file
 * Internal functions which localize the data needed for simulation.
 * That is, they determine whether performing a simulation requires
 * Qureg amplitudes from other distributed nodes and if so, invoke
 * the necessary communication, before finally calling the 
 * embarrassingly parallel subroutines in accelerator.cpp. This is
 * done agnostically of whether amplitudes of the Qureg are being
 * stored in RAM (CPU) or VRAM (GPU).
 */

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/initialisations.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <tuple>
#include <array>
#include <vector>
#include <complex>
#include <algorithm>

using std::vector;



/*
 * PRIVATE FUNCTIONS
 */


// some localiser functions cheekily spoof Quregs for code-reuse
extern Qureg qureg_populateNonHeapFields(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread);


void assertValidCtrlStates(vector<int> ctrls, vector<int> ctrlStates) {

    // providing no control states is always valid (to invoke default all-on-1)
    if (ctrlStates.empty())
        return;

    // otherwise a state must be explicitly given for each ctrl
    if (ctrlStates.size() != ctrls.size())
        error_localiserNumCtrlStatesInconsistentWithNumCtrls();
}


void setDefaultCtrlStates(vector<int> ctrls, vector<int> &states) {

    // no states necessary if there are no control qubits
    if (ctrls.empty())
        return;

    // default ctrl state is all-1
    if (states.empty())
        states.insert(states.end(), ctrls.size(), 1);
}


bool doesGateRequireComm(Qureg qureg, vector<int> targs) {

    // non-distributed quregs never communicate (duh)
    if (!qureg.isDistributed)
        return false;

    // communication necessary when any prefix qubit is targeted
    for (int targ : targs)
        if (!util_isQubitInSuffix(targ, qureg))
            return true;

    // sufix qubit targets need no communication
    return false;
}

bool doesGateRequireComm(Qureg qureg, int targ) {

    return doesGateRequireComm(qureg, vector{targ});
}


bool doesChannelRequireComm(Qureg qureg, vector<int> ketQubits) {
    if (!qureg.isDensityMatrix)
        error_localiserPassedStateVecToChannelComCheck();

    // ket-qubits are gauranteed to be in the suffix (because we distributed >=1 column per node),
    // so channels invoke communication if any corresponding bra-qubits are in prefix
    auto braQubits = util_getBraQubits(ketQubits, qureg);
    return doesGateRequireComm(qureg, braQubits);
}

bool doesChannelRequireComm(Qureg qureg, int ketQubit) {

    return doesChannelRequireComm(qureg, vector{ketQubit});
}


bool doAnyLocalStatesHaveQubitValues(Qureg qureg, vector<int> qubits, vector<int> states) {

    // this answers the generic question of "do any of the given qubits lie in the
    // prefix substate with node-fixed values inconsistent with the given states?"

    // non-distributed quregs always have amps satisfying ctrls
    if (!qureg.isDistributed)
        return true;

    // check each ctrl qubit
    for (size_t i=0; i<qubits.size(); i++) {

        // consider only ctrls which operate on the prefix substate
        if (util_isQubitInSuffix(qubits[i], qureg))
            continue;

        // abort if any prefix ctrl has wrong bit value
        if (util_getRankBitOfQubit(qubits[i], qureg) != states[i])
            return false;
    }

    // otherwise all prefix qubits have the specified values
    return true;
}


void removePrefixQubitsAndStates(Qureg qureg, vector<int> &qubits, vector<int> &states) {

    vector<int> suffixQubits(0);  suffixQubits.reserve(qubits.size());
    vector<int> suffixStates(0);  suffixStates.reserve(states.size());

    // collect suffix qubits/states
    for (size_t i=0; i<qubits.size(); i++)
        if (util_isQubitInSuffix(qubits[i], qureg)) {
            suffixQubits.push_back(qubits[i]);
            suffixStates.push_back(states[i]);
        }

    // overwrite given vectors
    qubits = suffixQubits;
    states = suffixStates;
}


auto getCtrlsAndTargsSwappedToMinSuffix(Qureg qureg, vector<int> ctrls, vector<int> targs) {

    // this function is called by multi-target dense matrix, and is used to find
    // targets in the prefix substate and where they can be swapped into the suffix
    // to enable subsequent embarrassingly parallel simulation. Note we seek the MIN
    // available indices in the suffix, since this minimises the stride of the local
    // simulation, improving caching performance.

    // nothing to do if all targs are already in suffix
    if (!doesGateRequireComm(qureg, targs))
        return std::tuple{ctrls, targs};

    // prepare masks to avoid quadratic nested looping
    qindex targMask = getBitMask(targs.data(), targs.size());
    qindex ctrlMask = getBitMask(ctrls.data(), ctrls.size());
    int minNonTarg = getIndOfNextRightmostZeroBit(targMask, -1);

    // prepare indices of ctrls in the given list (i.e. the inverse of ctrls), if any exist
    int maxCtrlInd = (ctrls.empty())? -1 : *std::max_element(ctrls.begin(), ctrls.end());
    vector<int> ctrlInds(maxCtrlInd+1); // bounded by ~64
    for (size_t i=0; i<ctrls.size(); i++)
        ctrlInds[ctrls[i]] = i;

    // check every target in arbitrary order, modifying our copies of targs and ctrls as we go
    for (size_t i=0; i<targs.size(); i++) {
        int targ = targs[i];

        // consider only targs in the prefix substate
        if (util_isQubitInSuffix(targ, qureg))
            continue;
            
        // if our replacement targ happens to be a ctrl... 
        if (getBit(ctrlMask, minNonTarg) == 1) {

            // find and swap that ctrl with the old targ
            int ctrlInd = ctrlInds[minNonTarg];
            ctrls[ctrlInd] = targ;

            // update our ctrl trackers
            ctrlInds[targ] = ctrlInd;
            ctrlInds[minNonTarg] = -1; // for clarity
            ctrlMask = flipTwoBits(ctrlMask, minNonTarg, targ);
        }

        // swap the prefix targ with the smallest available suffix targ
        targs[i] = minNonTarg;

        // update our targ trackers
        targMask = flipTwoBits(targMask, targ, minNonTarg);
        minNonTarg = getIndOfNextRightmostZeroBit(targMask, minNonTarg);
    }

    // the ordering in ctrls relative to the caller's ctrlStates is unchanged
    return std::tuple{ctrls, targs};
}


auto getQubitsSwappedToMaxSuffix(Qureg qureg, vector<int> qubits) {

    // this function is called by any-targ partial trace, and is used to find
    // targets in the prefix substate and where they can be swapped into the suffix
    // to enable subsequent embarrassingly parallel simulation. Note we seek the MAX
    // available indices in the suffix, since this heuristically reduces the
    // disordering of the surviving qubits after the trace, reduces the number of
    // subsequent order-restoring SWAPs

    // nothing to do if all qubits are already in suffix
    if (!doesGateRequireComm(qureg, qubits))
        return qubits;

    // prepare mask to avoid quadratic nested looping
    qindex qubitMask = getBitMask(qubits.data(), qubits.size());
    int maxFreeSuffixQubit = getIndOfNextLeftmostZeroBit(qubitMask, qureg.logNumAmpsPerNode);

    // enumerate qubits backward, modifying our copy of qubits as we go
    for (size_t i=qubits.size()-1; i-- != 0; ) {
        int qubit = qubits[i];

        // consider only qubits in the prefix substate
        if (util_isQubitInSuffix(qubit, qureg))
            continue;

        // swap the prefix qubit into the largest available suffix position
        qubits[i] = maxFreeSuffixQubit;

        // update trackers
        qubitMask = flipTwoBits(qubitMask, qubit, maxFreeSuffixQubit);
        maxFreeSuffixQubit = getIndOfNextLeftmostZeroBit(qubitMask, maxFreeSuffixQubit);
    }

    // return our modified copy
    return qubits;
}


auto getNonSwappedCtrlsAndStates(vector<int> oldCtrls, vector<int> oldStates, vector<int> newCtrls) {

    vector<int> sameCtrls(0);   sameCtrls .reserve(oldCtrls.size());
    vector<int> sameStates(0);  sameStates.reserve(oldStates.size());

    for (size_t i=0; i<oldCtrls.size(); i++)
        if (oldCtrls[i] == newCtrls[i]) {
            sameCtrls .push_back(oldCtrls[i]);
            sameStates.push_back(oldStates[i]);
        }

    return std::tuple{sameCtrls, sameStates};
}



/*
 * COMMUNICATION WRAPPERS
 */


void exchangeAmpsToBuffersWhereQubitsAreInStates(Qureg qureg, int pairRank, vector<int> qubits, vector<int> states) {

    // when there are no constraining qubits, all amps are exchanged; there is no need to pack the buffer.
    // this is typically triggered when a communicating localiser function is given no control qubits
    if (qubits.empty()) {
        comm_exchangeAmpsToBuffers(qureg, pairRank);
        return;
    }

    // otherwise, we pack and exchange only to-be-communicated amps between sub-buffers
    qindex numPacked = accel_statevec_packAmpsIntoBuffer(qureg, qubits, states);
    comm_exchangeSubBuffers(qureg, numPacked, pairRank);
}



/*
 * GETTERS
 */


qcomp localiser_statevec_getAmp(Qureg qureg, qindex globalInd) {

    // this custom routine exists, in lieu of merely invoking
    // getAmps() below, because it avoids the superfluous rank
    // enumeration and is therefore locally faster; there
    // may be applications where this is desirable (maybe)

    // only one (or all) node contains the target amp
    qcomp amp = 0;
    int sender = util_getRankContainingIndex(qureg, globalInd); // 0 if qureg not distributed

    // only node(s) containing the target amp bother to fetch 
    // it, avoiding superfluous cache costs on other nodes
    if (sender == qureg.rank) {
        qindex localInd = util_getLocalIndexOfGlobalIndex(qureg, globalInd);
        accel_statevec_getAmps_sub(&amp, qureg, localInd, 1);
    }

    // when distributed, sender shares amp to all nodes
    if (qureg.isDistributed)
        comm_broadcastAmp(sender, &amp);

    return amp;
}


void localiser_statevec_getAmps(qcomp* outAmps, Qureg qureg, qindex globalStartInd, qindex globalNumAmps) {
    
    // we do not assert state-vec, since the density matrix routine re-uses this function

    // when not distributed, all nodes merely perform direct local overwrite and finish
    if (!qureg.isDistributed) {
        accel_statevec_getAmps_sub(outAmps, qureg, globalStartInd, globalNumAmps);
        return;
    }

    // when distributed, each node will broadcast their overlap (which may be zero) with the global range
    int myRank = comm_getRank();
    int numNodes = comm_getNumNodes();

    // which they first overwrite into their local copies of out (may involve a GPU-to-CPU copy)
    if (util_areAnyVectorElemsWithinNode(myRank, qureg.numAmpsPerNode, globalStartInd, globalNumAmps)) {
        auto localInds = util_getLocalIndRangeOfVectorElemsWithinNode(myRank, qureg.numAmpsPerNode, globalStartInd, globalNumAmps);
        accel_statevec_getAmps_sub(&outAmps[localInds.localDuplicStartInd], qureg, localInds.localDistribStartInd, localInds.numElems);
    }

    // determine the overlap with each node (i.e. which amps, if any, they contribute)
    vector<qindex> globalRecvInds(numNodes);
    vector<qindex> localSendInds (numNodes);
    vector<qindex> numAmpsPerRank(numNodes, 0); // default = zero contributed amps

    for (int sendRank=0; sendRank<numNodes; sendRank++) {

        if (!util_areAnyVectorElemsWithinNode(sendRank, qureg.numAmpsPerNode, globalStartInd, globalNumAmps))
            continue;

        auto inds = util_getLocalIndRangeOfVectorElemsWithinNode(sendRank, qureg.numAmpsPerNode, globalStartInd, globalNumAmps);
        globalRecvInds[sendRank] = inds.localDuplicStartInd;
        localSendInds [sendRank] = inds.localDistribStartInd;
        numAmpsPerRank[sendRank] = inds.numElems;
    }

    // contributor nodes broadcast, all nodes receive, so that every node populates 'outAmps' fully
    comm_combineSubArrays(outAmps, globalRecvInds, localSendInds, numAmpsPerRank);
}


void localiser_densmatr_getAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols) {
    assert_localiserGivenDensMatr(qureg);

    // TODO: improve the performance!
    // this function simply serially invokes localiser_statevec_getAmps() upon 
    // every indicated column, for simplicity, and since we believe this function
    // will only ever be called upon tractably small sub-matrices. After all, the
    // user is likely to serially process outAmps themselves. Our method incurs the
    // below insignificant performance penalties:
    // - we must allocate temporary memory that is the same size as outAmps,
    //   but transposed, because we cannot make a pointer to a column of outAmps
    //   in order to invoke localiser_statevec_getAmps() thereupon. Thereafter, we
    //   serially populate outAmps with the transposed temp memory.
    // - every invocation of localiser_statevec_getAmps() invokes synchronous
    //   GPU-CPU copying (a total of #numCols), whereas a bespoke implementation
    //   could perform each non-contiguous copy asynchronously then wait
    // - every invocation invokes synchronous MPI broadcasting, whereas a bespoke
    //   method could asynch all per-col broadcasts before a final wait
    // A custom function to remedy these issues is complicated; it would involve
    // e.g. exposing MPI_Request outside of comm_routines.cpp (unacceptable for
    // compiler compatibility), or having comm_routines cache un-fulfilled asynch
    // requests, etc.

    vector<vector<qcomp>> tempOut;
    util_tryAllocMatrix(tempOut, numCols, numRows, error_localiserFailedToAllocTempMemory); // transposed dim of outAmps

    for (qindex c=0; c<numCols; c++) {
        qindex flatInd = util_getGlobalFlatIndex(qureg, startRow, startCol + c);
        localiser_statevec_getAmps(tempOut[c].data(), qureg, flatInd, numRows);
    }

    // serially overwrite outAmps = transpose(tempOut)
    for (qindex r=0; r<numRows; r++)
        for (qindex c=0; c<numCols; c++)
            outAmps[r][c] = tempOut[c][r]; // writes contiguous, reads strided
}


void localiser_fullstatediagmatr_getElems(qcomp* outElems, FullStateDiagMatr matr, qindex globalStartInd, qindex globalNumElems) {

    // printer.cpp sometimes needs to broadcast the elements of a FullStateDiagMatr;
    // since it has the same memory layout as a statevector Qureg, we spoof one!
    bool isDensMatr = false;
    bool useMultithr = false;
    bool useDistrib = matr.isDistributed;
    bool useGpuAccel = matr.isGpuAccelerated;
    Qureg qureg = qureg_populateNonHeapFields(matr.numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithr);

    // bind matr's existing CPU and GPU memory to Qureg. note that qureg's
    // comm-buffers remain null, which may be inconsistent with its distributed
    // status and would ergo fail validation; that's fine for our internal use
    qureg.cpuAmps = matr.cpuElems;
    qureg.gpuAmps = matr.gpuElems;

    localiser_statevec_getAmps(outElems, qureg, globalStartInd, globalNumElems);
}


/*
 * SETTERS
 */


void localiser_statevec_setAmps(qcomp* inAmps, Qureg qureg, qindex globalStartInd, qindex globalNumAmps) {

    // we do not assert Qureg is a state-vector, since the 
    // density matrix routine leverages this function

    // always embarrassingly parallel, since inAmps is duplicated on every node;
    // nodes simply determine which amps overlap their partition, and use them

    if (!qureg.isDistributed) {
        accel_statevec_setAmps_sub(inAmps, qureg, globalStartInd, globalNumAmps);
        return;
    }

    if (!util_areAnyVectorElemsWithinNode(qureg.rank, qureg.numAmpsPerNode, globalStartInd, globalNumAmps))
        return;

    auto range = util_getLocalIndRangeOfVectorElemsWithinNode(qureg.rank, qureg.numAmpsPerNode, globalStartInd, globalNumAmps);
    accel_statevec_setAmps_sub(&inAmps[range.localDuplicStartInd], qureg, range.localDistribStartInd, range.numElems);
}


void localiser_densmatr_setAmps(qcomp** inAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols) {
    assert_localiserGivenDensMatr(qureg);

    // TODO: improve the performance!
    // this function works by simply enumerating each column of inAmps
    // (which requires explicit preparation, because inAmps is passed
    // column-wise, rather than row-wise), passing each to the above
    // statevec routine. It is ergo similar to the naive method used by
    // localiser_densmatr_getAmps(), though is embarrassingly parallel.
    // This func allocates temporary memory as large as numRows*numCols 
    // (which could be an entire Qureg's worth), and serially computes 
    // the transpose (which can be as bad as serial iteration of the
    // whole qureg). This is grossly inefficient, and worse than merely
    // parallel-overwriting CPU memory then copying to GPU.

    vector<vector<qcomp>> tempAmps;
    util_tryAllocMatrix(tempAmps, numCols, numRows, error_localiserFailedToAllocTempMemory); // transpose of inAmps

    // serially overwrite tempAmps = transpose(inAmps)
    for (qindex c=0; c<numCols; c++)
        for (qindex r=0; r<numRows; r++)
            tempAmps[c][r] = inAmps[r][c]; // writes contiguous, reads strided

    // call the statevector function upon each column
    for (qindex c=0; c<numCols; c++) {
        qindex flatInd = util_getGlobalFlatIndex(qureg, startRow, startCol + c);
        localiser_statevec_setAmps(tempAmps[c].data(), qureg, flatInd, numRows);
    }
}


void localiser_statevec_initArbitraryPureState(Qureg qureg, qcomp* amps) {
    assert_localiserGivenStateVec(qureg);

    // always embarrassingly parallel, and merely invokes local copying
    qindex startInd = 0;
    localiser_statevec_setAmps(amps, qureg, startInd, qureg.numAmps);
}


void localiser_densmatr_initArbitraryPureState(Qureg qureg, qcomp* amps) {
    assert_localiserGivenDensMatr(qureg);

    // we cannot simply call a copy routine like the statevector case,
    // because amps will not be contiguously located in the density matrix.
    // Instead, we spoof a local CPU-only statevector, and bind amps to it
    int isDensMatr = 0;
    int useDistrib = 0;
    int useGpuAccel = 0;
    int useMultithread = 0;
    Qureg spoof = qureg_populateNonHeapFields(qureg.numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithread);

    spoof.cpuAmps = amps;
    localiser_densmatr_initPureState(qureg, spoof);
}


void localiser_densmatr_initArbitraryMixedState(Qureg qureg, qcomp** amps) {

    // TODO:
    // the current invoked implementation of setAmps() is extraordinarily
    // inefficient in the full-Qureg regime. Fix setAmps(), then update this!

    qindex startRow = 0;
    qindex startCol = 0;
    qindex numRows = powerOf2(qureg.numQubits);
    qindex numCols = numRows;
    localiser_densmatr_setAmps(amps, qureg, startRow, startCol, numRows, numCols);
}


void localiser_fullstatediagmatr_setElems(FullStateDiagMatr matr, qindex startInd, qcomp* in, qindex numElems) {

    // modification of a FullStateDiagMatr is identical to that of a 
    // statevector Qureg, so we spoof an identically-deployed Qureg
    int isDensMatr = 0;
    int useMultithread = 0; // not used by setAmps()
    Qureg spoof = qureg_populateNonHeapFields(matr.numQubits, isDensMatr, matr.isDistributed, matr.isGpuAccelerated, useMultithread);

    // we bind matr's memory to the spoofed Qureg. Note that spoof's
    // communication buffers remain nullptr which may be inconsistent
    // with its distribution status, which validation would detect
    spoof.cpuAmps = matr.cpuElems;
    spoof.gpuAmps = matr.gpuElems; // may be null

    // invoke the Qureg setter, which will modify matr's memory
    localiser_statevec_setAmps(in, spoof, startInd, numElems);

    // note that FullStateDiagMatr needs its CPU and GPU memories
    // to always be consistent with one another (like all matrix
    // types), since this is a precondition assumed by accelerator.cpp.
    // So if the above function modified GPU mem, we also trigger CPU.
    if (spoof.isGpuAccelerated) {
        spoof.isGpuAccelerated = 0;
        localiser_statevec_setAmps(in, spoof, startInd, numElems);
    }
}



/*
 * STATE INITIALISATION
 */


// defined later in this file, but repurposed here for density matrix initialisation
void mixDensityMatrixWithStatevector(qreal outProb, Qureg out, qreal inProb, Qureg in);


void localiser_statevec_initUniformState(Qureg qureg, qcomp amp) {

    // always embarrassingly parallel
    accel_statevec_initUniformState_sub(qureg, amp);
}


void localiser_statevec_initDebugState(Qureg qureg) {

    // always embarrassingly parallel
    accel_statevec_initDebugState_sub(qureg);
}


void localiser_statevec_initClassicalState(Qureg qureg, qindex globalInd) {

    // all nodes clear all amps
    accel_statevec_initUniformState_sub(qureg, 0);

    // one node (or all if qureg not distributed) modifies 1 amp
    qcomp amp = 1;
    localiser_statevec_setAmps(&amp, qureg, globalInd, 1);
}


void localiser_densmatr_initPureState(Qureg qureg, Qureg pure) {
    assert_localiserGivenDensMatr(qureg);
    assert_localiserGivenStateVec(pure);

    // we sneakily re-use the above mixing functions, since the
    // superfluous flops (multiplication of existing qureg amps
    // with zero) are completely eclipsed by the memory move costs
    qreal quregProb = 0;
    qreal pureProb = 1;
    mixDensityMatrixWithStatevector(quregProb, qureg, pureProb, pure);
}


void localiser_statevec_initUnnormalisedUniformlyRandomPureStateAmps(Qureg qureg) {
    assert_localiserGivenStateVec(qureg);

    // generation of unnormalised states is embarrassingly parallel;
    // subsequent normalisation will require reduction/communication
    accel_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(qureg);
}


Qureg createSpoofedLocalStateVecFromDensMatr(Qureg densmatr, bool &memWasAlloc) {
    assert_localiserGivenDensMatr(densmatr);

    // this function spoofs a non-distributed statevector Qureg with the 
    // same number of qubits as the given densmatr. It uses densmatr's
    // mutlithread and GPU status, and if they exist, re-uses densmatr's
    // CPU and GPU communication buffers for its main memory. If densmatr
    // has no communication buffers, temporary memory is created, which is
    // acceptable since it is expected much smaller than densmatr's local
    // memory (a factor densmatr.numColsPerNode). In that scenario,
    // memWasAlloc is overwritten to be true.

    bool isDensMatr = false;
    bool useDistrib = false;
    bool useMultithr = densmatr.isMultithreaded;  // unnecessary
    bool useGpuAccel = densmatr.isGpuAccelerated; // necessary
    Qureg pure = qureg_populateNonHeapFields(densmatr.numQubits, isDensMatr, useDistrib, useGpuAccel, useMultithr);

    // if it exists, we repurpose densmatr's existing buffer space to store
    // pure's statevector, since every node contains >=1 columns; it can fit!
    if (densmatr.isDistributed) {
        memWasAlloc = false;

        pure.cpuAmps = densmatr.cpuCommBuffer;
        pure.gpuAmps = densmatr.gpuCommBuffer; // may be null

    // otherwise, we must create temporary memory for pure. Even if we need only
    // create GPU memory, we create the matching CPU memory too for defensive design
    } else {
        memWasAlloc = true;

        pure.cpuAmps = cpu_allocArray(pure.numAmps);
        assert_localiserSuccessfullyAllocatedTempMemory(pure.cpuAmps, false); // CPU

        if (pure.isGpuAccelerated) {
            pure.gpuAmps = gpu_allocArray(pure.numAmps);
            assert_localiserSuccessfullyAllocatedTempMemory(pure.gpuAmps, true); // GPU
        }
    }

    return pure;
}


void freeSpoofedLocalStateVec(Qureg spoof, bool wasMemAlloc) {

    if (!wasMemAlloc)
        return;

    cpu_deallocArray(spoof.cpuAmps);

    if (spoof.isGpuAccelerated)
        gpu_deallocArray(spoof.gpuAmps);
}


void localiser_densmatr_initUniformlyRandomPureStateAmps(Qureg qureg) {
    assert_localiserGivenDensMatr(qureg);

    // we require a random, normalised statevector, duplicated on every
    // node, from which to initialise the density-matrix. We obtain this
    // by first spoofing a new non-distributed pure state, re-using qureg's
    // communication buffers if they exist, else creating temp memory
    bool wasMemAlloc = false;
    Qureg pure = createSpoofedLocalStateVecFromDensMatr(qureg, wasMemAlloc); // overwrites wasMemAlloc

    // initialise the spoofed pure Qureg to a normalised, uniformly random 
    // statevector; this is calling the same API function which invoked THIS 
    // very function, but instead passing a statevector
    initRandomPureState(pure); // harmlessly re-valdates

    // then, we simply initialise the density matrix in this pure state. 
    // Note that initPureState() calls mixDensityMatrixWithStatevector()
    // which writes to qureg's communication buffer only when the pure
    // qureg is distributed; we safely avoid that scenario, which would
    // conflict with our hijacking of qureg's buffer when allocated.
    localiser_densmatr_initPureState(qureg, pure);

    // if we allocated temporary memory, free it
    freeSpoofedLocalStateVec(pure, wasMemAlloc);
}


void localiser_densmatr_initMixtureOfUniformlyRandomPureStates(Qureg qureg, qindex numPureStates) {
    assert_localiserGivenDensMatr(qureg);

    // this function merely generates a random pure state, uniformly
    // mixes it into qureg (initially blank), and sequentially
    // repeats this process, increasingly mixing qureg.
    initBlankState(qureg);

    // spoof a statevector Qureg we can rnadomise, re-using
    // qureg's communication buffer memory if possible, else
    // creating temporary memory we must later free
    bool wasMemAlloc = false;
    Qureg pure = createSpoofedLocalStateVecFromDensMatr(qureg, wasMemAlloc); // overwrites wasMemAlloc

    // we wish to create qureg = sum_n^N P pure_n, where P = 1/N is
    // the uniform probability, which we achieve through repeated
    // mixing of qureg_n = a qureg_(n-1) + b_n pure_n, where
    // a = P^(1/N) = 1/N^(1/N) and b_n = a^n
    qreal a = 1. / pow(numPureStates, 1 / numPureStates);

    for (qindex n=0; n<numPureStates; n++) {
        
        qreal b = pow(a, n);
        initRandomPureState(pure); // harmlessly re-valdates
        mixDensityMatrixWithStatevector(a, qureg, b, pure);
    }

    // for large numPureStates, the above process is likely to
    // be numerically imprecise, so we finally renormalise
    setQuregToRenormalized(qureg);

    // free any allocated temp memory
    freeSpoofedLocalStateVec(pure, wasMemAlloc);
}



/*
 * SWAP
 */


void anyCtrlSwapBetweenPrefixAndPrefix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    int prefInd1 = util_getPrefixInd(targ1, qureg);
    int prefInd2 = util_getPrefixInd(targ2, qureg);

    // half of all nodes contain no to-be-swapped amps and immediately finish
    if (getBit(qureg.rank, prefInd1) == getBit(qureg.rank, prefInd2))
        return;

    // but the remaining half exchange the entirety of their amps which are in the ctrl states
    int pairRank = flipTwoBits(qureg.rank, prefInd1, prefInd2);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank, ctrls, ctrlStates);

    // and use them to overwrite their local amps satisfying ctrl states, then finish
    accel_statevec_anyCtrlSwap_subB(qureg, ctrls, ctrlStates);
}


void anyCtrlSwapBetweenPrefixAndSuffix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int suffixTarg, int prefixTarg) {

    // every node exchanges at most half its amps; those where suffixTarg bit differs from rank's fixed prefixTarg bit
    int pairRank = util_getRankWithQubitFlipped(prefixTarg, qureg);
    int suffixState =  ! util_getRankBitOfQubit(prefixTarg, qureg);

    // pack and exchange only to-be-communicated amps between sub-buffers
    vector<int> qubits = ctrls;
    vector<int> states = ctrlStates;
    qubits.push_back(suffixTarg);
    states.push_back(suffixState);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank, qubits, states);

    // we use the recevied buffer amplitudes to modify half of the local bits which satisfy ctrls
    accel_statevec_anyCtrlSwap_subC(qureg, ctrls, ctrlStates, suffixTarg, suffixState);
}


void localiser_statevec_anyCtrlSwap(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // ensure targ2 > targ1
    if (targ1 > targ2)
        std::swap(targ1, targ2);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits as relevant to communication and local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);

    // determine necessary communication
    bool comm1 = doesGateRequireComm(qureg, targ1);
    bool comm2 = doesGateRequireComm(qureg, targ2);

    if (comm2 && comm1)
        anyCtrlSwapBetweenPrefixAndPrefix(qureg, ctrls, ctrlStates, targ1, targ2);
    if (comm2 && !comm1)
        anyCtrlSwapBetweenPrefixAndSuffix(qureg, ctrls, ctrlStates, targ1, targ2);
    if (!comm2 && !comm1)
        accel_statevec_anyCtrlSwap_subA(qureg, ctrls, ctrlStates, targ1, targ2);
}



/*
 * MULTI-SWAP
 */


void anyCtrlMultiSwapBetweenPrefixAndSuffix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targsA, vector<int> targsB) {

    // this is an internal function called by the below routines which require
    // performing a sequence of SWAPs to reorder qubits, or move them into suffix.
    // the SWAPs act on unique qubit pairs and so commute.

    // TODO:
    //   - the sequence of pair-wise full-swaps should be more efficient as a
    //     "single" sequence of smaller messages sending amps directly to their
    //     final destination node. This could use a new "multiSwap" function.
    //   - if the user has compiled cuQuantum, and Qureg is GPU-accelerated, the
    //     multiSwap function should use custatevecSwapIndexBits() if local,
    //     or custatevecDistIndexBitSwapSchedulerSetIndexBitSwaps() if distributed,
    //     although the latter requires substantially more work like setting up
    //     a communicator which may be inelegant alongside our own distribution scheme.

    // perform necessary swaps to move all targets into suffix, each of which invokes communication
    for (size_t i=0; i<targsA.size(); i++)
        if (targsA[i] != targsB[i])
            anyCtrlSwapBetweenPrefixAndSuffix(qureg, ctrls, ctrlStates, targsA[i], targsB[i]);
}



/*
 * ONE-TARGET DENSE MATRIX
 */


void anyCtrlOneTargDenseMatrOnPrefix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
  
    int pairRank = util_getRankWithQubitFlipped(targ, qureg);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank, ctrls, ctrlStates);

    // extract relevant gate elements
    int bit = util_getRankBitOfQubit(targ, qureg);
    qcomp fac0 = matr.elems[bit][ bit];
    qcomp fac1 = matr.elems[bit][!bit];

    // update local amps using received amps in buffer
    accel_statevec_anyCtrlOneTargDenseMatr_subB(qureg, ctrls, ctrlStates, fac0, fac1);
}


void localiser_statevec_anyCtrlOneTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr, bool conj) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits as relevant to communication and local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);

    if (conj) 
        matr = util_getConj(matr);

    // perform embarrassingly parallel routine or communication-inducing swaps
    doesGateRequireComm(qureg, targ)?
        anyCtrlOneTargDenseMatrOnPrefix(qureg, ctrls, ctrlStates, targ, matr) :
        accel_statevec_anyCtrlOneTargDenseMatr_subA(qureg, ctrls, ctrlStates, targ, matr);
}



/*
 * TWO-TARGET & ANY-TARGET DENSE MATRIX
 *
 * which are intermixed, despite each having their own local backend 
 * implementations, because they use identical communication logic
 */


void anyCtrlTwoOrAnyTargDenseMatrOnSuffix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr2 matr, bool conj) {
    if (conj) matr = util_getConj(matr);
    accel_statevec_anyCtrlTwoTargDenseMatr_sub(qureg, ctrls, ctrlStates, targs[0], targs[1], matr);
}
void anyCtrlTwoOrAnyTargDenseMatrOnSuffix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr  matr, bool conj) {
    accel_statevec_anyCtrlAnyTargDenseMatr_sub(qureg, ctrls, ctrlStates, targs, matr, conj);
}


// T can be CompMatr2 or CompMatr
template <typename T>
void anyCtrlTwoOrAnyTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, T matr, bool conj) {

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    // skip straight to embarrasingly parallel simulation if possible
    if (!doesGateRequireComm(qureg, targs)) {

        // using only the suffix ctrls
        removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);
        anyCtrlTwoOrAnyTargDenseMatrOnSuffix(qureg, ctrls, ctrlStates, targs, matr, conj);
        return;
    }

    // find suffix positions for all prefix targs, moving colliding ctrls out of the way
    auto [newCtrls, newTargs] = getCtrlsAndTargsSwappedToMinSuffix(qureg, ctrls, targs);

    // only unmoved ctrls can be applied to the swaps, to accelerate them
    auto [unmovedCtrls, unmovedCtrlStates] = getNonSwappedCtrlsAndStates(newCtrls, ctrlStates, newTargs); 

    // perform necessary swaps to move all targets into suffix, invoking communication (swaps are real, so no need to conj)
    anyCtrlMultiSwapBetweenPrefixAndSuffix(qureg, unmovedCtrls, unmovedCtrlStates, targs, newTargs);

    // if the moved ctrls do not eliminate this node's need for local simulation...
    if (doAnyLocalStatesHaveQubitValues(qureg, newCtrls, ctrlStates)) {

        // perform embarrassingly parallel simulation using only the new suffix ctrls
        removePrefixQubitsAndStates(qureg, newCtrls, ctrlStates);
        anyCtrlTwoOrAnyTargDenseMatrOnSuffix(qureg, newCtrls, ctrlStates, newTargs, matr, conj);
    }

    // undo swaps, again invoking communication
    anyCtrlMultiSwapBetweenPrefixAndSuffix(qureg, unmovedCtrls, unmovedCtrlStates, targs, newTargs);
}


void localiser_statevec_anyCtrlTwoTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr, bool conj) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    anyCtrlTwoOrAnyTargDenseMatr(qureg, ctrls, ctrlStates, {targ1,targ2}, matr, conj);
}


void localiser_statevec_anyCtrlAnyTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr, bool conj) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // despite our use of compile-time templating, the bespoke one-targ routines are still faster 
    // than this any-targ routine when given a single target, because they can leverage a bespoke
    // communication pattern (rather than swapping qubits into suffix), and pass the matrix elems
    // to GPU kernels via arguments rather than global memory, which is faster for threads to read.
    // Callers may however still choose this function (rather than the one-qubit specific one) for 
    // its convenient generality, so we divert to the one-targ routine when possible, copying the 
    // heap CPU matrix (assumed consistent with GPU memory) into stack memory
    if (targs.size() == 1)
        localiser_statevec_anyCtrlOneTargDenseMatr(qureg, ctrls, ctrlStates, targs[0], getCompMatr1(matr.cpuElems), conj);
    
    // similarly, bespoke two-targ routines are preferable although they offer no communication
    // benefit because they call the same any-targ localiser, but still accelerate GPU memory access
    else if (targs.size() == 2)
        localiser_statevec_anyCtrlTwoTargDenseMatr(qureg, ctrls, ctrlStates, targs[0], targs[1], getCompMatr2(matr.cpuElems), conj);
    
    // call the any-targ routine when given 3 or more targs, which may still invoke bespoke,
    // fixed-targ instances of backend templated functions depending the number of targs
    else
        anyCtrlTwoOrAnyTargDenseMatr(qureg, ctrls, ctrlStates, targs, matr, conj);
}



/*
 * ANY-TARGET DIAGONAL MATRIX
 *
 * which have num-target specific implementations (e.g. for avoiding
 * GOU memory, if possible), but identical communication logic
 * because diagonals are always embarrassingly parallel
 */


void localiser_statevec_anyCtrlOneTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr, bool conj) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    if (conj)
        matr = util_getConj(matr);

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);
    accel_statevec_anyCtrlOneTargDiagMatr_sub(qureg, ctrls, ctrlStates, targ, matr);
}


void localiser_statevec_anyCtrlTwoTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr, bool conj) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    if (conj)
        matr = util_getConj(matr);

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);
    accel_statevec_anyCtrlTwoTargDiagMatr_sub(qureg, ctrls, ctrlStates, targ1, targ2, matr);
}


void localiser_statevec_anyCtrlAnyTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent, bool conj) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);
    accel_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, ctrls, ctrlStates, targs, matr, exponent, conj);
}



/*
 * ALL-TARGET DIAGONAL MATRIX
 */


void localAllTargDiagMatrOnDistribStatevector(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    // this is a strange scenario (why distribute qureg but not
    // an equally-sized matrix?) we support merely for defensive-
    // design. In this scenario, every node has the full matrix
    // amps and needs only to consult a contiguous subset. To
    // avoid bespoke but trivially different backend functions,
    // we spoof a distributed copy of 'matr' with pointers
    // offset to the elements needed by this node, and call the
    // equally-distributed backend function.

    qindex offset = util_getGlobalIndexOfFirstLocalAmp(qureg);

    // cheaply copy matr, spoofed as distributed (to pass internal checks)
    FullStateDiagMatr copy = matr;
    copy.isDistributed = true;
    copy.numElemsPerNode = qureg.numAmpsPerNode; // not consulted, but for safety

    // offset pointers to qureg's existing memory, avoiding de-referencing nullptr (illegal)
    copy.cpuElems = &copy.cpuElems[offset];
    copy.gpuElems = (copy.gpuElems == nullptr)? copy.gpuElems : &copy.gpuElems[offset];

    // invoke backend; copy is in stack and will be auto-freed
    accel_statevec_allTargDiagMatr_sub(qureg, copy, exponent);
}


void localiser_statevec_allTargDiagMatr(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {
    assert_localiserGivenStateVec(qureg);

    // matr and qureg are equal dimension, but may be distributed differently
    bool quregDist = qureg.isDistributed;
    bool matrDist = matr.isDistributed;

    // cannot distribute only matr; qureg has no buffer space to receive a broadcast.
    // allocating temporary buffer space is too dangerous, since it would be the same
    // size as matr and qureg, potentially increasing user memory by a factor x1.5.
    if (!quregDist && matrDist)
        error_localiserGivenDistribMatrixAndLocalQureg();

    // embarrassingly parallel when both distributed or both local
    if (quregDist == matrDist)
        accel_statevec_allTargDiagMatr_sub(qureg, matr, exponent);
    
    // embarrasingly parallel when only qureg is distributed (all nodes have all needed matr elems)
    if (quregDist && !matrDist)
        localAllTargDiagMatrOnDistribStatevector(qureg, matr, exponent);
}


void localiser_densmatr_allTargDiagMatr(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool multiplyOnly) {
    assert_localiserGivenDensMatr(qureg);

    // the diagonal matr has quadratically fewer elements than the density-matrix
    // qureg, so could in-theory be effected (without catastrophic performance) as
    // an N/2-qubit DiagMatr upon an N-qubit statevector. This requires O(N) bitwise
    // operations per-iteration which might cause a slowdown in non-memory-bandwidth
    // bound settings (e.g. 8 qubit Quregs). So we here use an O(1) bespoke method.

    // since Qureg is quadratically bigger than matr, it is likely they have different
    // distributions (matr is probably local). Because every column of qureg is 
    // dot-multiplied with the full matr, every node requires all matr elements
    bool quregDist = qureg.isDistributed;
    bool matrDist = matr.isDistributed;

    // cannot distribute only matr; qureg has no buffer space to receive a broadcast.
    // in theory, we could allocate temporary buffer space which would only be 
    // quadratically smaller than qureg; but this is a ludicrous scenario to support.
    if (!quregDist && matrDist) {
        error_localiserGivenDistribMatrixAndLocalQureg();
        return;
    }

    // when the matrix is not distributed, we call the same routine despite whether qureg 
    // is distributed or not; that merely changes how many qureg columns get updated
    if (!matrDist) {
        accel_densmatr_allTargDiagMatr_subA(qureg, matr, exponent, multiplyOnly);
        return;
    }

    // finally, when both are distributed, qureg has buffer space to receive all matr
    comm_combineElemsIntoBuffer(qureg, matr);

    // matr elems are inside qureg buffer, but we still pass matr struct along to
    // accelerator, because it is going to perform mischief to re-use subA().
    accel_densmatr_allTargDiagMatr_subB(qureg, matr, exponent, multiplyOnly); 
}



/*
 * ANY-TARGET ANY-TYPE MATRIX 
 * 
 * This is merely a convenient gateway for callers to automatically
 * dispatch to the above specific functions, based on matrix type
 */


template <class T>
void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, T matr, bool conj) {
    if constexpr (util_isDiagMatr <T>()) localiser_statevec_anyCtrlAnyTargDiagMatr(qureg,  ctrls, ctrlStates, targs, matr, 1, conj); // exponent=1
    if constexpr (util_isDiagMatr1<T>()) localiser_statevec_anyCtrlOneTargDiagMatr(qureg,  ctrls, ctrlStates, targs[0], matr, conj);
    if constexpr (util_isDiagMatr2<T>()) localiser_statevec_anyCtrlTwoTargDiagMatr(qureg,  ctrls, ctrlStates, targs[0], targs[1], matr, conj);
    if constexpr (util_isCompMatr <T>()) localiser_statevec_anyCtrlAnyTargDenseMatr(qureg, ctrls, ctrlStates, targs, matr, conj);
    if constexpr (util_isCompMatr1<T>()) localiser_statevec_anyCtrlOneTargDenseMatr(qureg, ctrls, ctrlStates, targs[0], matr, conj);
    if constexpr (util_isCompMatr2<T>()) localiser_statevec_anyCtrlTwoTargDenseMatr(qureg, ctrls, ctrlStates, targs[0], targs[1], matr, conj);
}

template void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg, vector<int>, vector<int>, vector<int>, DiagMatr,  bool);
template void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg, vector<int>, vector<int>, vector<int>, DiagMatr1, bool);
template void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg, vector<int>, vector<int>, vector<int>, DiagMatr2, bool);
template void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg, vector<int>, vector<int>, vector<int>, CompMatr,  bool);
template void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg, vector<int>, vector<int>, vector<int>, CompMatr1, bool);
template void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg, vector<int>, vector<int>, vector<int>, CompMatr2, bool);



/*
 * PAULI TENSORS AND GADGETS
 */


extern bool paulis_containsXOrY(PauliStr str);
extern vector<int> paulis_getInds(PauliStr str);
extern std::array<vector<int>,3> paulis_getSeparateInds(PauliStr str, Qureg qureg);
extern int paulis_getPrefixZSign(Qureg qureg, vector<int> prefixZ) ;
extern qcomp paulis_getPrefixPaulisElem(Qureg qureg, vector<int> prefixY, vector<int> prefixZ);


void anyCtrlZTensorOrGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, bool isGadget, qreal phase) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);

    // prefixZ merely applies a node-wide factor to fac0 and fac1
    auto [prefixZ, suffixZ] = util_getPrefixAndSuffixQubits(targs, qureg);
    int sign = paulis_getPrefixZSign(qureg, prefixZ);
    
    // tensor multiplies +-1, gadget multiplies exp(+- i phase)
    qcomp fac0 = (isGadget)? exp(+ phase * sign * 1_i) : +1 * sign;
    qcomp fac1 = (isGadget)? exp(- phase * sign * 1_i) : -1 * sign;

    // simulation is always embarrassingly parallel
    accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(qureg, ctrls, ctrlStates, suffixZ, fac0, fac1);
}


void anyCtrlPauliTensorOrGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qcomp ampFac, qcomp pairAmpFac) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // this routine is invalid for str=ZI
    if (!paulis_containsXOrY(str))
        error_localiserGivenPauliStrWithoutXorY();

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalStatesHaveQubitValues(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);

    // partition non-Id Paulis into prefix and suffix, since...
    // - prefix X,Y determine communication, because they apply bit-not to rank
    // - prefix Y,Z determine node-wide coefficient, because they contain rank-determined !=1 elements
    // - suffix X,Y,Z determine local amp coefficients
    auto [targsX, targsY, targsZ] = paulis_getSeparateInds(str, qureg);
    auto [prefixX, suffixX] = util_getPrefixAndSuffixQubits(targsX, qureg);
    auto [prefixY, suffixY] = util_getPrefixAndSuffixQubits(targsY, qureg);
    auto [prefixZ, suffixZ] = util_getPrefixAndSuffixQubits(targsZ, qureg);

    // scale pair amp's coefficient by node-wide coeff 
    pairAmpFac *= paulis_getPrefixPaulisElem(qureg, prefixY, prefixZ); // 1 when embarrassingly parallel

    // embarrassingly parallel when there is only Z's in prefix
    if (prefixX.empty() && prefixY.empty()) {
        accel_statevector_anyCtrlPauliTensorOrGadget_subA(qureg, ctrls, ctrlStates, suffixX, suffixY, suffixZ, ampFac, pairAmpFac);
        return;
    }

    // otherwise, we pair-wise communicate amps satisfying ctrls
    auto prefixXY = util_getConcatenated(prefixX, prefixY);
    int pairRank = util_getRankWithQubitsFlipped(prefixXY, qureg);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank, ctrls, ctrlStates);

    // ctrls reduce communicated amps, so received buffer is compacted;
    // we must ergo prepare a no-ctrl XY mask for accessing buffer elems
    auto sortedCtrls = util_getSorted(ctrls);
    auto suffixMaskXY = util_getBitMask(util_getConcatenated(suffixX, suffixY));
    auto bufferMaskXY = removeBits(suffixMaskXY, sortedCtrls.data(), sortedCtrls.size());

    accel_statevector_anyCtrlPauliTensorOrGadget_subB(qureg, ctrls, ctrlStates, suffixX, suffixY, suffixZ, ampFac, pairAmpFac, bufferMaskXY);
}


void localiser_statevec_anyCtrlPauliTensor(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qcomp factor) {

    // this function accepts a global factor, so that density matrices can effect conj(pauli)

    if (paulis_containsXOrY(str)) {

        // (X|Y)|0/1> ~ |1/0>
        qcomp ampFac     = 0 * factor;
        qcomp pairAmpFac = 1 * factor;
        anyCtrlPauliTensorOrGadget(qureg, ctrls, ctrlStates, str, ampFac, pairAmpFac);

    } else {
        // global factor is inapplicable to all-Z and is ignored
        if (factor != qcomp(1,0))
            error_localiserGivenNonUnityGlobalFactorToZTensor();

        bool isGadget = false;
        qreal phase = 0; // ignored
        anyCtrlZTensorOrGadget(qureg, ctrls, ctrlStates, paulis_getInds(str), isGadget, phase);
    }
}


void localiser_statevec_anyCtrlPhaseGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qreal phase) {

    bool isGadget = true;
    anyCtrlZTensorOrGadget(qureg, ctrls, ctrlStates, targs, isGadget, phase); 
}


void localiser_statevec_anyCtrlPauliGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qreal phase) {

    // when str=IZ, we must use the above bespoke algorithm
    if (!paulis_containsXOrY(str)) {
        localiser_statevec_anyCtrlPhaseGadget(qureg, ctrls, ctrlStates, paulis_getInds(str), phase);
        return;
    }

    qcomp ampFac     = cos(phase);
    qcomp pairAmpFac = sin(phase) * 1_i;
    anyCtrlPauliTensorOrGadget(qureg, ctrls, ctrlStates, str, ampFac, pairAmpFac);
}



/*
 * QUREG COMBINATION
 */


void localiser_statevec_setQuregToSuperposition(qcomp facOut, Qureg outQureg, qcomp fac1, Qureg inQureg1, qcomp fac2, Qureg inQureg2) {

    // TODO:
    // this function requires (as validated) distributions are identical.
    // It would be trivial to generalise this so that Qureg distributions
    // can differ (we merely spoof local Quregs, offsetting their memory).
    // They must still however be identically GPU-accelerated; this is a
    // low priority because this situation is non-sensical

    // given Qureg dimensions must match, this is always embarrassingly parallel
    accel_statevec_setQuregToSuperposition_sub(facOut, outQureg, fac1, inQureg1, fac2, inQureg2);
}


void mixDensityMatrixWithStatevector(qreal outProb, Qureg out, qreal inProb, Qureg in) {

    // we can handle 3 out of 4 possible combinations of distribution,
    // and accelerator.hpp will handle every combination of GPU-accel
    bool outDist = out.isDistributed;
    bool inDist = in.isDistributed;

    // illegal to distribute only the smaller Qureg; 'out' has no buffer space to receive it.
    // in theory, we could allocate a temporary receive buffer since it is merely quadratically
    // smaller than 'out' so will be a negligible memory overhead; but this scenario (having a
    // smaller distributed Qureg) is completely ludicrous and unworth supporting
    if (!outDist && inDist)
        error_mixQuregsAreLocalDensMatrAndDistribStatevec();

    // both non-distributed is trivial
    if (!outDist && !inDist)
        accel_densmatr_mixQureg_subB(outProb, out, inProb, in);

    // both distributed requires broadcasting 'in' into every node's 'out' buffer
    if (outDist && inDist) {
        comm_combineAmpsIntoBuffer(out, in); // uses same buffer that subC() consults
        accel_densmatr_mixQureg_subC(outProb, out, inProb);
    }

    // only 'out' being distributed means simulation is embarrasingly parallel, 
    // because the full 'in' is already known on every node
    if (outDist && !inDist)
        accel_densmatr_mixQureg_subD(outProb, out, inProb, in);
}


void localiser_densmatr_mixQureg(qreal outProb, Qureg out, qreal inProb, Qureg in) {
    assert_localiserGivenDensMatr(out);

    (in.isDensityMatrix)?
        accel_densmatr_mixQureg_subA(outProb, out, inProb, in): // trivial
        mixDensityMatrixWithStatevector(outProb, out, inProb, in);
}



/*
 * DEPHASING
 */


void localiser_densmatr_oneQubitDephasing(Qureg qureg, int qubit, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // both methods are embarrassingly parallel
    (util_isBraQubitInSuffix(qubit, qureg))?
        accel_densmatr_oneQubitDephasing_subA(qureg, qubit, prob):
        accel_densmatr_oneQubitDephasing_subB(qureg, qubit, prob);
}


void localiser_densmatr_twoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // relative size of qubit1 and qubit2 does not matter
    int leftQubit = std::max(qubit1, qubit2);

    // both methods are embarrassingly parallel
    (util_isBraQubitInSuffix(leftQubit, qureg))?
        accel_densmatr_twoQubitDephasing_subA(qureg, qubit1, qubit2, prob):
        accel_densmatr_twoQubitDephasing_subB(qureg, qubit1, qubit2, prob);
}



/*
 * ONE-QUBIT DEPOLARISING
 */


void oneQubitDepolarisingOnPrefix(Qureg qureg, int ketQubit, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // pack and exchange amps to buffers where local ket qubit and fixed-prefix-bra qubit agree
    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    int pairRank = util_getRankWithBraQubitFlipped(ketQubit, qureg);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank, {ketQubit}, {braBit});

    // use received sub-buffer to update local amps
    accel_densmatr_oneQubitDepolarising_subB(qureg, ketQubit, prob);
}


void localiser_densmatr_oneQubitDepolarising(Qureg qureg, int qubit, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // perform embarrassingly parallel routine or pairwise communication
    (doesChannelRequireComm(qureg, qubit))?
        oneQubitDepolarisingOnPrefix(qureg, qubit, prob):
        accel_densmatr_oneQubitDepolarising_subA(qureg, qubit, prob);
}



/*
 * TWO-QUBIT DEPOLARISING
 */


void twoQubitDepolarisingOnPrefixAndSuffix(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // scale 25% of amps; precisely those which are not communicated
    accel_densmatr_twoQubitDepolarising_subC(qureg, ketQb1, ketQb2, prob);

    // pack an eighth of the buffer with pair-summed amps
    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);
    qindex numPacked = accel_statevec_packPairSummedAmpsIntoBuffer(qureg, ketQb1, ketQb2, braQb1, braBit2);

    // exchange sub-buffers
    int pairRank = util_getRankWithBraQubitFlipped(ketQb2, qureg);
    comm_exchangeSubBuffers(qureg, numPacked, pairRank);

    // update 25% of local amps using received buffer amps
    accel_densmatr_twoQubitDepolarising_subD(qureg, ketQb1, ketQb2, prob);
}


void twoQubitDepolarisingOnPrefixAndPrefix(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    // scale 25% of (non-communicated) amps
    accel_densmatr_twoQubitDepolarising_subE(qureg, ketQb1, ketQb2, prob);

    // pack and swap 25% of buffer, and use it to modify 25% of local amps
    int pairRank1 = util_getRankWithBraQubitFlipped(ketQb1, qureg);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank1, {ketQb1,ketQb2}, {braBit1,braBit2});
    accel_densmatr_twoQubitDepolarising_subF(qureg, ketQb1, ketQb2, prob);

    // pack and swap another 25% of buffer (we could pack during subE, but we choose not to)
    int pairRank2 = util_getRankWithBraQubitFlipped(ketQb2, qureg);
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank2, {ketQb1,ketQb2}, {braBit1,braBit2});
    accel_densmatr_twoQubitDepolarising_subG(qureg, ketQb1, ketQb2, prob);
}


void localiser_densmatr_twoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // ensure qubit2 > qubit1
    if (qubit1 > qubit2)
        std::swap(qubit1, qubit2);

    // determine necessary communication
    bool comm1 = doesChannelRequireComm(qureg, qubit1);
    bool comm2 = doesChannelRequireComm(qureg, qubit2);

    if (comm2 && comm1)
        twoQubitDepolarisingOnPrefixAndPrefix(qureg, qubit1, qubit2, prob);
    if (comm2 && !comm1)
        twoQubitDepolarisingOnPrefixAndSuffix(qureg, qubit1, qubit2, prob);
    if (!comm2 && !comm1) {
        accel_densmatr_twoQubitDepolarising_subA(qureg, qubit1, qubit2, prob);
        accel_densmatr_twoQubitDepolarising_subB(qureg, qubit1, qubit2, prob);
    }
}



/*
 * PAULI CHANNEL
 */


void oneQubitPauliChannelOnPrefix(Qureg qureg, int ketQubit, qreal probI, qreal probX, qreal probY, qreal probZ) {
    assert_localiserGivenDensMatr(qureg);

    // exchange all amps with pair node
    int pairRank = util_getRankWithBraQubitFlipped(ketQubit, qureg);
    comm_exchangeAmpsToBuffers(qureg, pairRank);

    // use received buffer to update local amps
    accel_densmatr_oneQubitPauliChannel_subB(qureg, ketQubit, probI, probX, probY, probZ);
}


void localiser_densmatr_oneQubitPauliChannel(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {
    assert_localiserGivenDensMatr(qureg);

    // infer the no-error probability
    qreal probI = 1 - probX - probY - probZ;

    (doesChannelRequireComm(qureg, qubit))?
        oneQubitPauliChannelOnPrefix(qureg, qubit, probI, probX, probY, probZ):
        accel_densmatr_oneQubitPauliChannel_subA(qureg, qubit, probI, probX, probY, probZ);
}


// twoQubitPauliChannel() is regrettably too difficult; the communication model cannot be 
// simplified the way it was in twoQubitDepolarising() which levereaged the uniform
// coefficients. It is not clear whether arbitrary coefficients, which cause many more
// amplitudes to mix, can ever be performed in a sequence of pairwise communication



/*
 * AMPLITUDE DAMPING
 */


void oneQubitDampingOnPrefix(Qureg qureg, int ketQubit, qreal prob) {

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    int pairRank = util_getRankWithBraQubitFlipped(ketQubit, qureg);
    qindex numAmps = qureg.numAmpsPerNode / 2;

    // half of all nodes...
    if (braBit == 1) {

        // pack and async send half the buffer
        accel_statevec_packAmpsIntoBuffer(qureg, {ketQubit}, {1});
        comm_asynchSendSubBuffer(qureg, numAmps, pairRank);

        // scale the local amps which were just sent
        accel_densmatr_oneQubitDamping_subB(qureg, ketQubit, prob);
    }

    // all nodes scale the other half of their local amps
    accel_densmatr_oneQubitDamping_subC(qureg, ketQubit, prob);

    // the other remaining half of all nodes...
    if (braBit == 0) {

        // receive the async-sent buffer
        comm_receiveArrayToBuffer(qureg, numAmps, pairRank);
        accel_densmatr_oneQubitDamping_subD(qureg, ketQubit, prob);
    }

    // prevent asynch senders from proceeding so their buffer isn't prematurely modified
    comm_sync();
}


void localiser_densmatr_oneQubitDamping(Qureg qureg, int qubit, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    (doesChannelRequireComm(qureg, qubit))?
        oneQubitDampingOnPrefix(qureg, qubit, prob):
        accel_densmatr_oneQubitDamping_subA(qureg, qubit, prob);
}



/*
 * SUPEROPERATORS AND KRAUS
 */


CompMatr getCompMatrFromSuperOp(SuperOp op) {

    return (CompMatr) {
        // superoperator acts on twice as many qubits
        .numQubits = 2 * op.numQubits,
        .numRows = op.numRows,

        // isUnitary will not be consulted (we have passed validation)
        .isUnitary = nullptr,

        // copy pointers (noting cpuElems is 2D/nested)
        .cpuElems = op.cpuElems,
        .gpuElemsFlat = op.gpuElemsFlat
    };
}


void localiser_densmatr_superoperator(Qureg qureg, SuperOp op, vector<int> ketTargs) {
    assert_localiserGivenDensMatr(qureg);

    // effect the superoperator as a (non-conjugated) dense matrix on the ket + bra qubits
    bool conj = false;
    auto braTargs = util_getBraQubits(ketTargs, qureg);
    auto allTargs = util_getConcatenated(ketTargs, braTargs);
    CompMatr matr = getCompMatrFromSuperOp(op);
    localiser_statevec_anyCtrlAnyTargDenseMatr(qureg, {}, {}, allTargs, matr, conj);
}


void localiser_densmatr_krausMap(Qureg qureg, KrausMap map, vector<int> ketTargs) {
    
    // Kraus map is simulated through its existing superoperator
    localiser_densmatr_superoperator(qureg, map.superop, ketTargs);
}



/*
 * PARTIAL TRACE
 */


auto getNonTracedQubitOrder(Qureg qureg, vector<int> originalTargs, vector<int> revisedTargs) {

    // prepare a list of all the qureg's qubits when treated as a statevector
    vector<int> allQubits(2*qureg.numQubits);
    for (size_t q=0; q<allQubits.size(); q++)
        allQubits[q] = q;
    
    // determine the ordering of all the Qureg's qubits after swaps
    for (size_t i=0; i<originalTargs.size(); i++) {
        int qb1 = originalTargs[i];
        int qb2 = revisedTargs[i];
        if (qb1 != qb2)
            std::swap(allQubits[qb1], allQubits[qb2]);
    }

    // use a mask to avoid quadratic nested iteration below
    qindex revisedMask = util_getBitMask(revisedTargs);

    // retain only non-targeted qubits
    vector<int> remainingQubits(allQubits.size() - originalTargs.size());
    for (size_t q=0; q<allQubits.size(); q++)
        if (!getBit(revisedMask, q))
            remainingQubits.push_back(allQubits[q]);

    // shift down remaining qubits to be contiguous...
    qindex remainingMask = util_getBitMask(remainingQubits);
    for (int &qubit : remainingQubits) {
        int bound = qubit;

        // by subtracting the number of smaller un-targeted qubits from each qubit index
        for (int i=0; i<bound; i++)
            qubit -= ! getBit(remainingMask, i);
    }

    // return the ordering, i.e. a list [0, #final-qubits)
    return remainingQubits;
}


void reorderReducedQureg(Qureg inQureg, Qureg outQureg, vector<int> allTargs, vector<int> suffixTargs) {

    // TODO: 
    // this function performs a sequence of SWAPs which are NOT necessarily upon disjoint qubits,
    // and ergo do not commute. We still however may be able to effect this more efficiently in
    // a single communicating operation rather than this sequence of SWAP gates, and might still
    // even be able to use cuQuantum's distributed bit index swaps API. Check this!

    // determine the relative ordering of outQureg's remaining qubits
    auto remainingQubits = getNonTracedQubitOrder(inQureg, allTargs, suffixTargs);

   // perform additional swaps to re-order the remaining qubits (heuristically starting from back)
    for (int qubit=(int)remainingQubits.size(); qubit-- != 0; ) {

        // locate the next qubit which is out of its sorted position
        if (remainingQubits[qubit] == qubit)
            continue;

        // qubit is misplaced; locate its position among the remaining qubits
        int pair = 0;
        while (remainingQubits[pair] != qubit)
            pair++;
        
        // and swap it directly to its required position, triggering any communication scenario (I think)
        localiser_statevec_anyCtrlSwap(outQureg, {}, {}, qubit, pair);
        std::swap(remainingQubits[qubit], remainingQubits[pair]);
    }
}


void partialTraceOnSuffix(Qureg inQureg, Qureg outQureg, vector<int> ketTargs) {

    auto braTargs = util_getBraQubits(ketTargs, inQureg);
    accel_densmatr_partialTrace_sub(inQureg, outQureg, ketTargs, braTargs);
}


void partialTraceOnPrefix(Qureg inQureg, Qureg outQureg, vector<int> ketTargs) {

    // all ketTargs (pre-sorted) are in the suffix, but one or more braTargs are in the prefix
    auto braTargs = util_getBraQubits(ketTargs, inQureg); // sorted
    auto allTargs = util_getSorted(ketTargs, braTargs);   // sorted
    auto sufTargs = getQubitsSwappedToMaxSuffix(inQureg, allTargs); // arbitrarily ordered

    // swap iniQureg's prefix bra-qubits into suffix, invoking communication
    anyCtrlMultiSwapBetweenPrefixAndSuffix(inQureg, {}, {}, sufTargs, allTargs);

    // use the second half of sufTargs as the pair targs, which are now all in the suffix,
    // to perform embarrassingly parallel overwriting of outQureg
    vector<int> pairTargs(sufTargs.begin() + ketTargs.size(), sufTargs.end()); // arbitrarily ordered
    accel_densmatr_partialTrace_sub(inQureg, outQureg, ketTargs, pairTargs);

    // restore the relative order of outQureg's remaining qubits using SWAPs
    reorderReducedQureg(inQureg, outQureg, allTargs, sufTargs);

    // undo the swaps on inQureg
    anyCtrlMultiSwapBetweenPrefixAndSuffix(inQureg, {}, {}, sufTargs, allTargs);
}


void localiser_densmatr_partialTrace(Qureg inQureg, Qureg outQureg, vector<int> targs) {
    assert_localiserPartialTraceGivenCompatibleQuregs(inQureg, outQureg, targs.size());

    // TODO: 
    // this function requires inQureg and outQureg are both or neither distributed;
    // it does not support the (potentially reasonable) situation when only outQureg
    // is not-distributed, perhaps because it is too small. It is not simple to
    // extend our method to this case. We could force distribution of inQureg until
    // the end of the routine (where we might broadcast it back to non-distributed),
    // and even temporarily relax the condition that each node contains >=1 column,
    // but we would still be restricted by the requirement each node contains >=1 amp.
    // Think about whether we can relax this!

    // sorted targets needed by subsequent bitwise insertions
    targs = util_getSorted(targs);

    (doesChannelRequireComm(inQureg, targs.back()))?
        partialTraceOnPrefix(inQureg, outQureg, targs):
        partialTraceOnSuffix(inQureg, outQureg, targs);
}



/*
 * PROBABILITIES
 */


qreal localiser_statevec_calcTotalProb(Qureg qureg) {
    
    // not restricted to statevecs; density matrices use
    // a different routine for calcTotalProb, but they use
    // this routine for calcHilbertSchmidtDistance

    qreal prob = accel_statevec_calcTotalProb_sub(qureg);
    if (qureg.isDistributed)
        comm_reduceReal(&prob);

    return prob;
}


qreal localiser_densmatr_calcTotalProb(Qureg qureg) {
    assert_localiserGivenDensMatr(qureg);

    qreal prob = accel_densmatr_calcTotalProb_sub(qureg);
    if (qureg.isDistributed)
        comm_reduceReal(&prob);

    return prob;
}


qreal localiser_statevec_calcProbOfMultiQubitOutcome(Qureg qureg, vector<int> qubits, vector<int> outcomes) {
    assert_localiserGivenStateVec(qureg);

    qreal prob = 0;

    // only nodes containing amps with the specified outcomes need to compute probs
    if (doAnyLocalStatesHaveQubitValues(qureg, qubits, outcomes)) {

        // and do so using only the suffix qubits/outcomes
        removePrefixQubitsAndStates(qureg, qubits, outcomes);
        prob += accel_statevec_calcProbOfMultiQubitOutcome_sub(qureg, qubits, outcomes);
    }

    // but all nodes must sum their probabilities (unless qureg was cloned per-node), for conensus
    if (qureg.isDistributed)
        comm_reduceReal(&prob);

    return prob;
}


qreal localiser_densmatr_calcProbOfMultiQubitOutcome(Qureg qureg, vector<int> qubits, vector<int> outcomes) {
    assert_localiserGivenDensMatr(qureg);

    qreal prob = 0;

    // only some nodes contain columns wherein the diagonal element corresponds to
    // the given qubit configuration; this is determined by the bra-qubits (because
    // the ket-qubits are always in the suffix partition)
    auto braQubits = util_getBraQubits(qubits, qureg);

    if (doAnyLocalStatesHaveQubitValues(qureg, braQubits, outcomes)) {

        // such nodes need to know all ket qubits (which are all suffix)
        prob += accel_densmatr_calcProbOfMultiQubitOutcome_sub(qureg, qubits, outcomes);
    }

    // all nodes must sum their probabilities (unless qureg was cloned per-node), for consensus
    if (qureg.isDistributed)
        comm_reduceReal(&prob);

    return prob;
}


void localiser_statevec_calcProbsOfAllMultiQubitOutcomes(qreal* outProbs, Qureg qureg, vector<int> qubits) {
    assert_localiserGivenStateVec(qureg);

    // each node independently populates local outProbs
    accel_statevec_calcProbsOfAllMultiQubitOutcomes_sub(outProbs, qureg, qubits);

    // nodes sum their arrays
    if (qureg.isDistributed)
        comm_reduceReals(outProbs, powerOf2(qubits.size()));
}


void localiser_densmatr_calcProbsOfAllMultiQubitOutcomes(qreal* outProbs, Qureg qureg, vector<int> qubits) {
    assert_localiserGivenDensMatr(qureg);

    // each node independently populates local outProbs
    accel_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(outProbs, qureg, qubits);

    // nodes sum their arrays
    if (qureg.isDistributed)
        comm_reduceReals(outProbs, powerOf2(qubits.size()));
}



/*
 * EXPECTATION VALUES
 */


qcomp localiser_statevec_calcExpecPauliStr(Qureg qureg, PauliStr str) {
    assert_localiserGivenStateVec(qureg);

    qcomp value = 0;

    // partition non-Id Paulis into prefix and suffix, since...
    // - prefix X,Y determine communication, because they apply bit-not to rank
    // - prefix Y,Z determine node-wide coefficient, because they contain rank-determined !=1 elements
    // - suffix X,Y,Z determine local amp coefficients
    auto [targsX, targsY, targsZ] = paulis_getSeparateInds(str, qureg);
    auto [prefixX, suffixX] = util_getPrefixAndSuffixQubits(targsX, qureg);
    auto [prefixY, suffixY] = util_getPrefixAndSuffixQubits(targsY, qureg);
    auto [prefixZ, suffixZ] = util_getPrefixAndSuffixQubits(targsZ, qureg);

    // embarrassingly parallel when there is only Z's in prefix
    if (prefixX.empty() && prefixY.empty()) {

        // which has optimised local routines if suffix str is I, IZ, or IZXY
        if (suffixX.empty() && suffixY.empty() && suffixZ.empty())
            value = accel_statevec_calcTotalProb_sub(qureg);
        else if (suffixX.empty() && suffixY.empty())
            value = accel_statevec_calcExpecAnyTargZ_sub(qureg, suffixZ);
        else
            value = accel_statevec_calcExpecPauliStr_subA(qureg, suffixX, suffixY, suffixZ);

    // otherwise, communication is a pairwise exchange
    } else {
        auto prefixXY = util_getConcatenated(prefixX, prefixY);
        int pairRank = util_getRankWithQubitsFlipped(prefixXY, qureg);
        comm_exchangeAmpsToBuffers(qureg, pairRank);

        // and the subsequent process is optimal for all suffix paulis
        value = accel_statevec_calcExpecPauliStr_subB(qureg, suffixX, suffixY, suffixZ);
    }

    // apply this node's pre-factor
    value *= paulis_getPrefixPaulisElem(qureg, prefixY, prefixZ);

    // combine contributions from each node
    if (qureg.isDistributed)
        comm_reduceAmp(&value);

    return value;
}


qcomp localiser_densmatr_calcExpecPauliStr(Qureg qureg, PauliStr str) {
    assert_localiserGivenDensMatr(qureg);

    qcomp value = 0;

    // all ket-paulis are in the suffix state
    auto [targsX, targsY, targsZ] = paulis_getSeparateInds(str, qureg);

    // ergo all scenarios (I, IZ, IXYZ) are embarrassingly parallel
    if (targsX.empty() && targsY.empty() && targsZ.empty())
        value = accel_densmatr_calcTotalProb_sub(qureg);
    else if (targsX.empty() && targsY.empty())
        value = accel_densmatr_calcExpecAnyTargZ_sub(qureg, targsZ);
    else
        value = accel_densmatr_calcExpecPauliStr_sub(qureg, targsX, targsY, targsZ);

    // combine contributions from each node
    if (qureg.isDistributed)
        comm_reduceAmp(&value);

    return value;
}



/*
 * INNER PRODUCTS
 */


qcomp calcInnerProdOfSameDistribQuregs(Qureg quregA, Qureg quregB) {

    // both are statevectors/density-matrices, and both are local/distributed
    qcomp prod = accel_statevec_calcInnerProduct_sub(quregA, quregB);

    // unless quregs were cloned per-node, all nodes must sum their amp
    if (quregA.isDistributed)
        comm_reduceAmp(&prod);

    return prod;
}


Qureg spoofDistributedQuregFromLocalQureg(Qureg local, Qureg distrib) {

    // overwrite spoof's fields with distrib's, setting correct dimensions
    Qureg spoof = distrib;

    // set spoof's pointers to a rank-specific offset of local's
    qindex offset = util_getGlobalIndexOfFirstLocalAmp(distrib);
    spoof.cpuAmps = &local.cpuAmps[offset];
    spoof.gpuAmps = (local.isGpuAccelerated)? &local.gpuAmps[offset] : nullptr;

    return spoof;
}


qcomp localiser_statevec_calcInnerProduct(Qureg quregA, Qureg quregB) {

    // trivial when identically distributed
    if (quregA.isDistributed == quregB.isDistributed)
        return calcInnerProdOfSameDistribQuregs(quregA, quregB); // broadcasts

    // when only one qureg is distributed, we spoof a distributed
    // qureg from the local one, offsetting its amp pointers. In
    // defensive design, we modify a copy of quregA or quregB (even
    // though they themselves are mere copies of the user structs),
    // in case this function is modified in the future to mutate args
    Qureg copyA = (quregA.isDistributed)? quregA : spoofDistributedQuregFromLocalQureg(quregA, quregB);
    Qureg copyB = (quregB.isDistributed)? quregB : spoofDistributedQuregFromLocalQureg(quregB, quregA);

    return calcInnerProdOfSameDistribQuregs(copyA, copyB); // broadcasts
}


qcomp localiser_densmatr_calcFidelityWithPureState(Qureg rho, Qureg psi, bool conj) {
    assert_localiserGivenDensMatr(rho);
    assert_localiserGivenStateVec(psi);

    // each node will first compute their local fidelity contribution
    qcomp fid = 0;
    
    // rho and psi may have different distributions, though we
    // ultimately require psi to be duplicated on every node
    bool rhoDist = rho.isDistributed;
    bool psiDist = psi.isDistributed;

    // when psi is duplicated on every node, local eval is trivial
    if (!psiDist)
        fid = accel_densmatr_calcFidelityWithPureState_sub(rho, psi, conj);

    // psi distributed but rho duplicated is illegal
    else if (!rhoDist)
        error_calcFidStateVecDistribWhileDensMatrLocal();

    // when both are distributed...
    else {
        // we broadcast psi to every node's rho comm buffer...
        comm_combineAmpsIntoBuffer(rho, psi);

        // and spoof non-distributed statevector
        int isDense = 0;
        int isDistrib = 0;
        Qureg spoof = qureg_populateNonHeapFields(
            psi.numQubits, isDense, isDistrib, 
            rho.isGpuAccelerated, rho.isMultithreaded);

        // to which we bind rho's buffers containing psi's amps
        spoof.cpuAmps = rho.cpuCommBuffer;
        spoof.gpuAmps = rho.gpuCommBuffer;
        fid = accel_densmatr_calcFidelityWithPureState_sub(rho, spoof, conj);
    }

    // if rho was distributed, combine node contributions
    if (rho.isDistributed)
        comm_reduceAmp(&fid);

    return fid;
}


qreal localiser_densmatr_calcHilbertSchmidtDistance(Qureg quregA, Qureg quregB) {
    assert_localiserGivenDensMatr(quregA);
    assert_localiserGivenDensMatr(quregB);

    qreal dist = 0;

    // when distributions match, routine is embarrassingly parallel
    if (quregA.isDistributed == quregB.isDistributed)
        dist = accel_densmatr_calcHilbertSchmidtDistance_sub(quregA, quregB);

    // otherwise, we simply spoof a distributed qureg from the local 
    // one, offsetting its amp pointers, modifying copies in defensive design
    else {
        Qureg copyA = (quregA.isDistributed)? quregA : spoofDistributedQuregFromLocalQureg(quregA, quregB);
        Qureg copyB = (quregB.isDistributed)? quregB : spoofDistributedQuregFromLocalQureg(quregB, quregA);
        dist = accel_densmatr_calcHilbertSchmidtDistance_sub(copyA, copyB);
    }

    // combine node contributions unless both quregs were duplicated
    if (quregA.isDistributed || quregB.isDistributed)
        comm_reduceReal(&dist);

    return dist; // no sqrt
}



/*
 * PROJECTORS 
 */


void localiser_statevec_multiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {
    assert_localiserGivenStateVec(qureg);

    // we pass all qubits (including prefixes) to backend, which enumerates every
    // local amp and decides whether to multiply 0 or a renormalisation thereupon.
    // In distributed settings, it may be that some nodes zero their entire local
    // partition, while others renormalise; in theory, we could make bespoke routines
    // for these scenarios, eliminating superfluous flops. This is an unworthwhile
    // optimisation; there would very likely remain nodes that must do zero-and-renorm
    // iteration, for which all other nodes would later have to wait at next sync.

    // always embarrassingly parallel
    accel_statevec_multiQubitProjector_sub(qureg, qubits, outcomes, prob);
}


void localiser_densmatr_multiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {
    assert_localiserGivenDensMatr(qureg);

    // always embarrassingly parallel
    accel_densmatr_multiQubitProjector_sub(qureg, qubits, outcomes, prob);
}
