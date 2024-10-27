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

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <tuple>
#include <vector>
#include <complex>
#include <algorithm>

using std::vector;



/*
 * PRIVATE FUNCTIONS
 */


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


bool doAnyLocalAmpsSatisfyCtrls(Qureg qureg, vector<int> ctrls, vector<int> states) {

    // this answers the generic question of "do any of the given qubits lie in the
    // prefix substate with node-fixed values inconsistent with the given states?"

    // non-distributed quregs always have amps satisfying ctrls
    if (!qureg.isDistributed)
        return true;

    // check each ctrl qubit
    for (size_t i=0; i<ctrls.size(); i++) {

        // consider only ctrls which operate on the prefix substate
        if (util_isQubitInSuffix(ctrls[i], qureg))
            continue;

        // abort if any prefix ctrl has wrong bit value
        if (util_getRankBitOfQubit(ctrls[i], qureg) != states[i])
            return false;
    }

    // otherwise all prefix ctrls have the specified values
    return true;
}


auto getPrefixOrSuffixQubits(Qureg qureg, vector<int> qubits, bool getSuffix) {

    vector<int> subQubits(0);
    subQubits.reserve(qubits.size());

    for (int qubit : qubits)
        if (util_isQubitInSuffix(qubit, qureg) == getSuffix)
            subQubits.push_back(qubit);

    return subQubits;
}

auto getSuffixQubits(Qureg qureg, vector<int> qubits) {
    return getPrefixOrSuffixQubits(qureg, qubits, true);
}

auto getPrefixQubits(Qureg qureg, vector<int> qubits) {
    return getPrefixOrSuffixQubits(qureg, qubits, false);
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
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
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
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
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
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
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
    if (doAnyLocalAmpsSatisfyCtrls(qureg, newCtrls, ctrlStates)) {

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
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
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
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
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
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
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

    qindex offset = qureg.rank * qureg.numAmpsPerNode;

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
 * PAULI TENSORS AND GADGETS
 */


extern int  paulis_getPauliAt(PauliStr str, int ind);
extern bool paulis_containsXOrY(PauliStr str);
extern vector<int> paulis_getSortedIndsOfNonIdentityPaulis(PauliStr str);
extern vector<int> paulis_getTargsWithEitherPaulis(vector<int> targs, PauliStr str, int pauliA, int pauliB);


void anyCtrlAnyTargZOrPhaseGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qcomp fac0, qcomp fac1) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);

    // operation is diagonal and always embarrasingly parallel, regardless of prefix targs
    accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(qureg, ctrls, ctrlStates, targs, fac0, fac1);
}


void anyCtrlPauliTensorOrGadgetOnPrefix(
    Qureg qureg, int pairRank, vector<int> ctrls, vector<int> ctrlStates, 
    qindex suffixMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
) {
    // pack and exchange all amps satisfying ctrl
    exchangeAmpsToBuffersWhereQubitsAreInStates(qureg, pairRank, ctrls, ctrlStates);

    // we cannot use suffixMaskXY to determine the buffer indices corresponding to the subsequently processed local amps,
    // because it operates on the full local amp index, whereas the buffer (potentially) excluded many amps (which did not
    // satisfy ctrls) and has a reduced index space. So we compute a new mask which operates on the buffer indices by
    // removing all ctrl qubits from the full index mask.
    auto sortedCtrls = util_getSorted(ctrls);
    auto bufferMaskXY = removeBits(suffixMaskXY, sortedCtrls.data(), sortedCtrls.size());

    accel_statevector_anyCtrlPauliTensorOrGadget_subB(qureg, ctrls, ctrlStates, suffixMaskXY, bufferMaskXY, allMaskYZ, powI, fac0, fac1);
}


void anyCtrlPauliTensorOrGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, PauliStr str, qcomp fac0, qcomp fac1) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to local amp modification
    removePrefixQubitsAndStates(qureg, ctrls, ctrlStates);

    // readable flags
    const int X=1, Y=2, Z=3;

    // total number of Y determines a phase factor on all updated amps (because Y contains i)
    int numY = paulis_getTargsWithEitherPaulis(targs, str, Y, Y).size();

    // precision-agnostic calculation of i^numY
    vector<qcomp> powersOfI = {1, 1_i, -1, -1_i};
    qcomp powI = powersOfI[numY % 4];

    // parity of Y and Z on all qubits determines phase factor on updated amps (because Y and Z contain -1)
    auto allTargsYZ = paulis_getTargsWithEitherPaulis(targs, str, Y, Z);
    auto allMaskYZ = getBitMask(allTargsYZ.data(), allTargsYZ.size());
    
    // X and Y on suffix qubits determine local amp movement (because X and Y are anti-diagonal)
    auto suffixTargsXY = getSuffixQubits(qureg, paulis_getTargsWithEitherPaulis(targs, str, X, Y)); // sorted
    auto suffixMaskXY = getBitMask(suffixTargsXY.data(), suffixTargsXY.size());

    // X and Y on prefix qubits determines pair rank (because X and Y are anti-diagonal)
    auto prefixTargsXY = getPrefixQubits(qureg, paulis_getTargsWithEitherPaulis(targs, str, X, Y));
    auto pairRank = flipBits(qureg.rank, prefixTargsXY.data(), prefixTargsXY.size());

    // call embarrassingly parallel routine or induce pairwise communication
    (qureg.rank == pairRank)?
        accel_statevector_anyCtrlPauliTensorOrGadget_subA(
            qureg, ctrls, ctrlStates, suffixTargsXY, suffixMaskXY, allMaskYZ, powI, fac0, fac1):
        anyCtrlPauliTensorOrGadgetOnPrefix(
            qureg, pairRank, ctrls, ctrlStates, suffixMaskXY, allMaskYZ, powI, fac0, fac1);
}


void localiser_statevec_anyCtrlAnyTargZ(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs) {

    qcomp fac0 = qcomp(+1, 0); // even parity
    qcomp fac1 = qcomp(-1, 0); // odd parity
    anyCtrlAnyTargZOrPhaseGadget(qureg, ctrls, ctrlStates, targs, fac0, fac1);
}


void localiser_statevec_anyCtrlPhaseGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qreal angle) {

    qcomp fac0 = qcomp(cos(angle), +sin(angle)); // exp( i angle)
    qcomp fac1 = qcomp(cos(angle), -sin(angle)); // exp(-i angle)
    anyCtrlAnyTargZOrPhaseGadget(qureg, ctrls, ctrlStates, targs, fac0, fac1);
}


void localiser_statevec_anyCtrlPauliTensor(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str) {

    qcomp fac0 = 0; // even parity
    qcomp fac1 = 1; // odd parity
    auto targs = paulis_getSortedIndsOfNonIdentityPaulis(str);

    (paulis_containsXOrY(str))?
        anyCtrlPauliTensorOrGadget(qureg, ctrls, ctrlStates, targs, str, fac0, fac1):
        anyCtrlAnyTargZOrPhaseGadget(qureg, ctrls, ctrlStates, targs, fac0, fac1);
}


void localiser_statevec_anyCtrlPauliGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qreal angle) {

    qcomp fac0 = qcomp(cos(angle), 0); // even parity
    qcomp fac1 = qcomp(0, sin(angle)); // odd parity
    auto targs = paulis_getSortedIndsOfNonIdentityPaulis(str);

    (paulis_containsXOrY(str))?
        anyCtrlPauliTensorOrGadget(qureg, ctrls, ctrlStates, targs, str, fac0, fac1):
        anyCtrlAnyTargZOrPhaseGadget(qureg, ctrls, ctrlStates, targs, fac0, fac1);
}



/*
 * QUREG COMBINATION
 */


void mixDensityMatrixWithStatevector(qreal outProb, Qureg out, qreal inProb, Qureg in) {

    // we can handle 3 out of 4 possible combinations of distribution,
    // and accelerator.hpp will handle every combination of GPU-accel
    bool outDist = out.isDistributed;
    bool inDist = in.isDistributed;

    // illegal to distribute only the smaller Qureg ('out' has no buffer space to receive it)
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


void localiser_densmatr_oneQubitPauliChannel(Qureg qureg, int qubit, qreal probI, qreal probX, qreal probY, qreal probZ) {
    assert_localiserGivenDensMatr(qureg);

    (doesChannelRequireComm(qureg, qubit))?
        oneQubitPauliChannelOnPrefix(qureg, qubit, probI, probX, probY, probZ):
        accel_densmatr_oneQubitPauliChannel_subA(qureg, qubit, probI, probX, probY, probZ);
}


// twoQubitPauliChannel() is regrettably too difficult; the communication model cannot be 
// simplified the way it was in twoQubitDepolarising() which levereaged the uniform
// coefficients. It is not clear whether arbitrary coefficients, which cause many more
// amplitudes to mix, can eve be performed in a sequence of pairwise communication



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
 * GENERAL CHANNELS
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
