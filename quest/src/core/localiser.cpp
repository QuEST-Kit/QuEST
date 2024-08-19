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
#include "quest/src/comm/comm_routines.hpp"

#include <tuple>
#include <vector>
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
        if (targ >= qureg.logNumAmpsPerNode)
            return true;

    // sufix qubit targets need no communication
    return false;
}

bool doesGateRequireComm(Qureg qureg, int targ) {

    return doesGateRequireComm(qureg, vector{targ});
}


bool doAnyLocalAmpsSatisfyCtrls(Qureg qureg, vector<int> ctrls, vector<int> states) {

    // non-distributed quregs always have amps satisfying ctrls
    if (!qureg.isDistributed)
        return true;

    // check each ctrl qubit
    for (size_t i=0; i<ctrls.size(); i++) {

        // consider only ctrls which operate on the prefix substate
        if (ctrls[i] < qureg.logNumAmpsPerNode)
            continue;

        // abort if any prefix ctrl has wrong bit value
        int prefBit = getBit(qureg.rank, ctrls[i] - qureg.logNumAmpsPerNode);
        if (prefBit != states[i])
            return false;
    }

    // otherwise all prefix ctrls have the specified values
    return true;
}


auto getSuffixCtrlsAndStates(Qureg qureg, vector<int> ctrls, vector<int> states) {

    vector<int> suffixCtrls(0);   suffixCtrls .reserve(ctrls.size());
    vector<int> suffixStates(0);  suffixStates.reserve(states.size());

    for (size_t i=0; i<ctrls.size(); i++)
        if (ctrls[i] < qureg.logNumAmpsPerNode) {
            suffixCtrls.push_back(ctrls[i]);
            suffixStates.push_back(states[i]);
        }

    // return new vectors, leaving old unmodified
    return std::tuple{suffixCtrls, suffixStates};
}


auto getCtrlsAndTargsSwappedToSuffix(Qureg qureg, vector<int> ctrls, vector<int> targs) {

    // nothing to do if all targs are already in suffix
    if (!doesGateRequireComm(qureg, targs))
        return std::tuple{ctrls, targs};

    // prepare masks to avoid quadratic nested looping
    qindex targMask = getBitMask(targs.data(), targs.size());
    qindex ctrlMask = getBitMask(ctrls.data(), ctrls.size());
    int minNonTarg = getIndOfNextRightmostZeroBit(targMask, -1);

    // prepare indices of ctrls in the given list (i.e. the inverse of ctrls)
    int maxCtrlInd = *std::max_element(ctrls.begin(), ctrls.end());
    vector<int> ctrlInds(maxCtrlInd+1); // bounded by ~64
    for (size_t i=0; i<ctrls.size(); i++)
        ctrlInds[ctrls[i]] = i;

    // check every target in arbitrary order
    for (size_t i=0; i<targs.size(); i++) {

        // consider only targs in the prefix substate
        if (!doesGateRequireComm(qureg, targs[i]))
            continue;
            
        // if our replacement targ happens to be a ctrl... 
        if (getBit(ctrlMask, minNonTarg)) {

            // find and swap that ctrl with the old targ
            int ctrlInd = ctrlInds[minNonTarg];
            ctrls[ctrlInd] = targs[i];

            // update our ctrl trackers
            ctrlInds[targs[i]] = ctrlInd;
            ctrlInds[minNonTarg] = -1; // for clarity
            ctrlMask = flipTwoBits(ctrlMask, minNonTarg, targs[i]);
        }

        // swap the prefix targ with the smallest available suffix targ
        targs[i] = minNonTarg;

        // update our targ trackers
        targMask = flipTwoBits(targMask, targs[i], minNonTarg);
        minNonTarg = getIndOfNextRightmostZeroBit(targMask, minNonTarg);
    }

    // the ordering in ctrls relative to the caller's ctrlStates is unchanged
    return std::tuple{ctrls, targs};
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
 * MATRICES
 */


// DEBUG
void localiser_statevec_anyCtrlSwap(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {
}


void localiser_statevec_anyCtrlOneTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits as relevant to communication and local amp modification
    std::tie(ctrls, ctrlStates) = getSuffixCtrlsAndStates(qureg, ctrls, ctrlStates);

    // if the target permits embarrassingly parallel simulation, perform it and finish
    if (!doesGateRequireComm(qureg, targ)) {
        accel_statevec_anyCtrlOneTargDenseMatr_subA(qureg, ctrls, ctrlStates, targ, matr);
        return;
    }

    // otherwise we must communicate with a pair rank... 
    int rankTarg = targ - qureg.logNumAmpsPerNode;
    int pairRank = flipBit(qureg.rank, rankTarg);

    // to exchange all or some of our amps (those where ctrls are active) into buffer
    if (ctrls.empty())
        comm_exchangeAmpsToBuffers(qureg, pairRank);
    else {
        qindex numExch = qureg.numAmpsPerNode / powerOf2(ctrls.size());
        accel_statevec_packAmpsIntoBuffer(qureg, ctrls, ctrlStates);
        comm_exchangeBuffers(qureg, numExch, pairRank);
    }

    // extract relevant gate elements
    int bit = getBit(qureg.rank, rankTarg);
    qcomp fac0 = matr.elems[bit][ bit];
    qcomp fac1 = matr.elems[bit][!bit];

    // update local amps using received amps in buffer
    accel_statevec_anyCtrlOneTargDenseMatr_subB(qureg, ctrls, ctrlStates, fac0, fac1);
}


void localiser_statevec_anyCtrlAnyTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);
    
    // TODO:
    //   - the sequence of pair-wise full-swaps should be more efficient as a
    //     "single" sequence of smaller messages sending amps directly to their
    //     final destination node. This could use a new "multiSwap" function.
    //   - if the user has compiled cuQuantum, and Qureg is GPU-accelerated, the
    //     multiSwap function should use custatevecSwapIndexBits() if local,
    //     or custatevecDistIndexBitSwapSchedulerSetIndexBitSwaps() if distributed,
    //     although the latter requires substantially more work like setting up
    //     a communicator which may be inelegant alongside our own distribution scheme
    //   - the induced swap gates always trigger the targ2-prefix targ1-suffix scenario; 
    //     we should invoke that scenario directly to avoid all the extra faff

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // if all targets lie within the suffix node...
    if (!doesGateRequireComm(qureg, targs)) {

        // retain only suffix controls, and perform embarrassingly parallel simulation, then finish
        std::tie(ctrls, ctrlStates) = getSuffixCtrlsAndStates(qureg, ctrls, ctrlStates);
        accel_statevec_anyCtrlAnyTargDenseMatr_sub(qureg, ctrls, ctrlStates, targs, matr);
        return;
    }

    // else, find suffix positions for all prefix targs, moving colliding ctrls out of the way
    auto [newCtrls, newTargs] = getCtrlsAndTargsSwappedToSuffix(qureg, ctrls, targs);

    // only unmoved ctrls can be applied to the swaps, to accelerate them
    auto [unmovedCtrls, unmovedCtrlStates] = getNonSwappedCtrlsAndStates(newCtrls, ctrlStates, newTargs); 

    // perform necessary swaps to move all targets into suffix, each of which invokes communication
    for (size_t i=0; i<targs.size(); i++)
        if (targs[i] != newTargs[i])
            localiser_statevec_anyCtrlSwap(qureg, unmovedCtrls, unmovedCtrlStates, targs[i], newTargs[i]);

    // if the moved ctrls do not eliminate the need for local simulation...
    if (doAnyLocalAmpsSatisfyCtrls(qureg, newCtrls, ctrlStates)) {

        // then perform embarrassingly parallel simulation using only the new suffix ctrls
        auto [newSuffixCtrls, newSuffixStates] = getSuffixCtrlsAndStates(qureg, newCtrls, ctrlStates);
        accel_statevec_anyCtrlAnyTargDenseMatr_sub(qureg, newSuffixCtrls, newSuffixStates, newTargs, matr);
    }

    // undo swaps, each invoking communication
    for (size_t i=0; i<targs.size(); i++)
        if (targs[i] != newTargs[i])
            localiser_statevec_anyCtrlSwap(qureg, unmovedCtrls, unmovedCtrlStates, targs[i], newTargs[i]);
}


void localiser_statevec_anyCtrlAnyTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to local amp modification
    std::tie(ctrls, ctrlStates) = getSuffixCtrlsAndStates(qureg, ctrls, ctrlStates);
    
    // diagonal matrices are always embarrassingly parallel, regardless of whether any targs are in prefix
    return accel_statevec_anyCtrlAnyTargDiagMatr_sub(qureg, ctrls, ctrlStates, targs, matr);
}
