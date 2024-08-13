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

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/indexer.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <tuple>
#include <vector>

using std::vector;



/*
 * PRIVATE CONVENIENCE FUNCTIONS
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



/*
 * DETERMINING NECESSARY COMMUNICATION
 */


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

    // states will be empty or the same length as ctrls
    vector<int> suffixCtrls(0);
    vector<int> suffixStates(0);

    // lovely premature optimisation
    suffixCtrls.reserve(ctrls.size());
    suffixStates.reserve(states.size());

    for (size_t i=0; i<ctrls.size(); i++)
        if (ctrls[i] < qureg.logNumAmpsPerNode) {
            suffixCtrls.push_back(ctrls[i]);
            suffixStates.push_back(states[i]);
        }

    // return new vectors, leaving old unmodified
    return std::tuple{suffixCtrls, suffixStates};
}



/*
 * MATRICES
 */


void statevec_anyCtrlOneTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
    assertValidCtrlStates(ctrls, ctrlStates);
    setDefaultCtrlStates(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits as relevant to communication and local amp modification
    std::tie(ctrls, ctrlStates) = getSuffixCtrlsAndStates(qureg, ctrls, ctrlStates);

    // if the target permits embarrassingly parallel simulation, perform it and finish
    if (!doesGateRequireComm(qureg, {targ})) {
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


void statevec_anyCtrlManyTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr matr) {

    // TODO
    //  we nmay not even need above bespoke method - you'll have to benchmark!
}


void statevec_anyCtrlAnyTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {
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