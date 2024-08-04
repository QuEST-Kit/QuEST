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

using std::tuple;
using std::vector;
using namespace index_flags;



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
    for (int i=0; i<ctrls.size(); i++) {

        // consider only ctrls which operate on the prefix substate
        if (ctrls[i] < qureg.logNumAmpsPerNode)
            continue;

        // compare prefix ctrl to its goal bit (1 if unspecified)
        int prefBit = getBit(qureg.rank, ctrls[i] - qureg.logNumAmpsPerNode);
        int goalBit = (states.empty())? 1 : states[i];

        // abort if any prefix ctrl has wrong bit value
        if (prefBit != goalBit)
            return false;
    }

    // otherwise all prefix ctrls have the specified values
    return true;
}


tuple<vector<int>,vector<int>> getSuffixCtrls(Qureg qureg, vector<int> ctrls, vector<int> states) {

    // states will be empty or the same length as ctrls
    vector<int> suffixCtrls(0);
    vector<int> suffixStates(0);

    // lovely premature optimisation
    suffixCtrls.reserve(ctrls.size());
    suffixStates.reserve(states.size());

    for (int i=0; i<ctrls.size(); i++)
        if (ctrls[i] < qureg.logNumAmpsPerNode) {
            suffixCtrls.push_back(ctrls[i]);

            if (!states.empty())
                suffixStates.push_back(states[i]);
        }

    // return new vectors, leaving old unmodified
    return {suffixCtrls, suffixStates};
}



/*
 * ANY-CTRL ONE-TARG MATRIX
 */


template <class MatrType>
void inner_statevec_anyCtrlOneTargMatrix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, MatrType matr) {
    indexer_assertValidCtrls(ctrls, ctrlStates);

    // node has nothing to do if all local amps violate control condition
    if (!doAnyLocalAmpsSatisfyCtrls(qureg, ctrls, ctrlStates))
        return;

    // retain only suffix control qubits, as relevant to communication and local amp modification
    std::tie(ctrls, ctrlStates) = getSuffixCtrls(qureg, ctrls, ctrlStates);

    // if the target permits embarrassingly parallel simulation, perform it and finish
    if (!doesGateRequireComm(qureg, {targ}))
        return accel_statevec_anyCtrlOneTargMatrix_subA(qureg, ctrls, ctrlStates, targ, matr);
    
    // diagonal matrices are always embarrassingly parallel, so perform and finish
    if constexpr (util_isDiagonalMatrixType<MatrType>())
        return accel_statevec_anyCtrlOneTargMatrix_subA(qureg, ctrls, ctrlStates, targ, matr);

    // but for dense matrices, we must communicate with a pair rank... 
    int rankTarg = targ - qureg.logNumAmpsPerNode;
    int pairRank = flipBit(qureg.rank, rankTarg);

    // to exchange all or some of our amps (those where ctrls are active) into buffer
    if (ctrls.empty())
        comm_exchangeAmpsToBuffers(qureg, pairRank);
    else {
        qindex numPacked = accel_statevec_packAmpsIntoBuffer(qureg, ctrls, ctrlStates);
        comm_exchangeBuffers(qureg, numPacked, pairRank);
    }

    // extract relevant gate elements (but don't let the compiler see the unreachable diagonal matrix access)
    qcomp fac0 = 0;
    qcomp fac1 = 0;
    if constexpr (util_isDenseMatrixType<MatrType>()) {
        int bit = getBit(qureg.rank, rankTarg);
        qcomp fac0 = matr.elems[bit][ bit];
        qcomp fac1 = matr.elems[bit][!bit]; // compilers hate him
    }

    // update local amps using received amps in buffer
    accel_statevec_anyCtrlOneTargDenseMatrix_subB(qureg, ctrls, ctrlStates, fac0, fac1);
}

void statevec_anyCtrlOneTargMatrix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    inner_statevec_anyCtrlOneTargMatrix(qureg, ctrls, ctrlStates, targ, matr);
}

void statevec_anyCtrlOneTargMatrix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {

    inner_statevec_anyCtrlOneTargMatrix(qureg, ctrls, ctrlStates, targ, matr);
}
