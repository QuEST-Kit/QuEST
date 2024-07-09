/** @file
 * Internal functions which localize the data needed for simulation.
 * That is, they determine whether performing a simulation requires
 * Querg amplitudes from other distributed nodes and if so, invoke
 * the necessary communication, before finally calling the 
 * embarrassingly parallel subroutines in accelerator.cpp. This is
 * done agnostically of whether amplitudes of the Qureg are being
 * stored in RAM (CPU) or VRAM (GPU).
 */

#include "qureg.h"
#include "structures.h"

#include "../comm/comm_routines.hpp"
#include "../core/bitwise.hpp"
#include "../core/accelerator.hpp"



/*
 * PRIVATE
 */


bool oneTargetGateRequiresComm(Qureg qureg, int target) {
    if (!qureg.isDistributed)
        return false;
    
    return target >= qureg.logNumAmpsPerNode;
}



/*
 * OPERATORS
 */


void statevec_oneTargetGate(Qureg qureg, int target, CompMatr1 matrix) {

    if (!oneTargetGateRequiresComm(qureg, target))
        statevec_oneTargetGate_subA(qureg, target, matrix);

    else {
        // exchange all amps (receive to buffer)
        int rankTarget = target - qureg.logNumAmpsPerNode;
        int pairRank = flipBit(qureg.rank, rankTarget);

        comm_exchangeAmpsToBuffers(qureg, pairRank);

        // extract relevant gate elements
        int bit = getBit(qureg.rank, rankTarget);
        qcomp fac0 = matrix.elems[bit][bit];
        qcomp fac1 = matrix.elems[bit][!bit];

        statevec_oneTargetGate_subB(qureg, fac0, fac1);
    }
}
