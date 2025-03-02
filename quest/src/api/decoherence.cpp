/** @file
 * API definitions for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/channels.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/utilities.hpp"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions

// enable invocation by both C and C++ binaries
extern "C" {



void mixDephasing(Qureg qureg, int qubit, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDepashingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitDephasing(qureg, qubit, prob);
}


void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);
    validate_twoQubitDepashingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_twoQubitDephasing(qureg, qubit1, qubit2, prob);
}


void mixDepolarising(Qureg qureg, int qubit, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDepolarisingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitDepolarising(qureg, qubit, prob);
}


void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);
    validate_twoQubitDepolarisingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_twoQubitDepolarising(qureg, qubit1, qubit2, prob);
}


void mixDamping(Qureg qureg, int qubit, qreal prob) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitDampingProb(prob, __func__);

    // permit but do not change non-decohering statevecs
    if (prob == 0) 
        return;
    
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitDamping(qureg, qubit, prob);
}


void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_oneQubitPauliChannelProbs(probX, probY, probZ, __func__);

    // permit but do not change non-decohering statevecs
    if (probX == 0 && probY == 0 && probZ == 0)
        return;

    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_oneQubitPauliChannel(qureg, qubit, probX, probY, probZ);
}


void mixKrausMap(Qureg qureg, int* qubits, int numQubits, KrausMap map) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_mixedAmpsFitInNode(qureg, 2*numQubits, __func__); // superop acts on 2x
    validate_krausMapIsCPTP(map, __func__); // also checks fields and is-sync
    validate_krausMapMatchesTargets(map, numQubits, __func__);

    localiser_densmatr_krausMap(qureg, map, util_getVector(qubits, numQubits));
}


void mixQureg(Qureg outQureg, Qureg inQureg, qreal inProb) {
    validate_quregFields(outQureg, __func__);
    validate_quregFields(inQureg, __func__);
    validate_probability(inProb, __func__);
    validate_quregsCanBeMixed(outQureg, inQureg, __func__); // checks outQureg is densmatr

    qreal outProb = 1 - inProb;
    localiser_densmatr_mixQureg(outProb, outQureg, inProb, inQureg);
}



} // end de-mangler