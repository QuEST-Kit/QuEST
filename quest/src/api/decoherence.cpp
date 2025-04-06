/** @file
 * API definitions for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 * 
 * @author Tyson Jones
 * @author Balint Koczor (prototyped v3 mixKrausMap)
 * @author Nicolas Vogt (prototyped v3 mixDamping)
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/channels.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/utilities.hpp"

#include <vector>
using std::vector;



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

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


void mixSuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_superOpFields(superop, __func__);
    validate_superOpIsSynced(superop, __func__);
    validate_superOpDimMatchesTargs(superop, numTargets, __func__);
    validate_mixedAmpsFitInNode(qureg, 2*numTargets, __func__); // superop acts on 2x

    localiser_densmatr_superoperator(qureg, superop, util_getVector(targets, numTargets));
}



} // end de-mangler



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

void mixKrausMap(Qureg qureg, vector<int> targets, KrausMap map) {
    mixKrausMap(qureg, targets.data(), targets.size(), map);
}

void mixSuperOp(Qureg qureg, vector<int> targets, SuperOp superop) {
    mixSuperOp(qureg, targets.data(), targets.size(), superop);
}

#endif // __cplusplus
