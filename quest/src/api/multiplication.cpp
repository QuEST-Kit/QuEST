/** @file
 * API definitions for directly pre- and post-multiplying 
 * operators upon density matrices, likely constituting 
 * non-physical operations which break state normalisation.
 * 
 * @author Tyson Jones
 */

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/multiplication.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/paulilogic.hpp"

#include <vector>

using std::vector;



/*
 * CompMatr1
 */

extern "C" {

void leftapplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    bool transp = false;
    localiser_statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, target, matrix, conj, transp);
}

void rightapplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__);
    
    // rho matrix ~ transpose(rho) (x) I ||rho>>
    bool conj = false;
    bool transp = true;
    int qubit = util_getBraQubit(target, qureg);
    localiser_statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, qubit, matrix, conj, transp);
}

} // end de-mangler



/*
 * CompMatr2
 */

extern "C" {

void leftapplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__);
    validate_mixedAmpsFitInNode(qureg, 2, __func__);

    bool conj = false;
    bool transp = false;
    localiser_statevec_anyCtrlTwoTargDenseMatr(qureg, {}, {}, target1, target2, matrix, conj, transp);
}

void rightapplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__);
    validate_mixedAmpsFitInNode(qureg, 2, __func__);

    // rho matrix ~ transpose(rho) (x) I ||rho>>
    bool conj = false;
    bool transp = true;
    int qubit1 = util_getBraQubit(target1, qureg);
    int qubit2 = util_getBraQubit(target2, qureg);
    localiser_statevec_anyCtrlTwoTargDenseMatr(qureg, {}, {}, qubit1, qubit2, matrix, conj, transp);
}

} // end de-mangler



/*
 * CompMatr
 */

extern "C" {

void leftapplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync
    validate_mixedAmpsFitInNode(qureg, numTargets, __func__);

    bool conj = false;
    bool transp = false;
    localiser_statevec_anyCtrlAnyTargDenseMatr(qureg, {}, {}, util_getVector(targets, numTargets), matrix, conj, transp);
}

void rightapplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync
    validate_mixedAmpsFitInNode(qureg, numTargets, __func__);

    // rho matrix ~ transpose(rho) (x) I ||rho>>
    bool conj = false;
    bool transp = true;
    auto qubits = util_getBraQubits(util_getVector(targets, numTargets), qureg);
    localiser_statevec_anyCtrlAnyTargDenseMatr(qureg, {}, {}, qubits, matrix, conj, transp);
}

} // end de-mangler

void leftapplyCompMatr(Qureg qureg, vector<int> targets, CompMatr matr) {

    leftapplyCompMatr(qureg, targets.data(), targets.size(), matr);
}

void rightapplyCompMatr(Qureg qureg, vector<int> targets, CompMatr matr) {

    rightapplyCompMatr(qureg, targets.data(), targets.size(), matr);
}



/*
 * DiagMatr1
 */

extern "C" {

void leftapplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    localiser_statevec_anyCtrlOneTargDiagMatr(qureg, {}, {}, target, matrix, conj);
}

void rightapplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    int qubit = util_getBraQubit(target, qureg);
    localiser_statevec_anyCtrlOneTargDiagMatr(qureg, {}, {}, qubit, matrix, conj);
}

} // end de-mangler



/*
 * DiagMatr2
 */

extern "C" {

void leftapplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    localiser_statevec_anyCtrlTwoTargDiagMatr(qureg, {}, {}, target1, target2, matrix, conj);
}

void rightapplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    int qubit1 = util_getBraQubit(target1, qureg);
    int qubit2 = util_getBraQubit(target2, qureg);
    localiser_statevec_anyCtrlTwoTargDiagMatr(qureg, {}, {}, qubit1, qubit2, matrix, conj);
}

} // end de-mangler



/*
 * DiagMatr
 */

extern "C" {

void leftapplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync

    bool conj = false;
    qcomp exponent = 1;
    auto qubits = util_getVector(targets, numTargets);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, qubits, matrix, exponent, conj);
}

void rightapplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync

    bool conj = false;
    qcomp exponent = 1;
    auto qubits = util_getBraQubits(util_getVector(targets, numTargets), qureg);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, qubits, matrix, exponent, conj);
}

} // end de-mangler

void leftapplyDiagMatr(Qureg qureg, vector<int> targets, DiagMatr matrix) {

    leftapplyDiagMatr(qureg, targets.data(), targets.size(), matrix);
}

void rightapplyDiagMatr(Qureg qureg, vector<int> targets, DiagMatr matrix) {

    rightapplyDiagMatr(qureg, targets.data(), targets.size(), matrix);
}



/*
 * DiagMatrPower
 */

extern "C" {

void leftapplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync, but not unitarity
    validate_matrixExpIsNonDiverging(matrix, exponent, __func__); // harmlessly re-validates fields and is-sync

    bool conj = false;
    auto qubits = util_getVector(targets, numTargets);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, qubits, matrix, exponent, conj);
}

void rightapplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync, but not unitarity
    validate_matrixExpIsNonDiverging(matrix, exponent, __func__); // harmlessly re-validates fields and is-sync

    bool conj = false;
    auto qubits = util_getBraQubits(util_getVector(targets, numTargets), qureg);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, qubits, matrix, exponent, conj);
}

} // end de-mangler

void leftapplyDiagMatrPower(Qureg qureg, vector<int> targets, DiagMatr matrix, qcomp exponent) {

    leftapplyDiagMatrPower(qureg, targets.data(), targets.size(), matrix, exponent);
}

void rightapplyDiagMatrPower(Qureg qureg, vector<int> targets, DiagMatr matrix, qcomp exponent) {

    rightapplyDiagMatrPower(qureg, targets.data(), targets.size(), matrix, exponent);
}



/*
 * FullStateDiagMatr (and power)
 */

extern "C" {

void leftapplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__); // matrix can be non-unitary

    leftapplyFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}

void leftapplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__); // matrix can be non-unitary
    validate_matrixExpIsNonDiverging(matrix, exponent, __func__);

    // rho -> matrix^exponent rho
    bool leftMultiply = true;
    bool rightMultiply = false;
    bool rightConj = false;

    (qureg.isDensityMatrix)?
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, leftMultiply, rightMultiply, rightConj):
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent);
}

void rightapplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__); // matrix can be non-unitary

    rightapplyFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}

void rightapplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__); // matrix can be non-unitary
    validate_matrixExpIsNonDiverging(matrix, exponent, __func__);

    // rho -> rho matrix^exponent
    bool leftMultiply = false;
    bool rightMultiply = true;
    bool rightConj = false;
    localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, leftMultiply, rightMultiply, rightConj);
}

} // end de-mangler



/*
 * swap
 */

extern "C" {

void leftapplySwap(Qureg qureg, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);

    localiser_statevec_anyCtrlSwap(qureg, {}, {}, qubit1, qubit2);
}

void rightapplySwap(Qureg qureg, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);

    qubit1 = util_getBraQubit(qubit1, qureg);
    qubit2 = util_getBraQubit(qubit2, qureg);
    localiser_statevec_anyCtrlSwap(qureg, {}, {}, qubit1, qubit2);
}

} // end de-mangler



/*
 * individual Paulis
 */

extern "C" {

void leftapplyPauliX(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("X", {target});
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void leftapplyPauliY(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("Y", {target});
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void leftapplyPauliZ(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("Z", {target});
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void rightapplyPauliX(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("X", {target});
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void rightapplyPauliY(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);

    qcomp factor = -1; // undo transpose
    PauliStr str = getPauliStr("Y", {target});
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str, factor);
}

void rightapplyPauliZ(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("Z", {target});
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

} // end de-mangler



/*
 * Pauli strings
 */

extern "C" {

void leftapplyPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void rightapplyPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qcomp factor = paulis_getSignOfPauliStrConj(str); // undo transpose
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str, factor);
}

} // end de-mangler



/*
 * Pauli gadgets
 */

extern "C" {

void leftapplyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    localiser_statevec_anyCtrlPauliGadget(qureg, {}, {}, str, phase);
}

void rightapplyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qreal factor = paulis_getSignOfPauliStrConj(str);
    qreal phase = factor * util_getPhaseFromGateAngle(angle);
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliGadget(qureg, {}, {}, str, phase);
}

} // end de-mangler



/*
 * phase gadgets
 */

extern "C" {

void leftapplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    auto qubits = util_getVector(targets, numTargets);
    localiser_statevec_anyCtrlPhaseGadget(qureg, {}, {}, qubits, phase);
}

void rightapplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    auto qubits = util_getBraQubits(util_getVector(targets, numTargets), qureg);
    localiser_statevec_anyCtrlPhaseGadget(qureg, {}, {}, qubits, phase);
}

} // end de-mangler

void leftapplyPhaseGadget(Qureg qureg, vector<int> targets, qreal angle) {

    leftapplyPhaseGadget(qureg, targets.data(), targets.size(), angle);
}

void rightapplyPhaseGadget(Qureg qureg, vector<int> targets, qreal angle) {

    rightapplyPhaseGadget(qureg, targets.data(), targets.size(), angle);
}



/*
 * many-qubit NOTs
 */

extern "C" {

void leftapplyMultiQubitNot(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    PauliStr str = getPauliStr(std::string(numTargets, 'X'), targets, numTargets);
    leftapplyPauliStr(qureg, str);
}

void rightapplyMultiQubitNot(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    PauliStr str = getPauliStr(std::string(numTargets, 'X'), targets, numTargets);
    rightapplyPauliStr(qureg, str);
}

} // end de-mangler

void leftapplyMultiQubitNot(Qureg qureg, vector<int> targets) {

    leftapplyMultiQubitNot(qureg, targets.data(), targets.size());
}

void rightapplyMultiQubitNot(Qureg qureg, vector<int> targets) {

    rightapplyMultiQubitNot(qureg, targets.data(), targets.size());
}



/*
 * projectors
 */

extern "C" {

void leftapplyQubitProjector(Qureg qureg, int qubit, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__); 

    qreal prob = 1;
    localiser_statevec_multiQubitProjector(qureg, {qubit}, {outcome}, prob);
}

void leftapplyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesAreValid(outcomes, numQubits, __func__);

    qreal prob = 1;
    auto qubitVec = util_getVector(qubits, numQubits);
    auto outcomeVec = util_getVector(outcomes, numQubits);
    localiser_statevec_multiQubitProjector(qureg, qubitVec, outcomeVec, prob);
}

void rightapplyQubitProjector(Qureg qureg, int qubit, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__); 
    
    qreal prob = 1;
    localiser_statevec_multiQubitProjector(qureg, {util_getBraQubit(qubit,qureg)}, {outcome}, prob);
}

void rightapplyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesAreValid(outcomes, numQubits, __func__);

    qreal prob = 1;
    auto qubitVec = util_getBraQubits(util_getVector(qubits, numQubits), qureg);
    auto outcomeVec = util_getVector(outcomes, numQubits);
    localiser_statevec_multiQubitProjector(qureg, qubitVec, outcomeVec, prob);
}

} // end de-mangler

void leftapplyMultiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes) {
    validate_measurementOutcomesMatchTargets(qubits.size(), outcomes.size(), __func__);

    leftapplyMultiQubitProjector(qureg, qubits.data(), outcomes.data(), outcomes.size());
}

void rightapplyMultiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes) {
    validate_measurementOutcomesMatchTargets(qubits.size(), outcomes.size(), __func__);

    rightapplyMultiQubitProjector(qureg, qubits.data(), outcomes.data(), outcomes.size());
}



/*
 * Pauli string sums
 */

extern "C" {

void leftapplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(workspace, __func__);
    validate_quregCanBeWorkspace(qureg, workspace, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    // clone qureg to workspace, set qureg to blank
    localiser_statevec_setQuregToClone(workspace, qureg);
    localiser_statevec_initUniformState(qureg, 0);

    // left-multiply each term in-turn, mixing into output qureg, then undo using idempotency
    for (qindex i=0; i<sum.numTerms; i++) {
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, sum.strings[i]);
        localiser_statevec_setQuregToWeightedSum(qureg, {1, sum.coeffs[i]}, {qureg, workspace});
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, sum.strings[i]);
    }

    // workspace -> qureg, and qureg -> sum * qureg
}

void rightapplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(workspace, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_quregCanBeWorkspace(qureg, workspace, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    // clone qureg to workspace, set qureg to blank
    localiser_statevec_setQuregToClone(workspace, qureg);
    localiser_statevec_initUniformState(qureg, 0);

    // post-multiply each term in-turn, mixing into output qureg, then undo using idempotency
    for (qindex i=0; i<sum.numTerms; i++) {
        PauliStr str =  paulis_getShiftedPauliStr(sum.strings[i], qureg.numQubits);
        qcomp factor = paulis_getSignOfPauliStrConj(str); // undoes transpose

        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, str, factor);
        localiser_statevec_setQuregToWeightedSum(qureg, {1, sum.coeffs[i]}, {qureg, workspace});
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, str, factor);
    }

    // workspace -> qureg, and qureg -> sum * qureg
}

} // end de-mangler