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

#include <vector>

using std::vector;



/*
 * CompMatr1
 */

extern "C" {

void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    bool transp = false;
    localiser_statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, target, matrix, conj, transp);
}

void postMultiplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix) {
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

void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__);
    validate_mixedAmpsFitInNode(qureg, 2, __func__);

    bool conj = false;
    bool transp = false;
    localiser_statevec_anyCtrlTwoTargDenseMatr(qureg, {}, {}, target1, target2, matrix, conj, transp);
}

void postMultiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix) {
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

void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync
    validate_mixedAmpsFitInNode(qureg, numTargets, __func__);

    bool conj = false;
    bool transp = false;
    localiser_statevec_anyCtrlAnyTargDenseMatr(qureg, {}, {}, util_getVector(targets, numTargets), matrix, conj, transp);
}

void postMultiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix) {
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

void multiplyCompMatr(Qureg qureg, vector<int> targets, CompMatr matr) {

    multiplyCompMatr(qureg, targets.data(), targets.size(), matr);
}

void postMultiplyCompMatr(Qureg qureg, vector<int> targets, CompMatr matr) {

    postMultiplyCompMatr(qureg, targets.data(), targets.size(), matr);
}



/*
 * DiagMatr1
 */

extern "C" {

void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    localiser_statevec_anyCtrlOneTargDiagMatr(qureg, {}, {}, target, matrix, conj);
}

void postMultiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix) {
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

void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__);

    bool conj = false;
    localiser_statevec_anyCtrlTwoTargDiagMatr(qureg, {}, {}, target1, target2, matrix, conj);
}

void postMultiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix) {
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

void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync

    bool conj = false;
    qcomp exponent = 1;
    auto qubits = util_getVector(targets, numTargets);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, qubits, matrix, exponent, conj);
}

void postMultiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix) {
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

void multiplyDiagMatr(Qureg qureg, vector<int> targets, DiagMatr matrix) {

    multiplyDiagMatr(qureg, targets.data(), targets.size(), matrix);
}

void postMultiplyDiagMatr(Qureg qureg, vector<int> targets, DiagMatr matrix) {

    postMultiplyDiagMatr(qureg, targets.data(), targets.size(), matrix);
}



/*
 * DiagMatrPower
 */

extern "C" {

void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync, but not unitarity
    validate_matrixExpIsNonDiverging(matrix, exponent, __func__); // harmlessly re-validates fields and is-sync

    bool conj = false;
    auto qubits = util_getVector(targets, numTargets);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, qubits, matrix, exponent, conj);
}

void postMultiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
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

void multiplyDiagMatrPower(Qureg qureg, vector<int> targets, DiagMatr matrix, qcomp exponent) {

    multiplyDiagMatrPower(qureg, targets.data(), targets.size(), matrix, exponent);
}

void postMultiplyDiagMatrPower(Qureg qureg, vector<int> targets, DiagMatr matrix, qcomp exponent) {

    postMultiplyDiagMatrPower(qureg, targets.data(), targets.size(), matrix, exponent);
}



/*
 * FullStateDiagMatr (and power)
 */

extern "C" {

void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__); // matrix can be non-unitary

    multiplyFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}

void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
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

void postMultiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__); // matrix can be non-unitary

    postMultiplyFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}

void postMultiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
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

void multiplySwap(Qureg qureg, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);

    localiser_statevec_anyCtrlSwap(qureg, {}, {}, qubit1, qubit2);
}

void postMultiplySwap(Qureg qureg, int qubit1, int qubit2) {
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

extern PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift);

extern "C" {

void multiplyPauliX(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("X", {target});
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void multiplyPauliY(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("Y", {target});
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void multiplyPauliZ(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("Z", {target});
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void postMultiplyPauliX(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);

    PauliStr str = getPauliStr("X", {target});
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void postMultiplyPauliY(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, target, __func__);

    qcomp factor = -1; // undo transpose
    PauliStr str = getPauliStr("Y", {target});
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str, factor);
}

void postMultiplyPauliZ(Qureg qureg, int target) {
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

extern bool paulis_hasOddNumY(PauliStr str);

extern "C" {

void multiplyPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void postMultiplyPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qcomp factor = paulis_hasOddNumY(str)? -1 : 1; // undo transpose
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str, factor);
}

} // end de-mangler



/*
 * Pauli gadgets
 */

extern "C" {

void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    localiser_statevec_anyCtrlPauliGadget(qureg, {}, {}, str, phase);
}

void postMultiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qreal factor = paulis_hasOddNumY(str)? -1 : 1;
    qreal phase = factor * util_getPhaseFromGateAngle(angle);
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliGadget(qureg, {}, {}, str, phase);
}

} // end de-mangler



/*
 * phase gadgets
 */

extern "C" {

void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    auto qubits = util_getVector(targets, numTargets);
    localiser_statevec_anyCtrlPhaseGadget(qureg, {}, {}, qubits, phase);
}

void postMultiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    auto qubits = util_getBraQubits(util_getVector(targets, numTargets), qureg);
    localiser_statevec_anyCtrlPhaseGadget(qureg, {}, {}, qubits, phase);
}

} // end de-mangler

void multiplyPhaseGadget(Qureg qureg, vector<int> targets, qreal angle) {

    multiplyPhaseGadget(qureg, targets.data(), targets.size(), angle);
}

void postMultiplyPhaseGadget(Qureg qureg, vector<int> targets, qreal angle) {

    postMultiplyPhaseGadget(qureg, targets.data(), targets.size(), angle);
}



/*
 * many-qubit NOTs
 */

extern "C" {

void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    PauliStr str = getPauliStr(std::string(numTargets, 'X'), targets, numTargets);
    multiplyPauliStr(qureg, str);
}

void postMultiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    PauliStr str = getPauliStr(std::string(numTargets, 'X'), targets, numTargets);
    postMultiplyPauliStr(qureg, str);
}

} // end de-mangler

void multiplyMultiQubitNot(Qureg qureg, vector<int> targets) {

    multiplyMultiQubitNot(qureg, targets.data(), targets.size());
}

void postMultiplyMultiQubitNot(Qureg qureg, vector<int> targets) {

    postMultiplyMultiQubitNot(qureg, targets.data(), targets.size());
}



/*
 * projectors
 */

extern "C" {

void multiplyQubitProjector(Qureg qureg, int qubit, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__); 

    qreal prob = 1;
    localiser_statevec_multiQubitProjector(qureg, {qubit}, {outcome}, prob);
}

void multiplyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesAreValid(outcomes, numQubits, __func__);

    qreal prob = 1;
    auto qubitVec = util_getVector(qubits, numQubits);
    auto outcomeVec = util_getVector(outcomes, numQubits);
    localiser_statevec_multiQubitProjector(qureg, qubitVec, outcomeVec, prob);
}

void postMultiplyQubitProjector(Qureg qureg, int qubit, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__); 
    
    qreal prob = 1;
    localiser_statevec_multiQubitProjector(qureg, {util_getBraQubit(qubit,qureg)}, {outcome}, prob);
}

void postMultiplyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
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

void multiplyMultiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes) {
    validate_measurementOutcomesMatchTargets(qubits.size(), outcomes.size(), __func__);

    multiplyMultiQubitProjector(qureg, qubits.data(), outcomes.data(), outcomes.size());
}

void postMultiplyMultiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes) {
    validate_measurementOutcomesMatchTargets(qubits.size(), outcomes.size(), __func__);

    postMultiplyMultiQubitProjector(qureg, qubits.data(), outcomes.data(), outcomes.size());
}



/*
 * Pauli string sums
 */

extern "C" {

void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(workspace, __func__);
    validate_quregCanBeWorkspace(qureg, workspace, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    // clone qureg to workspace, set qureg to blank
    localiser_statevec_setQuregToSuperposition(0, workspace, 1, qureg, 0, qureg);
    localiser_statevec_initUniformState(qureg, 0);

    // left-multiply each term in-turn, mixing into output qureg, then undo using idempotency
    for (qindex i=0; i<sum.numTerms; i++) {
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, sum.strings[i]);
        localiser_statevec_setQuregToSuperposition(1, qureg, sum.coeffs[i], workspace, 0, workspace);
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, sum.strings[i]);
    }

    // workspace -> qureg, and qureg -> sum * qureg
}

void postMultiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(workspace, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_quregCanBeWorkspace(qureg, workspace, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    // clone qureg to workspace, set qureg to blank
    localiser_statevec_setQuregToSuperposition(0, workspace, 1, qureg, 0, qureg);
    localiser_statevec_initUniformState(qureg, 0);

    // post-multiply each term in-turn, mixing into output qureg, then undo using idempotency
    for (qindex i=0; i<sum.numTerms; i++) {
        PauliStr str =  paulis_getShiftedPauliStr(sum.strings[i], qureg.numQubits);
        qcomp factor = paulis_hasOddNumY(str)? -1 : 1; // undoes transpose

        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, str, factor);
        localiser_statevec_setQuregToSuperposition(1, qureg, sum.coeffs[i], workspace, 0, workspace);
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, str, factor);
    }

    // workspace -> qureg, and qureg -> sum * qureg
}

} // end de-mangler