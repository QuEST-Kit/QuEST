/** @file
 * API definitions for effecting operators (including
 * unitaries, projectors, channels, Hermitians, and
 * arbitrary matrices) upon Quregs, which can be
 * statevectors or density matrices
 */

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/operations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions



/*
 * PRVIATE UTILITIES
 */

extern bool paulis_hasOddNumY(PauliStr str);

extern PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift);

// T can be CompMatr, CompMatr1, CompMatr2, DiagMatr, DiagMatr1, DiagMatr2
template <class T>
void validateAndApplyAnyCtrlAnyTargUnitaryMatrix(Qureg qureg, int* ctrls, int* states, int numCtrls, int* targs, int numTargs, T matr, const char* caller) {
    validate_quregFields(qureg, caller);
    validate_controlsAndTargets(qureg, ctrls, numCtrls, targs, numTargs, caller);
    validate_controlStates(states, numCtrls, caller);
    validate_matrixDimMatchesTargets(matr, numTargs, caller); // also checks fields and is-synced
    validate_matrixIsUnitary(matr, caller); // harmlessly rechecks fields and is-synced

    auto ctrlVec  = util_getVector(ctrls,  numCtrls);
    auto stateVec = util_getVector(states, numCtrls);
    auto targVec  = util_getVector(targs,  numTargs);

    bool conj = false;
    localiser_statevec_anyCtrlAnyTargAnyMatr(qureg, ctrlVec, stateVec, targVec, matr, conj);

    if (!qureg.isDensityMatrix)
        return;

    conj = true;
    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    targVec = util_getBraQubits(targVec, qureg);
    localiser_statevec_anyCtrlAnyTargAnyMatr(qureg, ctrlVec, stateVec, targVec, matr, conj);
}



/*
 * API
 */

// enable invocation by both C and C++ binaries
extern "C" {



// DEBUG
#define _NOT_IMPLEMENTED_ERROR_DEF { error_functionNotImplemented(__func__); }




/*
 * CompMatr1
 */

void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__); // matrix can be non-unitary

    bool conj = false;
    localiser_statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, target, matrix, conj);
}

void applyCompMatr1(Qureg qureg, int target, CompMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, nullptr, nullptr, 0, &target, 1, matrix, __func__);
}

void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, &control, nullptr, 1, &target, 1, matrix, __func__);
}

void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, nullptr, numControls, &target, 1, matrix, __func__);
}

void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, &target, 1, matrix, __func__);
}



/*
 * CompMatr2
 */

void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__); // matrix can be non-unitary

    bool conj = false;
    localiser_statevec_anyCtrlTwoTargDenseMatr(qureg, {}, {}, target1, target2, matrix, conj);
}

void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, nullptr, nullptr, 0, targs, 2, matrix, __func__);
}

void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, &control, nullptr, 1, targs, 2, matrix, __func__);
}

void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, nullptr, numControls, targs, 2, matrix, __func__);
}

void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, targs, 2, matrix, __func__);
}



/*
 * CompMatr
 */

void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync, but not unitarity

    bool conj = false;
    localiser_statevec_anyCtrlAnyTargDenseMatr(qureg, {}, {}, util_getVector(targets, numTargets), matrix, conj);
}

void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, nullptr, nullptr, 0, targets, numTargets, matrix, __func__);
}

void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, &control, nullptr, 1, targets, numTargets, matrix, __func__);
}

void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, nullptr, numControls, targets, numTargets, matrix, __func__);
}

void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, targets, numTargets, matrix, __func__);
}



/*
 * DiagMatr1
 */

void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixFields(matrix, __func__); // matrix can be non-unitary

    bool conj = false;
    localiser_statevec_anyCtrlOneTargDiagMatr(qureg, {}, {}, target, matrix, conj);
}

void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, nullptr, nullptr, 0, &target, 1, matrix, __func__);
}

void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, &control, nullptr, 1, &target, 1, matrix, __func__);
}

void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, nullptr, numControls, &target, 1, matrix, __func__);
}

void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, &target, 1, matrix, __func__);
}



/*
 * DiagMatr2
 */

void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);
    validate_matrixFields(matrix, __func__); // matrix can be non-unitary

    bool conj = false;
    localiser_statevec_anyCtrlTwoTargDiagMatr(qureg, {}, {}, target1, target2, matrix, conj);
}

void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, nullptr, nullptr, 0, targs, 2, matrix, __func__);
}

void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, &control, nullptr, 1, targs, 2, matrix, __func__);
}

void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, nullptr, numControls, targs, 2, matrix, __func__);
}

void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matrix) {

    int targs[] = {target1, target2};
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, targs, 2, matrix, __func__);
}



/*
 * DiagMatr
 */

void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync, but not unitarity

    bool conj = false;
    qcomp exponent = 1;
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, util_getVector(targets, numTargets), matrix, exponent, conj);
}

void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, nullptr, nullptr, 0, targets, numTargets, matrix, __func__);
}

void applyControlledDiagMatr(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, &control, nullptr, 1, targets, numTargets, matrix, __func__);
}

void applyMultiControlledDiagMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, nullptr, numControls, targets, numTargets, matrix, __func__);
}

void applyMultiStateControlledDiagMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix) {

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, targets, numTargets, matrix, __func__);
}



/*
 * DiagMatrPower
 *
 * which still (except for multiply) assert unitarity,
 * even though a non-real exponent is possible
 */

void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also validates fields and is-sync, but not unitarity

    bool conj = false;
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, {}, {}, util_getVector(targets, numTargets), matrix, exponent, conj);
}

void applyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent)  {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced

    // harmlessly re-validates
    applyMultiStateControlledDiagMatrPower(qureg, nullptr, nullptr, 0, targets, numTargets, matrix, exponent);
}

void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTargets(qureg, control, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced

    // harmlessly re-validates
    applyMultiStateControlledDiagMatrPower(qureg, &control, nullptr, 1, targets, numTargets, matrix, exponent);
}

void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced

    // harmlessly re-validates
    applyMultiStateControlledDiagMatrPower(qureg, controls, nullptr, numControls, targets, numTargets, matrix, exponent);
}

void applyMultiStateControlledDiagMatrPower(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);
    validate_controlStates(states, numControls, __func__); // can be nullptr, ignoring numControls
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced

    bool conj = false;
    auto ctrlVec = util_getVector(controls, numControls);
    auto stateVec = util_getVector(states,  numControls); // empty if states==nullptr
    auto targVec = util_getVector(targets,  numTargets);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, ctrlVec, stateVec, targVec, matrix, exponent, conj);

    if (!qureg.isDensityMatrix)
        return;

    conj = true;
    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    targVec = util_getBraQubits(targVec, qureg);
    localiser_statevec_anyCtrlAnyTargAnyMatr(qureg, {}, {}, targVec, matrix, conj);
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, ctrlVec, stateVec, targVec, matrix, exponent, conj);
}



/*
 * FullStateDiagMatr
 */

void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, __func__); // matrix can be non-unitary

    bool onlyMultiply = true;
    qcomp exponent = qcomp(1, 0);
    (qureg.isDensityMatrix)?
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent) :
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, onlyMultiply);
}

void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, __func__); // matrix can be non-unitary

    bool onlyMultiply = true;
    (qureg.isDensityMatrix)?
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent) :
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, onlyMultiply);
}

void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, __func__);
    validate_matrixIsUnitary(matrix, __func__);

    bool onlyMultiply = false;
    qcomp exponent = qcomp(1, 0);
    (qureg.isDensityMatrix)?
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent) :
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, onlyMultiply);
}

void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, __func__);
    validate_matrixIsUnitary(matrix, __func__);

    bool onlyMultiply = false;
    (qureg.isDensityMatrix)?
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent) :
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, onlyMultiply);
}



/*
 * swaps
 */

void multiplySwap(Qureg qureg, int qubit1, int qubit2)
    _NOT_IMPLEMENTED_ERROR_DEF

void applySwap(Qureg qureg, int qubit1, int qubit2)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * individual Paulis
 *
 * where Y and Z are most efficiently effected as DiagMatr1,
 * but where X is best effected as a 1-qubit PauliStr.
 */

void applyPauliX(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliX(qureg, nullptr, nullptr, 0, target);
}

void applyPauliY(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliY(qureg, nullptr, nullptr, 0, target);
}

void applyPauliZ(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliZ(qureg, nullptr, nullptr, 0, target);
}

void applyControlledPauliX(Qureg qureg, int control, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliX(qureg, &control, nullptr, 1, target);
}

void applyControlledPauliY(Qureg qureg, int control, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliY(qureg, &control, nullptr, 1, target);
}

void applyControlledPauliZ(Qureg qureg, int control, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliZ(qureg, &control, nullptr, 1, target);
}

void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliX(qureg, controls, nullptr, numControls, target);
}

void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliY(qureg, controls, nullptr, numControls, target);
}

void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliZ(qureg, controls, nullptr, numControls, target);
}

void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target)  {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // harmlessly re-validates
    applyMultiStateControlledPauliStr(qureg, controls, states, numControls, getPauliStr("X", {target}));
}

void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // harmlessly re-validates
    DiagMatr1 matr = getDiagMatr1({-1_i, 1_i});
    applyMultiStateControlledDiagMatr1(qureg, controls, states, numControls, target, matr);
}

void applyMultiStateControlledPauliZ(Qureg qureg, int* controls, int* states, int numControls, int target)  {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // harmlessly re-validates
    DiagMatr1 matr = getDiagMatr1({1, -1});
    applyMultiStateControlledDiagMatr1(qureg, controls, states, numControls, target, matr);
}



/*
 * Pauli strings
 */

void multiplyPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    localiser_statevec_anyCtrlPauliTensor(qureg, {}, {}, str);
}

void applyPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliStr(qureg, nullptr, nullptr, 0, str);
}

void applyControlledPauliStr(Qureg qureg, int control, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_controlAndPauliStrTargets(qureg, control, str, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliStr(qureg, &control, nullptr, 1, str);
}

void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndPauliStrTargets(qureg, controls, numControls, str, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPauliStr(qureg, controls, nullptr, numControls, str);
}

void applyMultiStateControlledPauliStr(Qureg qureg, int* controls, int* states, int numControls, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndPauliStrTargets(qureg, controls, numControls, str, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    auto ctrlVec = util_getVector(controls, numControls);
    auto stateVec = util_getVector(states, numControls); // empty if states==nullptr
    localiser_statevec_anyCtrlPauliTensor(qureg, ctrlVec, stateVec, str);

    if (!qureg.isDensityMatrix)
        return;

    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliTensor(qureg, ctrlVec, stateVec, str); // excludes conj

    // effect conj by qureg *= -1
    if (paulis_hasOddNumY(str))
        localiser_statevec_setWeightedQureg(-1, qureg, 0, qureg, 0, qureg);
}



/*
 * Pauli string sums
 */

void multiplyPauliStrSum(Qureg qureg, PauliStrSum str)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyTrotterizedTimeEvol(Qureg qureg, PauliStrSum hamiltonian, qreal time, int order, int reps) {

    // validate that PauliStrSum is Hermitian
    
    // TODO
    error_functionNotImplemented(__func__);
}



/*
 * individual axis rotations
 */

// don't think users will ever want to left-multiply only

void applyRotateX(Qureg qureg, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateX(qureg, nullptr, nullptr, 0, target, angle);
}

void applyRotateY(Qureg qureg, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateY(qureg, nullptr, nullptr, 0, target, angle);
}

void applyRotateZ(Qureg qureg, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateZ(qureg, nullptr, nullptr, 0, target, angle);
}

void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateX(qureg, &control, nullptr, 1, target, angle);
}

void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateY(qureg, &control, nullptr, 1, target, angle);
}

void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateZ(qureg, &control, nullptr, 1, target, angle);
}

void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateX(qureg, controls, nullptr, numControls, target, angle);
}

void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateY(qureg, controls, nullptr, numControls, target, angle);
}

void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledRotateZ(qureg, controls, nullptr, numControls, target, angle);
}

void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // harmlessly re-validates
    applyMultiStateControlledPauliGadget(qureg, controls, states, numControls, getPauliStr("X", {target}), angle);
}

void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // harmlessly re-validates
    applyMultiStateControlledPauliGadget(qureg, controls, states, numControls, getPauliStr("Y", {target}), angle);
}

void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // harmlessly re-validates
    applyMultiStateControlledPauliGadget(qureg, controls, states, numControls, getPauliStr("Z", {target}), angle);
}



/*
 * arbitrary axis rotation
 */

void applyRotateAroundAxis(Qureg qureg, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledRotateAroundAxis(Qureg qureg, int control, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * Pauli gadgets
 */

void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    localiser_statevec_anyCtrlPauliGadget(qureg, {}, {}, str, angle);
}

void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);
    
    applyMultiStateControlledPauliGadget(qureg, nullptr, nullptr, 0, str, angle);
}

void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlAndPauliStrTargets(qureg, control, str, __func__);
    
    applyMultiStateControlledPauliGadget(qureg, &control, nullptr, 1, str, angle);
}

void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndPauliStrTargets(qureg, controls, numControls, str, __func__);
    
    applyMultiStateControlledPauliGadget(qureg, controls, nullptr, numControls, str, angle);
}

void applyMultiStateControlledPauliGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndPauliStrTargets(qureg, controls, numControls, str, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    // TODO:
    // CRUCIAL NOTE:
    // exp(theta I..I) might be algorithmically ok (I'm not sure), but it WILL NOT
    // effect a global phase change of theta (I think). Should validate against this
    // sitaution just in case, or make the doc extremely explicit

    auto ctrlVec = util_getVector(controls, numControls);
    auto stateVec = util_getVector(states, numControls); // empty if states==nullptr
    localiser_statevec_anyCtrlPauliGadget(qureg, ctrlVec, stateVec, str, angle);

    if (!qureg.isDensityMatrix)
        return;

    angle *= paulis_hasOddNumY(str) ? -1 : +1;
    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliGadget(qureg, ctrlVec, stateVec, str, angle);
}



/*
 * phase gadgets
 */

void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    localiser_statevec_anyCtrlPhaseGadget(qureg, {}, {}, util_getVector(targets,numTargets), angle);
}

void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPhaseGadget(qureg, nullptr, nullptr, 0, targets, numTargets, angle);
}

void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTargets(qureg, control, targets, numTargets, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPhaseGadget(qureg, &control, nullptr, 1, targets, numTargets, angle);
}

void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);

    // harmlessly re-validates
    applyMultiStateControlledPhaseGadget(qureg, controls, nullptr, numControls, targets, numTargets, angle);
}

void applyMultiStateControlledPhaseGadget(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);
    validate_controlStates(states, numControls, __func__);

    auto ctrlVec = util_getVector(controls, numControls);
    auto targVec = util_getVector(targets,  numTargets);
    auto stateVec = util_getVector(states,  numControls); // empty if states==nullptr
    localiser_statevec_anyCtrlPhaseGadget(qureg, ctrlVec, stateVec, targVec, angle);

    if (!qureg.isDensityMatrix)
        return;

    angle *= -1;
    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    targVec = util_getBraQubits(ctrlVec, qureg);
    localiser_statevec_anyCtrlPhaseGadget(qureg, ctrlVec, stateVec, targVec, angle);
}



/*
 * phase shift
 */

void applyPhaseShift(Qureg qureg, int target, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiQubitPhaseShift(qureg, &target, 1, angle);
}

void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);

    // harmlessly re-validates
    int targets[] = {target1, target2};
    applyMultiQubitPhaseShift(qureg, targets, 2, angle);
}

void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // treat as a (numTargets-1)-controlled 1-target diagonal matrix
    DiagMatr1 matr = getDiagMatr1({1, exp(1_i * angle)});

    // harmlessly re-validates
    applyMultiStateControlledDiagMatr1(qureg, &targets[1], nullptr, numTargets-1, targets[0], matr);
}



/*
 * phase flips
 */

void applyPhaseFlip(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiQubitPhaseFlip(qureg, &target, 1);
}

void applyTwoQubitPhaseFlip(Qureg qureg, int target1, int target2) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, target1, target2, __func__);

    // harmlessly re-validates
    int targets[] = {target1, target2};
    applyMultiQubitPhaseFlip(qureg, targets, 2);
}

void applyMultiQubitPhaseFlip(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // treat as a (numTargets-1)-controlled 1-target Pauli Z
    DiagMatr1 matr = getDiagMatr1({1, -1});

    // harmlessly re-validates
    applyMultiStateControlledDiagMatr1(qureg, &targets[1], nullptr, numTargets-1, targets[0], matr);
}



/*
 * many-qubit CNOTs
 */

void applyNot(Qureg qureg, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledNot(Qureg qureg, int control, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledNot(Qureg qureg, int* controls, int numControls, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledNot(Qureg qureg, int* controls, int* states, int numControls, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiQubitNot(Qureg, int* targets, int numTargets)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledMultiQubitNot(Qureg, int control, int* targets, int numTargets)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledMultiQubitNot(Qureg, int* controls, int numControls, int* targets, int numTargets)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledMultiQubitNot(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * superoperator
 */

void applySuperOp(Qureg qureg, SuperOp superop, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);
    validate_superOpFields(superop, __func__);
    validate_superOpIsSynced(superop, __func__);
    validate_superOpDimMatchesTargs(superop, numTargets, __func__);

    localiser_densmatr_superoperator(qureg, superop, util_getVector(targets, numTargets));
}



/*
 * measurement
 */

int applyQubitMeasurement(Qureg qureg, int target) {

    // // TODO
    error_functionNotImplemented(__func__);
    return -1;
}

int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability) {

    // // TODO
    error_functionNotImplemented(__func__);
    return -1;
}

qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome) {

    // // TODO
    error_functionNotImplemented(__func__);
    return -1;
}

void applyQubitProjector(Qureg qureg, int target, int outcome)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * QFT
 */

void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyFullQuantumFourierTransform(Qureg qureg)
    _NOT_IMPLEMENTED_ERROR_DEF



} // end de-mangler