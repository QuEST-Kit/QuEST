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
 */

void multiplyPauliX(Qureg qureg, int target) _NOT_IMPLEMENTED_ERROR_DEF
void multiplyPauliY(Qureg qureg, int target) _NOT_IMPLEMENTED_ERROR_DEF
void multiplyPauliZ(Qureg qureg, int target) _NOT_IMPLEMENTED_ERROR_DEF

void applyPauliX(Qureg qureg, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyPauliY(Qureg qureg, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyPauliZ(Qureg qureg, int target) _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledPauliX(Qureg qureg, int control, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyControlledPauliY(Qureg qureg, int control, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyControlledPauliZ(Qureg qureg, int control, int target) _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target) _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiStateControlledPauliZ(Qureg qureg, int* controls, int* states, int numControls, int target) _NOT_IMPLEMENTED_ERROR_DEF



/*
 * Pauli strings
 */

void multiplyPauliStr(Qureg qureg, PauliStr str)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyPauliStr(Qureg qureg, PauliStr str)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledPauliStr(Qureg qureg, int control, PauliStr str)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledPauliStr(Qureg qureg, int* controls, int* states, int numControls, PauliStr str)
    _NOT_IMPLEMENTED_ERROR_DEF



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

void applyRotateX(Qureg qureg, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyRotateY(Qureg qureg, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyRotateZ(Qureg qureg, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF
void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle) _NOT_IMPLEMENTED_ERROR_DEF



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

// don't think users will ever want to left-multiply only

void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledPauliGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStr str, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * phase gadgets
 */

void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledPhaseGadget(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * phase shift
 */

void applyPhaseShift(Qureg qureg, int target, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * phase flips
 */

void applyPhaseFlip(Qureg qureg, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyTwoQubitPhaseFlip(Qureg qureg, int target1, int target2)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiQubitPhaseFlip(Qureg qureg, int* targets, int numTargets)
    _NOT_IMPLEMENTED_ERROR_DEF



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