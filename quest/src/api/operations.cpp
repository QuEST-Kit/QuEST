/** @file
 * API definitions for effecting operators (including
 * unitaries, projectors, channels, Hermitians, and
 * arbitrary matrices) upon Quregs, which can be
 * statevectors or density matrices.
 * 
 * @author Tyson Jones
 */

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/operations.h"
#include "quest/include/calculations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/randomiser.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/constants.hpp"

#include <vector>

using std::vector;



/*
 * PRVIATE UTILITIES
 */

extern bool paulis_hasOddNumY(PauliStr str);
extern PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift);
extern PauliStr paulis_getKetAndBraPauliStr(PauliStr str, Qureg qureg);

// T can be CompMatr, CompMatr1, CompMatr2, DiagMatr, DiagMatr1, DiagMatr2
template <class T>
void validateAndApplyAnyCtrlAnyTargUnitaryMatrix(Qureg qureg, int* ctrls, int* states, int numCtrls, int* targs, int numTargs, T matr, const char* caller) {
    validate_quregFields(qureg, caller);
    validate_controlsAndTargets(qureg, ctrls, numCtrls, targs, numTargs, caller);
    validate_controlStates(states, numCtrls, caller);
    validate_matrixDimMatchesTargets(matr, numTargs, caller); // also checks fields and is-synced
    validate_matrixIsUnitary(matr, caller); // harmlessly rechecks fields and is-synced
    if (util_isDenseMatrixType<T>())
        validate_mixedAmpsFitInNode(qureg, numTargs, caller);

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
    validate_mixedAmpsFitInNode(qureg, 2, __func__);

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
    validate_mixedAmpsFitInNode(qureg, numTargets, __func__);

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
    validate_exponentIsReal(exponent, __func__); // checks matrix^exponent is unitary (abs=1)

    // harmlessly re-validates
    applyMultiStateControlledDiagMatrPower(qureg, nullptr, nullptr, 0, targets, numTargets, matrix, exponent);
}

void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTargets(qureg, control, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced
    validate_exponentIsReal(exponent, __func__); // checks matrix^exponent is unitary (abs=1)

    // harmlessly re-validates
    applyMultiStateControlledDiagMatrPower(qureg, &control, nullptr, 1, targets, numTargets, matrix, exponent);
}

void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced
    validate_exponentIsReal(exponent, __func__); // checks matrix^exponent is unitary (abs=1)

    // harmlessly re-validates
    applyMultiStateControlledDiagMatrPower(qureg, controls, nullptr, numControls, targets, numTargets, matrix, exponent);
}

void applyMultiStateControlledDiagMatrPower(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);
    validate_controlStates(states, numControls, __func__); // can be nullptr, ignoring numControls
    validate_matrixDimMatchesTargets(matrix, numTargets, __func__); // also checks fields and is-synced
    validate_matrixIsUnitary(matrix, __func__); // harmlessly rechecks fields and is-synced
    validate_exponentIsReal(exponent, __func__); // checks matrix^exponent is unitary (abs=1)

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
    localiser_statevec_anyCtrlAnyTargDiagMatr(qureg, ctrlVec, stateVec, targVec, matrix, exponent, conj);
}



/*
 * FullStateDiagMatr (and power)
 */

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

    bool onlyMultiply = true;
    (qureg.isDensityMatrix)?
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, onlyMultiply):
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent);
}

void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__);
    validate_matrixIsUnitary(matrix, __func__);

    applyFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}

void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, false, __func__);
    validate_matrixIsUnitary(matrix, __func__);  // checks abs=1
    validate_exponentIsReal(exponent, __func__); // checks matrix^exponent is unitary (abs=1)

    bool onlyMultiply = false;
    (qureg.isDensityMatrix)?
        localiser_densmatr_allTargDiagMatr(qureg, matrix, exponent, onlyMultiply):
        localiser_statevec_allTargDiagMatr(qureg, matrix, exponent);
}



/*
 * S gate
 */

void applyS(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledS(qureg, nullptr, nullptr, 0, target);
}

void applyControlledS(Qureg qureg, int control, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledS(qureg, &control, nullptr, 1, target);
}

void applyMultiControlledS(Qureg qureg, int* controls, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledS(qureg, controls, nullptr, numControls, target);
}

void applyMultiStateControlledS(Qureg qureg, int* controls, int* states, int numControls, int target) {

    DiagMatr1 matr = getDiagMatr1({1, 1_i});
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, &target, 1, matr, __func__);
}



/*
 * T gate
 */

void applyT(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledT(qureg, nullptr, nullptr, 0, target);
}

void applyControlledT(Qureg qureg, int control, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledT(qureg, &control, nullptr, 1, target);
}

void applyMultiControlledT(Qureg qureg, int* controls, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledT(qureg, controls, nullptr, numControls, target);
}

void applyMultiStateControlledT(Qureg qureg, int* controls, int* states, int numControls, int target) {

    DiagMatr1 matr = getDiagMatr1({1, 1/std::sqrt(2) + 1_i/std::sqrt(2)});
    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, &target, 1, matr, __func__);
}



/*
 * Hadamard 
 */

void applyHadamard(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledHadamard(qureg, nullptr, nullptr, 0, target);
}

void applyControlledHadamard(Qureg qureg, int control, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, control, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledHadamard(qureg, &control, nullptr, 1, target);
}

void applyMultiControlledHadamard(Qureg qureg, int* controls, int numControls, int target) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, controls, numControls, target, __func__);

    // harmlessly re-validates
    applyMultiStateControlledHadamard(qureg, controls, nullptr, numControls, target);
}

void applyMultiStateControlledHadamard(Qureg qureg, int* controls, int* states, int numControls, int target) {

    qcomp a = 1/std::sqrt(2);
    CompMatr1 matr = getCompMatr1({
        {a, a}, 
        {a,-a}});

    validateAndApplyAnyCtrlAnyTargUnitaryMatrix(qureg, controls, states, numControls, &target, 1, matr, __func__);
}



/*
 * swap
 */

void multiplySwap(Qureg qureg, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);

    localiser_statevec_anyCtrlSwap(qureg, {}, {}, qubit1, qubit2);
}

void applySwap(Qureg qureg, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg, qubit1, qubit2, __func__);

    // harmlessly re-valdiates
    applyMultiStateControlledSwap(qureg, nullptr, nullptr, 0, qubit1, qubit2);
}

void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTwoTargets(qureg, control, qubit1, qubit2, __func__);

    // harmlessly re-valdiates
    applyMultiStateControlledSwap(qureg, &control, nullptr, 1, qubit1, qubit2);
}

void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTwoTargets(qureg, controls, numControls, qubit1, qubit2, __func__);

    // harmlessly re-valdiates
    applyMultiStateControlledSwap(qureg, controls, nullptr, numControls, qubit1, qubit2);
}

void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTwoTargets(qureg, controls, numControls, qubit1, qubit2, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    auto ctrlVec = util_getVector(controls, numControls);
    auto stateVec = util_getVector(states, numControls); // empty if states==nullptr
    localiser_statevec_anyCtrlSwap(qureg, ctrlVec, stateVec, qubit1, qubit2);

    if (!qureg.isDensityMatrix)
        return;

    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    qubit1 = util_getBraQubit(qubit1, qureg);
    qubit2 = util_getBraQubit(qubit2, qureg);
    localiser_statevec_anyCtrlSwap(qureg, ctrlVec, stateVec, qubit1, qubit2);
}



/*
 * sqrt swap
 */

void applySqrtSwap(Qureg qureg, int target1, int target2) {
    validate_quregFields(qureg, __func__);
    validate_twoTargets(qureg,target1, target2, __func__);

    // harmlessly re-validates
    applyMultiStateControlledSqrtSwap(qureg, nullptr, nullptr, 0, target1, target2);
}

void applyControlledSqrtSwap(Qureg qureg, int control, int target1, int target2) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTwoTargets(qureg, control, target1, target2, __func__);

    // harmlessly re-validates
    applyMultiStateControlledSqrtSwap(qureg, &control, nullptr, 1, target1, target2);
}

void applyMultiControlledSqrtSwap(Qureg qureg, int* controls, int numControls, int target1, int target2) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTwoTargets(qureg, controls, numControls, target1, target2, __func__);

    // harmlessly re-validates
    applyMultiStateControlledSqrtSwap(qureg, controls, nullptr, numControls, target1, target2);
}

void applyMultiStateControlledSqrtSwap(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTwoTargets(qureg, controls, numControls, target1, target2, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr

    /// @todo
    /// this function exacts sqrtSwap as a dense 2-qubit matrix,
    /// where as bespoke communication and simulation strategy is
    /// clearly possible which we have not supported because the gate
    /// is somewhat esoteric. As such, we must validate mixed-amps fit

    validate_mixedAmpsFitInNode(qureg, 2, __func__); // to throw SqrtSwap error, not generic CompMatr2 error

    CompMatr2 matr = getCompMatr2({
        {1, 0, 0, 0},
        {0, .5+.5_i, .5-.5_i, 0},
        {0, .5-.5_i, .5+.5_i, 0},
        {0, 0, 0, 1}});

    applyMultiStateControlledCompMatr2(qureg, controls, states, numControls, target1, target2, matr);
}



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
    applyMultiStateControlledPauliStr(qureg, controls, states, numControls, getPauliStr("Y", {target}));
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

    qcomp factor = 1;
    auto ctrlVec = util_getVector(controls, numControls);
    auto stateVec = util_getVector(states, numControls); // empty if states==nullptr

    // when there are no control qubits, we can merge the density matrix's 
    // operation sinto a single tensor, i.e. +- (shift(str) (x) str), to 
    // avoid superfluous re-enumeration of the state
    if (qureg.isDensityMatrix && numControls == 0) {
        factor = paulis_hasOddNumY(str)? -1 : 1;
        ctrlVec = util_getConcatenated(ctrlVec, util_getBraQubits(ctrlVec, qureg));
        stateVec = util_getConcatenated(stateVec, stateVec); 
        str = paulis_getKetAndBraPauliStr(str, qureg);
    }

    localiser_statevec_anyCtrlPauliTensor(qureg, ctrlVec, stateVec, str, factor);

    // but density-matrix control qubits require two distinct operations
    if (qureg.isDensityMatrix && numControls > 0) {
        factor = paulis_hasOddNumY(str)? -1 : 1;
        ctrlVec = util_getBraQubits(ctrlVec, qureg);
        str = paulis_getShiftedPauliStr(str, qureg.numQubits);
        localiser_statevec_anyCtrlPauliTensor(qureg, ctrlVec, stateVec, str, factor);
    }
}



/*
 * Pauli string sums
 */

void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(workspace, __func__);
    validate_quregCanBeWorkspace(qureg, workspace, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    // clone qureg to workspace, set qureg to blank
    localiser_statevec_setQuregToSuperposition(0, workspace, 1, qureg, 0, qureg);
    localiser_statevec_initUniformState(qureg, 0);

    // apply each term in-turn, mixing into output qureg, then undo using idempotency
    for (qindex i=0; i<sum.numTerms; i++) {
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, sum.strings[i]);
        localiser_statevec_setQuregToSuperposition(1, qureg, sum.coeffs[i], workspace, 0, workspace);
        localiser_statevec_anyCtrlPauliTensor(workspace, {}, {}, sum.strings[i]);
    }

    // workspace -> qureg, and qureg -> sum * qureg
}

void applyFirstOrderTrotter(Qureg qureg, PauliStrSum sum, qreal angle, bool reverse) {

    // (internal, invoked by applyTrotterizedPauliStrSumGadget)

    for (qindex i=0; i<sum.numTerms; i++) {
        int j = reverse? sum.numTerms - i - 1 : i;
        qreal arg = 2 * angle * std::real(sum.coeffs[j]);  // 2 undoes Gadget convention
        applyPauliGadget(qureg, sum.strings[j], arg); // re-validates, grr
    }
}

void applyHigherOrderTrotter(Qureg qureg, PauliStrSum sum, qreal angle, int order) {

    // (internal, invoked by applyTrotterizedPauliStrSumGadget)

    if (order == 1) {
        applyFirstOrderTrotter(qureg, sum, angle, false);
    
    } else if (order == 2) {
        applyFirstOrderTrotter(qureg, sum, angle/2, false);
        applyFirstOrderTrotter(qureg, sum, angle/2, true);
    
    } else {
        qreal p = 1. / (4 - std::pow(4, 1./(order-1)));
        qreal a = p * angle;
        qreal b = (1-4*p) * angle;

        int lower = order - 2;
        applyFirstOrderTrotter(qureg, sum, a, lower);
        applyFirstOrderTrotter(qureg, sum, a, lower);
        applyFirstOrderTrotter(qureg, sum, b, lower);
        applyFirstOrderTrotter(qureg, sum, a, lower);
        applyFirstOrderTrotter(qureg, sum, a, lower);
    }
}

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    /// @todo
    /// the accuracy of Trotterisation is greatly improved by randomisation
    /// or (even sub-optimal) grouping into commuting terms. Should we 
    /// implement these here or into another function?

    if (angle == 0)
        return;

    for (int r=0; r<reps; r++)
        applyHigherOrderTrotter(qureg, sum, angle/reps, order);
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

void applyRotateAroundAxis(Qureg qureg, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, targ, __func__);
    validate_rotationAxisNotZeroVector(axisX, axisY, axisZ, __func__);

    applyMultiStateControlledRotateAroundAxis(qureg, nullptr, nullptr, 0, targ, angle, axisX, axisY, axisZ);
}

void applyControlledRotateAroundAxis(Qureg qureg, int ctrl, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTarget(qureg, ctrl, targ, __func__);
    validate_rotationAxisNotZeroVector(axisX, axisY, axisZ, __func__);

    applyMultiStateControlledRotateAroundAxis(qureg, &ctrl, nullptr, 1, targ, angle, axisX, axisY, axisZ);
}

void applyMultiControlledRotateAroundAxis(Qureg qureg, int* ctrls, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, ctrls, numCtrls, targ, __func__);
    validate_rotationAxisNotZeroVector(axisX, axisY, axisZ, __func__);

    applyMultiStateControlledRotateAroundAxis(qureg, ctrls, nullptr, numCtrls, targ, angle, axisX, axisY, axisZ);
}

void applyMultiStateControlledRotateAroundAxis(Qureg qureg, int* ctrls, int* states, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTarget(qureg, ctrls, numCtrls, targ, __func__);
    validate_controlStates(states, numCtrls, __func__); // permits states==nullptr
    validate_rotationAxisNotZeroVector(axisX, axisY, axisZ, __func__);

    // defer division of vector norm to improve numerical accuracy
    qreal norm = std::sqrt(std::pow(axisX,2) + std::pow(axisY,2) + std::pow(axisZ,2)); // != 0

    // treat as generic 1-qubit matrix
    qreal c = std::cos(angle/2);
    qreal s = std::sin(angle/2);
    qcomp u11 = c - (s * axisZ * 1_i) / norm;
    qcomp u12 =   - (s * (axisY + axisX * 1_i)) / norm;
    qcomp u21 =     (s * (axisY - axisX * 1_i)) / norm;
    qcomp u22 = c + (s * axisZ * 1_i) / norm;
    auto matr = getCompMatr1({{u11,u12},{u21,u22}});

    // harmlessly re-validates, and checks unitarity of matr
    applyMultiStateControlledCompMatr1(qureg, ctrls, states, numCtrls, targ, matr);
}



/*
 * Pauli gadgets
 */

void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    localiser_statevec_anyCtrlPauliGadget(qureg, {}, {}, str, phase);
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

    /// @todo
    /// CRUCIAL NOTE:
    /// exp(theta I..I) might be algorithmically ok (I'm not sure), but it WILL NOT
    /// effect a global phase change of theta (I think). Should validate against this
    /// sitaution just in case, or make the doc extremely explicit

    qreal phase = util_getPhaseFromGateAngle(angle);
    auto ctrlVec = util_getVector(controls, numControls);
    auto stateVec = util_getVector(states, numControls); // empty if states==nullptr
    localiser_statevec_anyCtrlPauliGadget(qureg, ctrlVec, stateVec, str, phase);

    if (!qureg.isDensityMatrix)
        return;

    // conj(e^iXZ) = e^(-iXZ), but conj(Y)=-Y, so odd-Y undoes phase negation
    phase *= paulis_hasOddNumY(str) ? 1 : -1;
    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    str = paulis_getShiftedPauliStr(str, qureg.numQubits);
    localiser_statevec_anyCtrlPauliGadget(qureg, ctrlVec, stateVec, str, phase);
}



/*
 * phase gadgets
 */

void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    qreal phase = util_getPhaseFromGateAngle(angle);
    localiser_statevec_anyCtrlPhaseGadget(qureg, {}, {}, util_getVector(targets,numTargets), phase);
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

    qreal phase = util_getPhaseFromGateAngle(angle);
    auto ctrlVec = util_getVector(controls, numControls);
    auto targVec = util_getVector(targets,  numTargets);
    auto stateVec = util_getVector(states,  numControls); // empty if states==nullptr
    localiser_statevec_anyCtrlPhaseGadget(qureg, ctrlVec, stateVec, targVec, phase);

    if (!qureg.isDensityMatrix)
        return;

    phase *= -1;
    ctrlVec = util_getBraQubits(ctrlVec, qureg);
    targVec = util_getBraQubits(targVec, qureg);
    localiser_statevec_anyCtrlPhaseGadget(qureg, ctrlVec, stateVec, targVec, phase);
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
    DiagMatr1 matr = getDiagMatr1({1, std::exp(1_i * angle)});

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

void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    multiplyPauliStr(qureg, getPauliStr(std::string(numTargets, 'X'), targets, numTargets));
}

void applyMultiQubitNot(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    // harmlessly re-validates
    applyMultiStateControlledMultiQubitNot(qureg, nullptr, nullptr, 0, targets, numTargets);
}

void applyControlledMultiQubitNot(Qureg qureg, int control, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_controlAndTargets(qureg, control, targets, numTargets, __func__);

    // harmlessly re-validates
    applyMultiStateControlledMultiQubitNot(qureg, &control, nullptr, 1, targets, numTargets);
}

void applyMultiControlledMultiQubitNot(Qureg qureg, int* controls, int numControls, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);

    // harmlessly re-validates
    applyMultiStateControlledMultiQubitNot(qureg, controls, nullptr, numControls, targets, numTargets);
}

void applyMultiStateControlledMultiQubitNot(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_controlsAndTargets(qureg, controls, numControls, targets, numTargets, __func__);
    validate_controlStates(states, numControls, __func__);

    // treat as an all-X PauliStr
    PauliStr str = getPauliStr(std::string(numTargets, 'X'), targets, numTargets);

    // harmlessly re-validates
    applyMultiStateControlledPauliStr(qureg, controls, states, numControls, str);
}



/*
 * measurement
 */

int applyQubitMeasurement(Qureg qureg, int target) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    qreal prob = 0; // ignored
    return applyQubitMeasurementAndGetProb(qureg, target, &prob); // harmlessly re-validates
}

int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);

    // we do not assume state normalisation (that is posteriori checked),
    // so we must perform two reductions; one for each outcome. We choose
    // to re-enumerate the state (potentially doubling caching costs) to
    // avoid the nuisances/race-cons of parallel adding to two scalars.
    vector<qreal> probs(2);
    probs[0] = calcProbOfQubitOutcome(qureg, target, 0); // harmlessly re-validates
    probs[1] = calcProbOfQubitOutcome(qureg, target, 1); // " "
    validate_measurementProbsAreNormalised(probs, __func__);

    // randomly choose the outcome
    int outcome = rand_getRandomSingleQubitOutcome(probs[0]);
    *probability = probs[outcome];

    // collapse to the outcome
    (qureg.isDensityMatrix)?
        localiser_densmatr_multiQubitProjector(qureg, {target}, {outcome}, *probability):
        localiser_statevec_multiQubitProjector(qureg, {target}, {outcome}, *probability);

    return outcome;
}

qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__);

    // note that we do not merely invoke applyForcedMultiQubitMeasurement()
    // because we must validate the renormalising probability and
    // report this function's name during the error message
    qreal prob = calcProbOfQubitOutcome(qureg, target, outcome); // harmlessly re-validates
    validate_measurementOutcomeProbNotZero(outcome, prob, __func__);

    // project to the outcome, renormalising the surviving states
    (qureg.isDensityMatrix)?
        localiser_densmatr_multiQubitProjector(qureg, {target}, {outcome}, prob):
        localiser_statevec_multiQubitProjector(qureg, {target}, {outcome}, prob);

    return prob;
}

void applyQubitProjector(Qureg qureg, int target, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__); 
    
    // we permit the outcome to be negligibly likely, leaving state = null
    qreal prob = 1;
    (qureg.isDensityMatrix)?
        localiser_densmatr_multiQubitProjector(qureg, {target}, {outcome}, prob):
        localiser_statevec_multiQubitProjector(qureg, {target}, {outcome}, prob);
}

qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);

    qreal prob = 0; // ignored

    // below validates post-measurement and would report 'AndGetProb' function suffix. Eh!
    return applyMultiQubitMeasurementAndGetProb(qureg, qubits, numQubits, &prob);
}

qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);

    // find the probability of all possible outcomes
    qindex numProbs = powerOf2(numQubits);
    vector<qreal> probs(numProbs);
    calcProbsOfAllMultiQubitOutcomes(probs.data(), qureg, qubits, numQubits); // harmlessly re-validates

    // we cannot meaningfully sample these probs if not normalised
    validate_measurementProbsAreNormalised(probs, __func__);

    // randomly choose an outcome
    qindex outcome = rand_getRandomMultiQubitOutcome(probs);
    *probability = probs[outcome];

    // map outcome to individual qubit outcomes
    auto qubitVec = util_getVector(qubits, numQubits);
    auto outcomeVec = vector<int>(numQubits);
    getBitsFromInteger(outcomeVec.data(), outcome, numQubits);

    // project to the outcomes, renormalising the surviving states
    (qureg.isDensityMatrix)?
        localiser_densmatr_multiQubitProjector(qureg, qubitVec, outcomeVec, *probability):
        localiser_statevec_multiQubitProjector(qureg, qubitVec, outcomeVec, *probability);

    return outcome;
}

qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesAreValid(outcomes, numQubits, __func__);

    auto qubitVec = util_getVector(qubits, numQubits);
    auto outcomeVec = util_getVector(outcomes, numQubits);

    // ensure probability of the forced measurement outcome is not negligible
    qreal prob = calcProbOfMultiQubitOutcome(qureg, qubits, outcomes, numQubits); // harmlessly re-validates
    validate_measurementOutcomesProbNotZero(outcomes, numQubits, prob, __func__);

    // project to the outcome, renormalising the surviving states
    (qureg.isDensityMatrix)?
        localiser_densmatr_multiQubitProjector(qureg, qubitVec, outcomeVec, prob):
        localiser_statevec_multiQubitProjector(qureg, qubitVec, outcomeVec, prob);

    return prob;
}

void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesAreValid(outcomes, numQubits, __func__);

    qreal prob = 1;
    auto qubitVec = util_getVector(qubits, numQubits);
    auto outcomeVec = util_getVector(outcomes, numQubits);

    (qureg.isDensityMatrix)?
        localiser_densmatr_multiQubitProjector(qureg, qubitVec, outcomeVec, prob):
        localiser_statevec_multiQubitProjector(qureg, qubitVec, outcomeVec, prob);
}



/*
 * QFT
 */

void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, targets, numTargets, __func__);

    /// @todo
    /// change this placeholder implementation to the bespoke, optimised routine,
    /// wherein each contiguous controlled-phase gate is merged

    for (int n=numTargets-1; n>=0; n--) {
        applyHadamard(qureg, targets[n]);
        for (int m=0; m<n; m++) {
            qreal arg = const_PI / powerOf2(m+1);
            applyTwoQubitPhaseShift(qureg, targets[n], targets[n-m-1], arg);
        }
    }

    int mid = numTargets/2; // floors
    for (int n=0; n<mid; n++)
        applySwap(qureg, targets[n], targets[numTargets-1-n]);
}

void applyFullQuantumFourierTransform(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    vector<int> targets(qureg.numQubits);
    for (size_t i=0; i<targets.size(); i++)
        targets[i] = i;

    applyQuantumFourierTransform(qureg, targets.data(), targets.size());
}



} // end de-mangler