/** @file
 * API definitions for effecting operators (including
 * unitaries, projectors, channels, Hermitians, and
 * arbitrary matrices) upon Quregs, which can be
 * statevectors or density matrices
 */

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "../core/validation.hpp"
#include "../core/utilities.hpp"
#include "../core/localiser.hpp"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions


// enable invocation by both C and C++ binaries
extern "C" {


#define _NOT_IMPLEMENTED_ERROR_DEF { error_functionNotImplemented(__func__); }
/*
 * Named gates
 */

void applyS(Qureg qureg, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyT(Qureg qureg, int target)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyHadamard(Qureg qureg, int target)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * CompMatr1
 */

void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyCompMatr1(Qureg qureg, int target, CompMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * CompMatr2
 */

void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * CompMatr
 */

void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matr)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * DiagMatr1
 */

void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matr)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * DiagMatr2
 */

void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matr)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * DiagMatr
 */

void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyControlledDiagMatr(Qureg, int control, int* targets, int numTargets, DiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiControlledDiagMatr(Qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF

void applyMultiStateControlledDiagMatr(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF



/*
 * FullStateDiagMatr
 */

void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF

void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix)
    _NOT_IMPLEMENTED_ERROR_DEF



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

void applySuperOp(Qureg qureg, SuperOp superop) { 

    // we only ever left-apply it (of course), so it
    // has no equivalent multiplySuperOp
    
    // TODO
    error_functionNotImplemented(__func__);
}



/*
 * measurement
 */

int applyQubitMeasurement(Qureg qureg, int target) {

    // // TODO
    // error_functionNotImplemented(__func__);
    return -1;
}

int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability) {

    // // TODO
    // error_functionNotImplemented(__func__);
    return -1;
}

void applyForcedQubitMeasurement(Qureg qureg, int target, int outcome)
    _NOT_IMPLEMENTED_ERROR_DEF

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