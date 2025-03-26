/** @file
 * API signatures for effecting operators (such as gates and unitaries) 
 * upon Quregs which are instantiated as either statevectors or 
 * density matrices. This excludes decoherence channels which are
 * instead exposed in decoherence.h
 * 
 * @author Tyson Jones
 * 
 * @defgroup operations Operations
 * @ingroup api
 * @brief Functions for effecting operators upon Quregs.
 * @{
 */

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/*
 * CompMatr1
 */

/// @notvalidated
void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matr);

/// @notvalidated
void applyCompMatr1(Qureg qureg, int target, CompMatr1 matr);

/// @notvalidated
void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matr);

/// @notvalidated
void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matr);

/// @notvalidated
void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matr);



/*
 * CompMatr2
 */

/// @notvalidated
void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);

/// @notvalidated
void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);

/// @notvalidated
void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matr);

/// @notvalidated
void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matr);

/// @notvalidated
void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matr);



/*
 * CompMatr
 */

/// @notvalidated
void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);

/// @notvalidated
void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);

/// @notvalidated
void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matr);

/// @notvalidated
void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matr);

/// @notvalidated
void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matr);



/*
 * DiagMatr1
 */

/// @notvalidated
void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);

/// @notvalidated
void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);

/// @notvalidated
void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matr);

/// @notvalidated
void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matr);

/// @notvalidated
void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matr);



/*
 * DiagMatr2
 */

/// @notvalidated
void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);

/// @notvalidated
void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);

/// @notvalidated
void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matr);

/// @notvalidated
void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matr);

/// @notvalidated
void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matr);



/*
 * DiagMatr
 */

/// @notvalidated
void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);

/// @notvalidated
void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);

/// @notvalidated
void applyControlledDiagMatr(Qureg, int control, int* targets, int numTargets, DiagMatr matrix);

/// @notvalidated
void applyMultiControlledDiagMatr(Qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix);

/// @notvalidated
void applyMultiStateControlledDiagMatr(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix);



/*
 * DiagMatrPower
 */

/// @notvalidated
void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notvalidated
void applyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notvalidated
void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notvalidated
void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notvalidated
void applyMultiStateControlledDiagMatrPower(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);



/*
 * FullStateDiagMatr
 */

/// @notvalidated
void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);

/// @notvalidated
void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);

/// @notvalidated
void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);

/// @notvalidated
void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);



/*
 * S gate
 */

/// @notvalidated
void applyS(Qureg qureg, int target);

/// @notvalidated
void applyControlledS(Qureg qureg, int control, int target);

/// @notvalidated
void applyMultiControlledS(Qureg qureg, int* controls, int numControls, int target);

/// @notvalidated
void applyMultiStateControlledS(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * T gate
 */

/// @notvalidated
void applyT(Qureg qureg, int target);

/// @notvalidated
void applyControlledT(Qureg qureg, int control, int target);

/// @notvalidated
void applyMultiControlledT(Qureg qureg, int* controls, int numControls, int target);

/// @notvalidated
void applyMultiStateControlledT(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * Hadamard 
 */

/// @notvalidated
void applyHadamard(Qureg qureg, int target);

/// @notvalidated
void applyControlledHadamard(Qureg qureg, int control, int target);

/// @notvalidated
void applyMultiControlledHadamard(Qureg qureg, int* controls, int numControls, int target);

/// @notvalidated
void applyMultiStateControlledHadamard(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * swaps
 */

/// @notvalidated
void multiplySwap(Qureg qureg, int qubit1, int qubit2);

/// @notvalidated
void applySwap(Qureg qureg, int qubit1, int qubit2);

/// @notvalidated
void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2);

/// @notvalidated
void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);

/// @notvalidated
void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);



/*
 * sqrt-swap
 */

/// @notvalidated
void applySqrtSwap(Qureg qureg, int qubit1, int qubit2);

/// @notvalidated
void applyControlledSqrtSwap(Qureg qureg, int control, int qubit1, int qubit2);

/// @notvalidated
void applyMultiControlledSqrtSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);

/// @notvalidated
void applyMultiStateControlledSqrtSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);



/*
 * individual Paulis
 */

/// @notvalidated
void multiplyPauliX(Qureg qureg, int target);
/// @notvalidated
void multiplyPauliY(Qureg qureg, int target);
/// @notvalidated
void multiplyPauliZ(Qureg qureg, int target);

/// @notvalidated
void applyPauliX(Qureg qureg, int target);
/// @notvalidated
void applyPauliY(Qureg qureg, int target);
/// @notvalidated
void applyPauliZ(Qureg qureg, int target);

/// @notvalidated
void applyControlledPauliX(Qureg qureg, int control, int target);
/// @notvalidated
void applyControlledPauliY(Qureg qureg, int control, int target);
/// @notvalidated
void applyControlledPauliZ(Qureg qureg, int control, int target);

/// @notvalidated
void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target);
/// @notvalidated
void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target);
/// @notvalidated
void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target);

/// @notvalidated
void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target);
/// @notvalidated
void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target);
/// @notvalidated
void applyMultiStateControlledPauliZ(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * Pauli strings
 */

/// @notvalidated
void multiplyPauliStr(Qureg qureg, PauliStr str);

/// @notvalidated
void applyPauliStr(Qureg qureg, PauliStr str);

/// @notvalidated
void applyControlledPauliStr(Qureg qureg, int control, PauliStr str);

/// @notvalidated
void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str);

/// @notvalidated
void applyMultiStateControlledPauliStr(Qureg qureg, int* controls, int* states, int numControls, PauliStr str);



/*
 * Pauli string sums
 */

/// @notvalidated
void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);

/// @nottested
void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);



/*
 * individual axis rotations
 */

// don't think users will ever want to left-multiply only

/// @notvalidated
void applyRotateX(Qureg qureg, int target, qreal angle);
/// @notvalidated
void applyRotateY(Qureg qureg, int target, qreal angle);
/// @notvalidated
void applyRotateZ(Qureg qureg, int target, qreal angle);

/// @notvalidated
void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle);
/// @notvalidated
void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle);
/// @notvalidated
void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle);

/// @notvalidated
void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle);
/// @notvalidated
void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle);
/// @notvalidated
void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle);

/// @notvalidated
void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);
/// @notvalidated
void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);
/// @notvalidated
void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);



/*
 * arbitrary axis rotation
 */

/// @notvalidated
void applyRotateAroundAxis(Qureg qureg, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

/// @notvalidated
void applyControlledRotateAroundAxis(Qureg qureg, int ctrl, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

/// @notvalidated
void applyMultiControlledRotateAroundAxis(Qureg qureg, int* ctrls, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

/// @notvalidated
void applyMultiStateControlledRotateAroundAxis(Qureg qureg, int* ctrls, int* states, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);



/*
 * Pauli gadgets
 */

/// @notvalidated
void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);

/// @notvalidated
void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle);

/// @notvalidated
void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle);

/// @notvalidated
void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle);

/// @notvalidated
void applyMultiStateControlledPauliGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStr str, qreal angle);



/*
 * phase gadgets
 */

/// @notvalidated
void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);

/// @notvalidated
void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);

/// @notvalidated
void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle);

/// @notvalidated
void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle);

/// @notvalidated
void applyMultiStateControlledPhaseGadget(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, qreal angle);



/*
 * phase shifts and flips
 */

/// @notvalidated
void applyPhaseFlip (Qureg qureg, int target);
/// @notvalidated
void applyPhaseShift(Qureg qureg, int target, qreal angle);

/// @notvalidated
void applyTwoQubitPhaseFlip( Qureg qureg, int target1, int target2);
/// @notvalidated
void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle);

/// @notvalidated
void applyMultiQubitPhaseFlip (Qureg qureg, int* targets, int numTargets);
/// @notvalidated
void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle);



/*
 * many-qubit CNOTs (aliases for X)
 */

/// @notvalidated
void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);

/// @notvalidated
void applyMultiQubitNot(Qureg, int* targets, int numTargets);

/// @notvalidated
void applyControlledMultiQubitNot(Qureg, int control, int* targets, int numTargets);

/// @notvalidated
void applyMultiControlledMultiQubitNot(Qureg, int* controls, int numControls, int* targets, int numTargets);

/// @notvalidated
void applyMultiStateControlledMultiQubitNot(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets);



/*
 * measurement
 */

/// @notvalidated
int applyQubitMeasurement(Qureg qureg, int target);

/// @notvalidated
int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability);

/// @notvalidated
qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome);

/// @notvalidated
void applyQubitProjector(Qureg qureg, int target, int outcome);

/// @notvalidated
qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);

/// @notvalidated
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);

/// @notvalidated
qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);

/// @notvalidated
void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);



/*
 * QFT
 */

/// @notvalidated
void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets);

/// @notvalidated
void applyFullQuantumFourierTransform(Qureg qureg);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // OPERATIONS_H

/** @} (end doxygen defgroup) */
