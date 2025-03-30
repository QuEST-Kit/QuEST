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

/// @notdoced
/// @notvalidated
void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matr);

/// @notdoced
/// @notvalidated
void applyCompMatr1(Qureg qureg, int target, CompMatr1 matr);

/// @notdoced
/// @notvalidated
void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matr);

/// @notdoced
/// @notvalidated
void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matr);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matr);



/*
 * CompMatr2
 */

/// @notdoced
/// @notvalidated
void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);

/// @notdoced
/// @notvalidated
void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);

/// @notdoced
/// @notvalidated
void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matr);

/// @notdoced
/// @notvalidated
void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matr);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matr);



/*
 * CompMatr
 */

/// @notdoced
/// @notvalidated
void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);

/// @notdoced
/// @notvalidated
void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);

/// @notdoced
/// @notvalidated
void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matr);

/// @notdoced
/// @notvalidated
void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matr);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matr);



/*
 * DiagMatr1
 */

/// @notdoced
/// @notvalidated
void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);

/// @notdoced
/// @notvalidated
void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);

/// @notdoced
/// @notvalidated
void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matr);

/// @notdoced
/// @notvalidated
void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matr);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matr);



/*
 * DiagMatr2
 */

/// @notdoced
/// @notvalidated
void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);

/// @notdoced
/// @notvalidated
void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);

/// @notdoced
/// @notvalidated
void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matr);

/// @notdoced
/// @notvalidated
void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matr);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matr);



/*
 * DiagMatr
 */

/// @notdoced
/// @notvalidated
void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);

/// @notdoced
/// @notvalidated
void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);

/// @notdoced
/// @notvalidated
void applyControlledDiagMatr(Qureg, int control, int* targets, int numTargets, DiagMatr matrix);

/// @notdoced
/// @notvalidated
void applyMultiControlledDiagMatr(Qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledDiagMatr(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix);



/*
 * DiagMatrPower
 */

/// @notdoced
/// @notvalidated
void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notdoced
/// @notvalidated
void applyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notdoced
/// @notvalidated
void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notdoced
/// @notvalidated
void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledDiagMatrPower(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);



/*
 * FullStateDiagMatr
 */

/// @notdoced
/// @notvalidated
void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);

/// @notdoced
/// @notvalidated
void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);

/// @notdoced
/// @notvalidated
void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);

/// @notdoced
/// @notvalidated
void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);



/*
 * S gate
 */

/// @notdoced
/// @notvalidated
void applyS(Qureg qureg, int target);

/// @notdoced
/// @notvalidated
void applyControlledS(Qureg qureg, int control, int target);

/// @notdoced
/// @notvalidated
void applyMultiControlledS(Qureg qureg, int* controls, int numControls, int target);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledS(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * T gate
 */

/// @notdoced
/// @notvalidated
void applyT(Qureg qureg, int target);

/// @notdoced
/// @notvalidated
void applyControlledT(Qureg qureg, int control, int target);

/// @notdoced
/// @notvalidated
void applyMultiControlledT(Qureg qureg, int* controls, int numControls, int target);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledT(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * Hadamard 
 */

/// @notdoced
/// @notvalidated
void applyHadamard(Qureg qureg, int target);

/// @notdoced
/// @notvalidated
void applyControlledHadamard(Qureg qureg, int control, int target);

/// @notdoced
/// @notvalidated
void applyMultiControlledHadamard(Qureg qureg, int* controls, int numControls, int target);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledHadamard(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * swaps
 */

/// @notdoced
/// @notvalidated
void multiplySwap(Qureg qureg, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applySwap(Qureg qureg, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);



/*
 * sqrt-swap
 */

/// @notdoced
/// @notvalidated
void applySqrtSwap(Qureg qureg, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applyControlledSqrtSwap(Qureg qureg, int control, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applyMultiControlledSqrtSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledSqrtSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);



/*
 * individual Paulis
 */

/// @notdoced
/// @notvalidated
void multiplyPauliX(Qureg qureg, int target);
/// @notdoced
/// @notvalidated
void multiplyPauliY(Qureg qureg, int target);
/// @notdoced
/// @notvalidated
void multiplyPauliZ(Qureg qureg, int target);

/// @notdoced
/// @notvalidated
void applyPauliX(Qureg qureg, int target);
/// @notdoced
/// @notvalidated
void applyPauliY(Qureg qureg, int target);
/// @notdoced
/// @notvalidated
void applyPauliZ(Qureg qureg, int target);

/// @notdoced
/// @notvalidated
void applyControlledPauliX(Qureg qureg, int control, int target);
/// @notdoced
/// @notvalidated
void applyControlledPauliY(Qureg qureg, int control, int target);
/// @notdoced
/// @notvalidated
void applyControlledPauliZ(Qureg qureg, int control, int target);

/// @notdoced
/// @notvalidated
void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target);
/// @notdoced
/// @notvalidated
void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target);
/// @notdoced
/// @notvalidated
void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target);
/// @notdoced
/// @notvalidated
void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target);
/// @notdoced
/// @notvalidated
void applyMultiStateControlledPauliZ(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * Pauli strings
 */

/// @notdoced
/// @notvalidated
void multiplyPauliStr(Qureg qureg, PauliStr str);

/// @notdoced
/// @notvalidated
void applyPauliStr(Qureg qureg, PauliStr str);

/// @notdoced
/// @notvalidated
void applyControlledPauliStr(Qureg qureg, int control, PauliStr str);

/// @notdoced
/// @notvalidated
void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledPauliStr(Qureg qureg, int* controls, int* states, int numControls, PauliStr str);



/*
 * Pauli string sums
 */

/// @notdoced
/// @notvalidated
void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);

/// @nottested
void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);



/*
 * individual axis rotations
 */

// don't think users will ever want to left-multiply only

/// @notdoced
/// @notvalidated
void applyRotateX(Qureg qureg, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyRotateY(Qureg qureg, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyRotateZ(Qureg qureg, int target, qreal angle);

/// @notdoced
/// @notvalidated
void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);
/// @notdoced
/// @notvalidated
void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);



/*
 * arbitrary axis rotation
 */

/// @notdoced
/// @notvalidated
void applyRotateAroundAxis(Qureg qureg, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

/// @notdoced
/// @notvalidated
void applyControlledRotateAroundAxis(Qureg qureg, int ctrl, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

/// @notdoced
/// @notvalidated
void applyMultiControlledRotateAroundAxis(Qureg qureg, int* ctrls, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledRotateAroundAxis(Qureg qureg, int* ctrls, int* states, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);



/*
 * Pauli gadgets
 */

/// @notdoced
/// @notvalidated
void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);

/// @notdoced
/// @notvalidated
void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle);

/// @notdoced
/// @notvalidated
void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledPauliGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStr str, qreal angle);



/*
 * phase gadgets
 */

/// @notdoced
/// @notvalidated
void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);

/// @notdoced
/// @notvalidated
void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);

/// @notdoced
/// @notvalidated
void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledPhaseGadget(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, qreal angle);



/*
 * phase shifts and flips
 */

/// @notdoced
/// @notvalidated
void applyPhaseFlip (Qureg qureg, int target);
/// @notdoced
/// @notvalidated
void applyPhaseShift(Qureg qureg, int target, qreal angle);

/// @notdoced
/// @notvalidated
void applyTwoQubitPhaseFlip( Qureg qureg, int target1, int target2);
/// @notdoced
/// @notvalidated
void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle);

/// @notdoced
/// @notvalidated
void applyMultiQubitPhaseFlip (Qureg qureg, int* targets, int numTargets);
/// @notdoced
/// @notvalidated
void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle);



/*
 * many-qubit CNOTs (aliases for X)
 */

/// @notdoced
/// @notvalidated
void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);

/// @notdoced
/// @notvalidated
void applyMultiQubitNot(Qureg, int* targets, int numTargets);

/// @notdoced
/// @notvalidated
void applyControlledMultiQubitNot(Qureg, int control, int* targets, int numTargets);

/// @notdoced
/// @notvalidated
void applyMultiControlledMultiQubitNot(Qureg, int* controls, int numControls, int* targets, int numTargets);

/// @notdoced
/// @notvalidated
void applyMultiStateControlledMultiQubitNot(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets);



/*
 * measurement
 */

/// @notdoced
/// @notvalidated
int applyQubitMeasurement(Qureg qureg, int target);

/// @notdoced
/// @notvalidated
int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability);

/// @notdoced
/// @notvalidated
qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome);

/// @notdoced
/// @notvalidated
void applyQubitProjector(Qureg qureg, int target, int outcome);

/// @notdoced
/// @notvalidated
qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);

/// @notdoced
/// @notvalidated
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);

/// @notdoced
/// @notvalidated
qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);

/// @notdoced
/// @notvalidated
void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);



/*
 * QFT
 */

/// @notdoced
/// @notvalidated
void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets);

/// @notdoced
/// @notvalidated
void applyFullQuantumFourierTransform(Qureg qureg);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // OPERATIONS_H

/** @} (end doxygen defgroup) */
