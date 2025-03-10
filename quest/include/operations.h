/** @file
 * API signatures for effecting operators (such as gates and unitaries) 
 * upon Quregs which are instantiated as either statevectors or 
 * density matrices. This excludes decoherence channels which are
 * instead exposed in decoherence.h
 * 
 * @author Tyson Jones
 * @defgroup operations Operations
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

void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matr);

void applyCompMatr1(Qureg qureg, int target, CompMatr1 matr);

void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matr);

void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matr);

void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matr);



/*
 * CompMatr2
 */

void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);

void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);

void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matr);

void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matr);

void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matr);



/*
 * CompMatr
 */

void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);

void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);

void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matr);

void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matr);

void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matr);



/*
 * DiagMatr1
 */

void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);

void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);

void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matr);

void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matr);

void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matr);



/*
 * DiagMatr2
 */

void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);

void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);

void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matr);

void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matr);

void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matr);



/*
 * DiagMatr
 */

void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);

void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);

void applyControlledDiagMatr(Qureg, int control, int* targets, int numTargets, DiagMatr matrix);

void applyMultiControlledDiagMatr(Qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix);

void applyMultiStateControlledDiagMatr(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix);



/*
 * DiagMatrPower
 */

void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

void applyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);

void applyMultiStateControlledDiagMatrPower(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);



/*
 * FullStateDiagMatr
 */

void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);

void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);

void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);

void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);



/*
 * S gate
 */

void applyS(Qureg qureg, int target);

void applyControlledS(Qureg qureg, int control, int target);

void applyMultiControlledS(Qureg qureg, int* controls, int numControls, int target);

void applyMultiStateControlledS(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * T gate
 */

void applyT(Qureg qureg, int target);

void applyControlledT(Qureg qureg, int control, int target);

void applyMultiControlledT(Qureg qureg, int* controls, int numControls, int target);

void applyMultiStateControlledT(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * Hadamard 
 */

void applyHadamard(Qureg qureg, int target);

void applyControlledHadamard(Qureg qureg, int control, int target);

void applyMultiControlledHadamard(Qureg qureg, int* controls, int numControls, int target);

void applyMultiStateControlledHadamard(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * swaps
 */

void multiplySwap(Qureg qureg, int qubit1, int qubit2);

void applySwap(Qureg qureg, int qubit1, int qubit2);

void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2);

void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);

void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);



/*
 * sqrt-swap
 */

void applySqrtSwap(Qureg qureg, int qubit1, int qubit2);

void applyControlledSqrtSwap(Qureg qureg, int control, int qubit1, int qubit2);

void applyMultiControlledSqrtSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);

void applyMultiStateControlledSqrtSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);



/*
 * individual Paulis
 */

void multiplyPauliX(Qureg qureg, int target);
void multiplyPauliY(Qureg qureg, int target);
void multiplyPauliZ(Qureg qureg, int target);

void applyPauliX(Qureg qureg, int target);
void applyPauliY(Qureg qureg, int target);
void applyPauliZ(Qureg qureg, int target);

void applyControlledPauliX(Qureg qureg, int control, int target);
void applyControlledPauliY(Qureg qureg, int control, int target);
void applyControlledPauliZ(Qureg qureg, int control, int target);

void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target);
void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target);
void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target);

void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target);
void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target);
void applyMultiStateControlledPauliZ(Qureg qureg, int* controls, int* states, int numControls, int target);



/*
 * Pauli strings
 */

void multiplyPauliStr(Qureg qureg, PauliStr str);

void applyPauliStr(Qureg qureg, PauliStr str);

void applyControlledPauliStr(Qureg qureg, int control, PauliStr str);

void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str);

void applyMultiStateControlledPauliStr(Qureg qureg, int* controls, int* states, int numControls, PauliStr str);



/*
 * Pauli string sums
 */

void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);



/*
 * individual axis rotations
 */

// don't think users will ever want to left-multiply only

void applyRotateX(Qureg qureg, int target, qreal angle);
void applyRotateY(Qureg qureg, int target, qreal angle);
void applyRotateZ(Qureg qureg, int target, qreal angle);

void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle);
void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle);
void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle);

void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle);
void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle);
void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle);

void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);
void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);
void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);



/*
 * arbitrary axis rotation
 */

void applyRotateAroundAxis(Qureg qureg, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

void applyControlledRotateAroundAxis(Qureg qureg, int ctrl, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

void applyMultiControlledRotateAroundAxis(Qureg qureg, int* ctrls, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);

void applyMultiStateControlledRotateAroundAxis(Qureg qureg, int* ctrls, int* states, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);



/*
 * Pauli gadgets
 */

void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);

void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle);

void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle);

void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle);

void applyMultiStateControlledPauliGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStr str, qreal angle);



/*
 * phase gadgets
 */

void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);

void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);

void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle);

void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle);

void applyMultiStateControlledPhaseGadget(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, qreal angle);



/*
 * phase shifts and flips
 */

void applyPhaseFlip (Qureg qureg, int target);
void applyPhaseShift(Qureg qureg, int target, qreal angle);

void applyTwoQubitPhaseFlip( Qureg qureg, int target1, int target2);
void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle);

void applyMultiQubitPhaseFlip (Qureg qureg, int* targets, int numTargets);
void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle);



/*
 * many-qubit CNOTs (aliases for X)
 */

void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);

void applyMultiQubitNot(Qureg, int* targets, int numTargets);

void applyControlledMultiQubitNot(Qureg, int control, int* targets, int numTargets);

void applyMultiControlledMultiQubitNot(Qureg, int* controls, int numControls, int* targets, int numTargets);

void applyMultiStateControlledMultiQubitNot(Qureg, int* controls, int* states, int numControls, int* targets, int numTargets);



/*
 * superoperator
 */

void applySuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop);



/*
 * measurement
 */

int applyQubitMeasurement(Qureg qureg, int target);

int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability);

qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome);

void applyQubitProjector(Qureg qureg, int target, int outcome);

qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);

qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);

qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);

void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);



/*
 * QFT
 */

void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets);

void applyFullQuantumFourierTransform(Qureg qureg);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // OPERATIONS_H

/** @} (end doxygen defgroup) */
