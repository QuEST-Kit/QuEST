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



/** 
 * @defgroup op_compmatr1 CompMatr1
 * @brief Functions for applying general one-qubit dense matrices, as CompMatr1.
 * @{
 */


/** Multiplies a general one-qubit dense @p matrix upon the specified @p target 
 * qubit of @p qureg.
 *  
 * @formulae
 * Let @f$ \hat{M} = @f$ @p matrix and @f$ t = @f$ @p target, and notate 
 * @f$\hat{M}_t@f$ as per applyCompMatr1(). Unlike applyCompMatr1() however,
 * this function only ever left-multiplies @p matrix upon @p qureg, regardless
 * of whether it is a statevector or density matrix.
 * 
 * Explicitly,
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function effects
 *   @f[ 
        \svpsi \rightarrow \hat{M}_t \, \svpsi.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function effects
 *   @f[ 
        \dmrho \rightarrow \hat{M}_t \, \dmrho.
 *   @f]
 *
 * There are no additional constraints like unitarity.
 *
 * @myexample
 * ```
    Qureg qureg = createDensityQureg(5);

    CompMatr1 matrix = getInlineCompMatr1({
        {0.1, 0.2},
        {0.3i, 0.4i}
    });

    multiplyCompMatr1(qureg, 2, matrix); 
 * ```
 *
 * @param[in,out] qureg  the state to modify.
 * @param[in]     target the index of the target qubit.
 * @param[in]     matrix the Z-basis matrix to multiply.
 * @throws invalidQuESTInputError()
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p target is an invalid qubit index.
 * @notvalidated
 * @see
 * - getCompMatr1()
 * - getInlineCompMatr1()
 * - applyCompMatr1()
 * - applyQubitProjector()
 * - multiplyCompMatr()
 * @author Tyson Jones
 */
void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matr);


/** Applies a general one-qubit dense unitary @p matrix to the specified @p target 
 * qubit of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  wireL [shape=plaintext, label="target"];
  wireR [shape=plaintext, label=""];
  gate  [shape=box,    label="matrix"];

  wireL -> gate -> wireR
}
 * @enddot
 * 
 * @formulae
 * Let @f$ \hat{U} = @f$ @p matrix, @f$ t = @f$ @p target, and let @f$\hat{U}_t@f$
 * notate operating @f$\hat{U}@f$ upon the @f$ t @f$-th qubit among@f$ N @f$, i.e.
 * @f[ 
        \hat{U}_t \equiv \id^{N-t} \otimes \hat{U} \otimes \id^{t-1}.
 * @f]
 * Then,
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function effects
 *   @f[ 
        \svpsi \rightarrow \hat{U}_t \, \svpsi.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function effects
 *   @f[ 
        \dmrho \rightarrow \hat{U}_t \, \dmrho \, {\hat{U}_t}^\dagger.
 *   @f]
 *
 * @constraints
 * - Unitarity of @f$ \hat{U} = @f$ @p matrix requires that 
 *   @f$ \hat{U} \hat{U}^\dagger = \id @f$. Validation will check that @p matrix is
 *   approximately unitarity via
 *   @f[ 
        \max\limits_{ij} \Big|\left(\hat{U} \hat{U}^\dagger - \id\right)_{ij}\Big|^2 \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setValidationEpsilon().
 * 
 * @myexample
 * ```
    Qureg qureg = createQureg(5);

    CompMatr1 matrix = getInlineCompMatr1({
        {-1i/sqrt(2), 1i/sqrt(2)},
        {(1i-1)/2,    (1i-1)/2}
    });

    applyCompMatr1(qureg, 2, matrix); 
 * ```
 *
 * @param[in,out] qureg  the state to modify.
 * @param[in]     target the index of the target qubit.
 * @param[in]     matrix the Z-basis unitary matrix to effect.
 * @throws invalidQuESTInputError()
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix is not approximately unitary.
 * - if @p target is an invalid qubit index.
 * @notvalidated
 * @see
 * - getCompMatr1()
 * - getInlineCompMatr1()
 * - multiplyCompMatr1()
 * - applyControlledCompMatr1()
 * - applyCompMatr2()
 * - applyCompMatr()
 * @author Tyson Jones
 */
void applyCompMatr1(Qureg qureg, int target, CompMatr1 matrix);


/** Applies a singly-controlled one-qubit dense unitary @p matrix to the specified 
 * @p target qubit of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  topWireL [shape=plaintext, label="control"];
  topWireR [shape=plaintext, label=""];
  ctrl  [shape=circle, label="", width=.12, style=filled, fillcolor=black];

  topWireL -> ctrl -> topWireR;

  botWireL [shape=plaintext, label="target"];
  botWireR [shape=plaintext, label=""];
  gate  [shape=box,    label="matrix"];

  botWireL -> gate -> botWireR;
  ctrl -> gate;

  {rank=same; topWireL; botWireL};
  {rank=same; ctrl;     gate};
  {rank=same; topWireR; botWireR};
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matrix);


/** Applies a multiply-controlled one-qubit dense unitary @p matrix to the specified 
 * @p target qubit of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  trailingCtrl [shape=plaintext, label="..."];

  topWireL [shape=plaintext, label="controls[1]"];
  topWireR [shape=plaintext, label=""];
  topCtrl  [shape=circle, label="", width=.12, style=filled, fillcolor=black];

  topWireL -> topCtrl -> topWireR;

  midWireL [shape=plaintext, label="controls[0]"];
  midWireR [shape=plaintext, label=""];
  midCtrl  [shape=circle, label="", width=.12, style=filled, fillcolor=black];

  midWireL -> midCtrl -> midWireR;

  botWireL [shape=plaintext, label="target"];
  botWireR [shape=plaintext, label=""];
  gate  [shape=box,    label="matrix"];

  botWireL -> gate -> botWireR;
  trailingCtrl -> topCtrl -> midCtrl -> gate;

  {rank=same; topWireL; midWireL; botWireL};
  {rank=same; trailingCtrl; topCtrl; midCtrl; gate};
  {rank=same; topWireR; midWireR; botWireR};
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matrix);


/** Applies an arbitrarily-controlled one-qubit dense unitary @p matrix to the specified 
 * @p target qubit of @p qureg, conditioned upon the @p controls being in the given @p states.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  trailingCtrl [shape=plaintext, label="..."];

  topWireL [shape=plaintext, label="controls[1]"];
  topWireR [shape=plaintext, label=""];
  topCtrl  [shape=circle, label="", width=.12, style=filled, fillcolor=black];

  topWireL -> topCtrl -> topWireR;

  midWireL [shape=plaintext, label="controls[0]"];
  midWireR [shape=plaintext, label=""];
  midCtrl  [shape=circle, label="", width=.12, style=filled, fillcolor=white];

  midWireL -> midCtrl -> midWireR;

  botWireL [shape=plaintext, label="target"];
  botWireR [shape=plaintext, label=""];
  gate  [shape=box,    label="matrix"];

  botWireL -> gate -> botWireR;
  trailingCtrl -> topCtrl -> midCtrl -> gate;

  {rank=same; topWireL; midWireL; botWireL};
  {rank=same; trailingCtrl; topCtrl; midCtrl; gate};
  {rank=same; topWireR; midWireR; botWireR};
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matrix);


/** @} */



/** 
 * @defgroup op_compmatr2 CompMatr2
 * @brief Functions for applying general two-qubit dense matrices, as CompMatr2.
 * @{
 */


/// @notdoced
/// @notvalidated
void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);


/** Applies a general two-qubit dense unitary @p matrix to qubits @p target1 and
 * @p target2 (treated as increasing significance) of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  layout=neato;
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  topWireL [shape=plaintext, pos="0,0!", label="target2"];
  topWireR [shape=plaintext, pos="2.5,0!", label=""];

  botWireL [shape=plaintext, pos="0,.5!", label="target1"];
  botWireR [shape=plaintext, pos="2.5,.5!", label=""];

  gate  [shape=rectangle, label="matrix", style=filled, fillcolor=white, height=1, pos="1.25,.25!"];

  topWireL -> topWireR;
  botWireL -> botWireR;
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix);


/** Applies a singly-controlled two-qubit dense unitary @p matrix to qubits 
 * @p target1 and @p target2 (treated as increasing significance) of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  layout=neato;
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  topWireL [shape=plaintext, pos="0,1!", label="control"];
  topWireR [shape=plaintext, pos="2.5,1!", label=""];

  midWireL [shape=plaintext, pos="0,0.5!", label="target2"];
  midWireR [shape=plaintext, pos="2.5,0.5!", label=""];

  botWireL [shape=plaintext, pos="0,0!", label="target1"];
  botWireR [shape=plaintext, pos="2.5,0!", label=""];

  gate [shape=rectangle, label="matrix", style=filled, fillcolor=white, height=1, pos="1.25,0.25!"];
  ctrl [shape=circle, label="", width=.12, style=filled, fillcolor=black, pos="1.25,1!"];

  topWireL -> ctrl -> topWireR;
  midWireL -> midWireR;
  botWireL -> botWireR;
  ctrl -> gate;
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matr);


/** Applies a multiply-controlled two-qubit dense unitary @p matrix to qubits 
 * @p target1 and @p target2 (treated as increasing significance) of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  layout=neato;
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  tippytopWireL [shape=plaintext, pos="0,1.5!", label="controls[1]"];
  tippytopWireR [shape=plaintext, pos="2.5,1.5!", label=""];

  topWireL [shape=plaintext, pos="0,1!", label="controls[0]"];
  topWireR [shape=plaintext, pos="2.5,1!", label=""];

  midWireL [shape=plaintext, pos="0,0.5!", label="target2"];
  midWireR [shape=plaintext, pos="2.5,0.5!", label=""];

  botWireL [shape=plaintext, pos="0,0!", label="target1"];
  botWireR [shape=plaintext, pos="2.5,0!", label=""];

  gate [shape=rectangle, label="matrix", style=filled, fillcolor=white, height=1, pos="1.25,0.25!"];
  ctrl1 [shape=circle, label="", width=.12, style=filled, fillcolor=black, pos="1.25,1!"];
  ctrl2 [shape=circle, label="", width=.12, style=filled, fillcolor=black, pos="1.25,1.5!"];
  trailingCtrl [shape=plaintext, label="...", pos="1.25,2!"];

  tippytopWireL -> ctrl2 -> tippytopWireR;
  topWireL -> ctrl1 -> topWireR;
  midWireL -> midWireR;
  botWireL -> botWireR;
  trailingCtrl -> ctrl2 -> ctrl1 -> gate;
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matr);


/** Applies an arbitrarily-controlled two-qubit dense unitary @p matrix to qubits 
 * @p target1 and @p target2 (treated as increasing significance) of @p qureg,
 * conditioned upon the @p controls being in the given @p states.
 * 
 * @diagram
 * @dot
digraph {
  layout=neato;
  rankdir=LR;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  tippytopWireL [shape=plaintext, pos="0,1.5!", label="controls[1]"];
  tippytopWireR [shape=plaintext, pos="2.5,1.5!", label=""];

  topWireL [shape=plaintext, pos="0,1!", label="controls[0]"];
  topWireR [shape=plaintext, pos="2.5,1!", label=""];

  midWireL [shape=plaintext, pos="0,0.5!", label="target2"];
  midWireR [shape=plaintext, pos="2.5,0.5!", label=""];

  botWireL [shape=plaintext, pos="0,0!", label="target1"];
  botWireR [shape=plaintext, pos="2.5,0!", label=""];

  gate [shape=rectangle, label="matrix", style=filled, fillcolor=white, height=1, pos="1.25,0.25!"];
  ctrl1 [shape=circle, label="", width=.12, style=filled, fillcolor=white, pos="1.25,1!"];
  ctrl2 [shape=circle, label="", width=.12, style=filled, fillcolor=black, pos="1.25,1.5!"];
  trailingCtrl [shape=plaintext, label="...", pos="1.25,2!"];

  tippytopWireL -> ctrl2 -> tippytopWireR;
  topWireL -> ctrl1 -> topWireR;
  midWireL -> midWireR;
  botWireL -> botWireR;
  trailingCtrl -> ctrl2 -> ctrl1 -> gate;
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matr);


/** @} */



/** 
 * @defgroup op_compmatr CompMatr
 * @brief Functions for applying general many-target dense matrices, as CompMatr.
 * @{
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


/** @} */



/** 
 * @defgroup op_diagmatr1 DiagMatr1
 * @brief Functions for applying general one-qubit diagonal matrices, as DiagMatr1.
 * @{
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


/** @} */



/** 
 * @defgroup op_diagmatr2 DiagMatr2
 * @brief Functions for applying general two-qubit diagonal matrices, as DiagMatr2.
 * @{
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


/** @} */



/** 
 * @defgroup op_diagmatr DiagMatr
 * @brief Functions for applying general many-qubit diagonal matrices, as DiagMatr.
 * @{
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


/** @} */



/** 
 * @defgroup op_fullstatediagmatr FullStateDiagMatr
 * @brief Functions for applying general all-qubit diagonal matrices, as FullStateDiagMatr.
 * @{
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


/** @} */



/** 
 * @defgroup op_fixed Fixed
 * @brief Functions for applying the one-qubit S, T and Hadamard gates.
 * @{
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


/** @} */



/** 
 * @defgroup op_swap Swap
 * @brief Functions for applying the two-qubit SWAP and related gates.
 * @{
 */


/// @notdoced
/// @notvalidated
void multiplySwap(Qureg qureg, int qubit1, int qubit2);


/** Applies a SWAP gate between @p qubit1 and @p qubit2 of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  layout=neato;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  topWireL [shape=plaintext, label="qubit2", pos="0,.5!"];
  topWireM [shape=point, label="", width=0, pos=".75,.5!"];
  topWireR [shape=plaintext, label="", pos="1.5,.5!"];

  botWireL [shape=plaintext, label="qubit1", pos="0,0!"];
  botWireM [shape=point, label="", width=0, pos=".75,0!"];
  botWireR [shape=plaintext, label="", pos="1.5,0!"];

  topWireL -> topWireR;
  botWireL -> botWireR;
  botWireM -> topWireM;

  topX [shape=plaintext, label="✕", pos=".75,.5!", fontsize=15];
  botX [shape=plaintext, label="✕", pos=".75,0!",  fontsize=15];
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
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


/** @} */



/** 
 * @defgroup op_pauli Pauli
 * @brief Functions for applying the individual one-qubit Pauli operators.
 * @{
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


/** @} */



/** 
 * @defgroup op_paulistr PauliStr
 * @brief Functions for applying a tensor product of Pauli operators, as a PauliStr
 * @{
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


/** @} */



/** 
 * @defgroup op_rotation Rotations
 * @brief Functions for applying one-qubit rotations around Pauli and arbitrary axis.
 * @{
 */


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


/** @} */



/** 
 * @defgroup op_pauligadget Pauli gadgets
 * @brief Functions for applying many-qubit rotations around arbitrary PauliStr.
 * @{
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


/** @} */



/** 
 * @defgroup op_phasegadget Phase gates
 * @brief Functions for applying many-qubit rotations around Pauli Z axis, and phase flips and shifts.
 * @{
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


/// @notdoced
/// @notvalidated
void applyPhaseFlip (Qureg qureg, int target);


/// @notdoced
/// @notvalidated
void applyPhaseShift(Qureg qureg, int target, qreal angle);


/** Applies a two-qubit phase flip upon @p qubit1 and @p qubit2 of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  layout=neato;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  topWireL [shape=plaintext, label="target1", pos="0,.5!"];
  topWireM [shape=point, label="", width=.1, pos=".75,.5!"]
  topWireR [shape=plaintext, label="", pos="1.5,.5!"];

  botWireL [shape=plaintext, label="target2", pos="0,0!"];
  botWireM [shape=point, label="", width=.1, pos=".75,0!"];
  botWireR [shape=plaintext, label="", pos="1.5,0!"];

  topWireL -> topWireR;
  botWireL -> botWireR;
  botWireM -> topWireM;
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyTwoQubitPhaseFlip( Qureg qureg, int target1, int target2);


/** Applies a two-qubit phase flip upon @p qubit1 and @p qubit2 of @p qureg.
 * 
 * @diagram
 * @dot
digraph {
  rankdir=LR;
  layout=neato;
  node [fontsize=10, fontname="Menlo"];
  edge [dir=none];

  topWireL [shape=plaintext, label="target1", pos="0,.5!"];
  topWireM [shape=point, label="", width=.1, pos=".75,.5!"]
  topWireR [shape=plaintext, label="", pos="1.5,.5!"];

  botWireL [shape=plaintext, label="target2", pos="0,0!"];
  botWireM [shape=point, label="", width=.1, pos=".75,0!"];
  botWireR [shape=plaintext, label="", pos="1.5,0!"];

  topWireL -> topWireR;
  botWireL -> botWireR;
  botWireM -> topWireM;

  angle [shape=plaintext, label="θ", pos=".85,-.2!"];
}
 * @enddot
 *
 * @notdoced
 * @notvalidated
 */
void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle);


/// @notdoced
/// @notvalidated
void applyMultiQubitPhaseFlip (Qureg qureg, int* targets, int numTargets);


/// @notdoced
/// @notvalidated
void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle);


/** @} */



/** 
 * @defgroup op_paulistrsum PauliStrSum
 * @brief Functions for applying, exponentiating or Trotterising a weigthed sum of Pauli tensors.
 * @{
 */


/// @notdoced
/// @notvalidated
void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);


/// @notdoced
/// @nottested
void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);


/** @} */



/** 
 * @defgroup op_nots Many-not gates
 * @brief Functions for effecting many-qubit NOT gates
 * @{
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


/** @} */



/** 
 * @defgroup op_measurement Measurements
 * @brief Functions for effecting destructive measurements.
 * @{
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
qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);


/// @notdoced
/// @notvalidated
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);


/// @notdoced
/// @notvalidated
qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);


/** @} */



/** 
 * @defgroup op_projectors Projectors
 * @brief Functions for effecting projectors which break the state normalisation.
 * @{
 */


/// @notdoced
/// @notvalidated
void applyQubitProjector(Qureg qureg, int target, int outcome);


/// @notdoced
/// @notvalidated
void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);


/** @} */



/** 
 * @defgroup op_qft QFT
 * @brief Functions for applying the Quantum Fourier Transform.
 * @{
 */


/// @notdoced
/// @notvalidated
void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets);


/// @notdoced
/// @notvalidated
void applyFullQuantumFourierTransform(Qureg qureg);


/** @} */



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // OPERATIONS_H

/** @} (end doxygen defgroup) */
