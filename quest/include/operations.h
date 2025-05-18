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

#ifdef __cplusplus
    #include <vector>
#endif



/*
 * unlike some other headers, we here intermix the C and C++-only
 * signatures, grouping them semantically & by their doc groups
 */



/** 
 * @defgroup op_compmatr1 CompMatr1
 * @brief Functions for applying general one-qubit dense matrices, as CompMatr1.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


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
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p target is an invalid qubit index.
 * @see
 * - getCompMatr1()
 * - getInlineCompMatr1()
 * - applyCompMatr1()
 * - applyQubitProjector()
 * - multiplyCompMatr()
 * @author Tyson Jones
 */
void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix);


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
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix is not approximately unitary.
 * - if @p target is an invalid qubit index.
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


/** @notyetdoced
 * 
 * Applies a singly-controlled one-qubit dense unitary @p matrix to the specified 
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
 * @see
 * - applyCompMatr1()
 */
void applyControlledCompMatr1(Qureg qureg, int control, int target, CompMatr1 matrix);


/** @notyetdoced
 * 
 * Applies a multiply-controlled one-qubit dense unitary @p matrix to the specified 
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
 * @see
 * - applyCompMatr1()
 */
void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matrix);


/** @notyetdoced
 * 
 * Applies an arbitrarily-controlled one-qubit dense unitary @p matrix to the specified 
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
 * @see
 * - applyCompMatr1()
 */
void applyMultiStateControlledCompMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, CompMatr1 matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledCompMatr1()
void applyMultiControlledCompMatr1(Qureg qureg, std::vector<int> controls, int target, CompMatr1 matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledCompMatr1(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, CompMatr1 matrix);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_compmatr2 CompMatr2
 * @brief Functions for applying general two-qubit dense matrices, as CompMatr2.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see
/// - applyCompMatr2()
/// - multiplyCompMatr1()
void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);


/** @notyetdoced
 * 
 * Applies a general two-qubit dense unitary @p matrix to qubits @p target1 and
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
 * @see
 * - applyCompMatr1()
 */
void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix);


/** @notyetdoced
 * 
 * Applies a singly-controlled two-qubit dense unitary @p matrix to qubits 
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
 * @see 
 * - applyCompMatr2()
 */
void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matr);


/** @notyetdoced
 * 
 * Applies a multiply-controlled two-qubit dense unitary @p matrix to qubits 
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
 * @see
 * - applyCompMatr2()
 */
void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matr);


/** @notyetdoced
 * 
 * Applies an arbitrarily-controlled two-qubit dense unitary @p matrix to qubits 
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
 * @see
 * - applyCompMatr2()
 * - applyMultiStateControlledCompMatr1()
 */
void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matr);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledCompMatr2()
void applyMultiControlledCompMatr2(Qureg qureg, std::vector<int> controls, int target1, int target2, CompMatr2 matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledCompMatr2()
void applyMultiStateControlledCompMatr2(Qureg qureg, std::vector<int> controls, std::vector<int> states, int numControls, int target1, int target2, CompMatr2 matr);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_compmatr CompMatr
 * @brief Functions for applying general many-target dense matrices, as CompMatr.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 * 
 * @see
 * - applyCompMatr()
 * - multiplyCompMatr1()
 */
void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix);


/** @notyetdoced
 * 
 * @formulae
 * 
 * Let @f$ M = @f$ @p matrix.
 * The qubits within @p targets are treated to be ordered least to most significant with respect
 * to @f$ M @f$. That is, if @f$ M @f$ was hypothetically separable single-qubit matrices
 * @f[
      M \equiv A \otimes B \otimes C \otimes \dots 
 * @f]
 * then this function would effect
 * @f[
      \hat{M}_{\text{targets}} \equiv A_{\text{targets}[0]} B_{\text{targets}[1]} C_{\text{targets}[2]} \dots
 * @f]
 *
 * @see
 * - applyCompMatr1()
 */
void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matr);


/// @notyetdoced
/// @see
/// - applyControlledCompMatr1()
void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matr);


/// @notyetdoced
/// @see
/// - applyMultiControlledCompMatr1()
void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matr);


/// @notyetdoced
/// @see
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matr);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see multiplyCompMatr()
void multiplyCompMatr(Qureg qureg, std::vector<int> targets, CompMatr matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyCompMatr()
void applyCompMatr(Qureg qureg, std::vector<int> targets, CompMatr matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyControlledCompMatr()
void applyControlledCompMatr(Qureg qureg, int control, std::vector<int> targets, CompMatr matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledCompMatr()
void applyMultiControlledCompMatr(Qureg qureg, std::vector<int> controls, std::vector<int> targets, CompMatr matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledCompMatr()
void applyMultiStateControlledCompMatr(Qureg qureg, std::vector<int> controls, std::vector<int> states, std::vector<int> targets, CompMatr matr);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_diagmatr1 DiagMatr1
 * @brief Functions for applying general one-qubit diagonal matrices, as DiagMatr1.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);


/// @notyetdoced
/// @see applyCompMatr1()
void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);


/// @notyetdoced
/// @see applyControlledCompMatr1()
void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matr);


/// @notyetdoced
/// @see applyMultiControlledCompMatr1()
void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matr);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matr);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledDiagMatr1()
void applyMultiControlledDiagMatr1(Qureg qureg, std::vector<int> controls, int target, DiagMatr1 matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledDiagMatr1()
void applyMultiStateControlledDiagMatr1(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, DiagMatr1 matr);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_diagmatr2 DiagMatr2
 * @brief Functions for applying general two-qubit diagonal matrices, as DiagMatr2.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);


/// @notyetdoced
/// @see applyCompMatr1()
void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);


/// @notyetdoced
/// @see applyControlledCompMatr1()
void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matr);


/// @notyetdoced
/// @see applyMultiControlledCompMatr1()
void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matr);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matr);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledDiagMatr2()
void applyMultiControlledDiagMatr2(Qureg qureg, std::vector<int> controls, int target1, int target2, DiagMatr2 matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledDiagMatr2()
void applyMultiStateControlledDiagMatr2(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target1, int target2, DiagMatr2 matr);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_diagmatr DiagMatr
 * @brief Functions for applying general many-qubit diagonal matrices, as DiagMatr.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see applyCompMatr1()
void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see applyControlledCompMatr1()
void applyControlledDiagMatr(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see applyMultiControlledCompMatr1()
void applyMultiControlledDiagMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledDiagMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see
/// - multiplyCompMatr1()
/// - applyDiagMatrPower()
void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/** @notyetdoced
 *
 * @formulae
 * 
 * This function is equivalent to applyDiagMatr() except that @p matrix is raised to the given @p exponent.
 */
void applyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @see
/// - applyDiagMatrPower()
/// - applyControlledCompMatr1()
void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @see
/// - applyDiagMatrPower()
/// - applyMultiControlledCompMatr1()
void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @see
/// - applyDiagMatrPower()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledDiagMatrPower(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see multiplyDiagMatr()
void multiplyDiagMatr(Qureg qureg, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyDiagMatr()
void applyDiagMatr(Qureg qureg, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyControlledDiagMatr()
void applyControlledDiagMatr(Qureg qureg, int control, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledDiagMatr()
void applyMultiControlledDiagMatr(Qureg qureg, std::vector<int> controls, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledDiagMatr()
void applyMultiStateControlledDiagMatr(Qureg qureg, std::vector<int> controls, std::vector<int> states, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see multiplyDiagMatrPower()
void multiplyDiagMatrPower(Qureg qureg, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyDiagMatrPower()
void applyDiagMatrPower(Qureg qureg, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyControlledDiagMatrPower()
void applyControlledDiagMatrPower(Qureg qureg, int control, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledDiagMatrPower()
void applyMultiControlledDiagMatrPower(Qureg qureg, std::vector<int> controls, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledDiagMatrPower()
void applyMultiStateControlledDiagMatrPower(Qureg qureg, std::vector<int> controls, std::vector<int> states, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_fullstatediagmatr FullStateDiagMatr
 * @brief Functions for applying general all-qubit diagonal matrices, as FullStateDiagMatr.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - multiplyCompMatr1
void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - multiplyCompMatr1
/// - applyDiagMatrPower
void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @notyetvalidated
void applyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - applyDiagMatrPower
void applyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


// end de-mangler
#ifdef __cplusplus
}
#endif


/** @} */



/** 
 * @defgroup op_fixed Fixed
 * @brief Functions for applying the one-qubit S, T and Hadamard gates.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
void applyS(Qureg qureg, int target);


/// @notyetdoced
void applyControlledS(Qureg qureg, int control, int target);


/// @notyetdoced
void applyMultiControlledS(Qureg qureg, int* controls, int numControls, int target);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledS(Qureg qureg, int* controls, int* states, int numControls, int target);


/// @notyetdoced
void applyT(Qureg qureg, int target);


/// @notyetdoced
void applyControlledT(Qureg qureg, int control, int target);


/// @notyetdoced
void applyMultiControlledT(Qureg qureg, int* controls, int numControls, int target);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledT(Qureg qureg, int* controls, int* states, int numControls, int target);


/// @notyetdoced
void applyHadamard(Qureg qureg, int target);


/// @notyetdoced
void applyControlledHadamard(Qureg qureg, int control, int target);


/// @notyetdoced
void applyMultiControlledHadamard(Qureg qureg, int* controls, int numControls, int target);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledHadamard(Qureg qureg, int* controls, int* states, int numControls, int target);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledS()
void applyMultiControlledS(Qureg qureg, std::vector<int> controls, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledS()
void applyMultiStateControlledS(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledT()
void applyMultiControlledT(Qureg qureg, std::vector<int> controls, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledT()
void applyMultiStateControlledT(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledHadamard()
void applyMultiControlledHadamard(Qureg qureg, std::vector<int> controls, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledHadamard()
void applyMultiStateControlledHadamard(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_swap Swap
 * @brief Functions for applying the two-qubit SWAP and related gates.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
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
 * @notyetdoced
 */
void applySwap(Qureg qureg, int qubit1, int qubit2);


/// @notyetdoced
void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2);


/// @notyetdoced
void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);


/// @notyetdoced
void applySqrtSwap(Qureg qureg, int qubit1, int qubit2);


/// @notyetdoced
void applyControlledSqrtSwap(Qureg qureg, int control, int qubit1, int qubit2);


/// @notyetdoced
void applyMultiControlledSqrtSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledSqrtSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledSwap()
void applyMultiControlledSwap(Qureg qureg, std::vector<int> controls, int qubit1, int qubit2);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledSwap()
void applyMultiStateControlledSwap(Qureg qureg, std::vector<int> controls, std::vector<int> states, int qubit1, int qubit2);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledSqrtSwap()
void applyMultiControlledSqrtSwap(Qureg qureg, std::vector<int> controls, int qubit1, int qubit2);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledSqrtSwap()
void applyMultiStateControlledSqrtSwap(Qureg qureg, std::vector<int> controls, std::vector<int> states, int numControls, int qubit1, int qubit2);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_pauli Pauli
 * @brief Functions for applying the individual one-qubit Pauli operators.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyPauliX(Qureg qureg, int target);


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyPauliY(Qureg qureg, int target);


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyPauliZ(Qureg qureg, int target);


/// @notyetdoced
void applyPauliX(Qureg qureg, int target);


/// @notyetdoced
void applyPauliY(Qureg qureg, int target);


/// @notyetdoced
void applyPauliZ(Qureg qureg, int target);


/// @notyetdoced
void applyControlledPauliX(Qureg qureg, int control, int target);


/// @notyetdoced
void applyControlledPauliY(Qureg qureg, int control, int target);


/// @notyetdoced
void applyControlledPauliZ(Qureg qureg, int control, int target);


/// @notyetdoced
void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target);


/// @notyetdoced
void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target);


/// @notyetdoced
void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledPauliZ(Qureg qureg, int* controls, int* states, int numControls, int target);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledPauliX()
void applyMultiControlledPauliX(Qureg qureg, std::vector<int> controls, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledPauliY()
void applyMultiControlledPauliY(Qureg qureg, std::vector<int> controls, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledPauliZ()
void applyMultiControlledPauliZ(Qureg qureg, std::vector<int> controls, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledPauliX()
void applyMultiStateControlledPauliX(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledPauliY()
void applyMultiStateControlledPauliY(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledPauliZ()
void applyMultiStateControlledPauliZ(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_paulistr PauliStr
 * @brief Functions for applying a tensor product of Pauli operators, as a PauliStr
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyPauliStr(Qureg qureg, PauliStr str);


/// @notyetdoced
void applyPauliStr(Qureg qureg, PauliStr str);


/// @notyetdoced
void applyControlledPauliStr(Qureg qureg, int control, PauliStr str);


/// @notyetdoced
void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str);


/// @notyetdoced
/// @see applyMultiStateControlledCompMatr1()
void applyMultiStateControlledPauliStr(Qureg qureg, int* controls, int* states, int numControls, PauliStr str);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledPauliStr()
void applyMultiControlledPauliStr(Qureg qureg, std::vector<int> controls, PauliStr str);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledPauliStr()
void applyMultiStateControlledPauliStr(Qureg qureg, std::vector<int> controls, std::vector<int> states, PauliStr str);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_rotation Rotations
 * @brief Functions for applying one-qubit rotations around Pauli and arbitrary axis.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 * 
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle.
 * This function effects unitary
 * @f[
      \hat{R}_{x}(\theta) 
        = 
        \exp \left(
          - \iu \frac{\theta}{2} 
            \hat{\sigma}_x
        \right)
 * @f]
 * upon the @p target qubit, where @f$ \hat{\sigma}_x @f$ is the Pauli X matrix.
 *
 * @equivalences
 * - This function is entirely equivalent to calling applyPauliGadget() with a single-site PauliStr.
 *   ```
     applyPauliGadget(qureg, getInlinePauliStr("X", {target}), angle);
 *   ```
 * - This function is faster than, but otherwise equivalent to, invoking applyRotateAroundAxis()
 *   with an axis vector equal to the X-axis.
 *   ```
     applyRotateAroundAxis(qureg, target, qreal angle, 1,0,0);
 *   ```
 * - This function is faster than, but otherwise equivalent to, effecting @f$ \hat{R}_{x}(\theta) @f$ as a CompMatr1.
 *   ```
     qcomp c = cos(angle/2);
     qcomp s = sin(angle/2) * (-1.i);
     CompMatr1 matr = getInlineCompMatr1({{c, s}, {s, c}});
     applyCompMatr1(qureg, target, matr);
 *   ```
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyRotateX(Qureg qureg, int target, qreal angle);


/** @notyetdoced
 * 
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle.
 * This function effects unitary
 * @f[
      \hat{R}_{y}(\theta) 
        = 
        \exp \left(
          - \iu \frac{\theta}{2} 
            \hat{\sigma}_y
        \right)
 * @f]
 * upon the @p target qubit, where @f$ \hat{\sigma}_y @f$ is the Pauli Y matrix.
 *
 * @equivalences
 * - This function is entirely equivalent to calling applyPauliGadget() with a single-site PauliStr.
 *   ```
     applyPauliGadget(qureg, getInlinePauliStr("Y", {target}), angle);
 *   ```
 * - This function is faster than, but otherwise equivalent to, invoking applyRotateAroundAxis()
 *   with an axis vector equal to the Y-axis.
 *   ```
     applyRotateAroundAxis(qureg, target, qreal angle, 0,1,0);
 *   ```
 * - This function is faster than, but otherwise equivalent to, effecting @f$ \hat{R}_{y}(\theta) @f$ as a CompMatr1.
 *   ```
     qcomp c = cos(angle/2);
     qcomp s = sin(angle/2);
     CompMatr1 matr = getInlineCompMatr1({{c, -s}, {s, c}});
     applyCompMatr1(qureg, target, matr);
 *   ```
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyRotateY(Qureg qureg, int target, qreal angle);


/** @notyetdoced
 * 
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle.
 * This function effects unitary
 * @f[
      \hat{R}_{z}(\theta) 
        = 
        \exp \left(
          - \iu \frac{\theta}{2} 
            \hat{\sigma}_z
        \right)
 * @f]
 * upon the @p target qubit, where @f$ \hat{\sigma}_z @f$ is the Pauli Z matrix.
 *
 * @equivalences
 * - This function is entirely equivalent to calling applyPauliGadget() with a single-site PauliStr.
 *   ```
     applyPauliGadget(qureg, getInlinePauliStr("Z", {target}), angle);
 *   ```
 * - This function is faster than, but otherwise equivalent to, invoking applyRotateAroundAxis()
 *   with an axis vector equal to the Z-axis.
 *   ```
     applyRotateAroundAxis(qureg, target, qreal angle, 0,0,1);
 *   ```
 * - This function is faster than, but otherwise equivalent to, effecting @f$ \hat{R}_{z}(\theta) @f$ as a DiagMatr1.
 *   ```
     qcomp a = cexp(- angle / 2 * 1.i);
     qcomp b = cexp(  angle / 2 * 1.i);
     DiagMatr1 matr = getInlineDiagMatr1({a, b});
     applyDiagMatr1(qureg, target, matr);
 *   ```
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyRotateZ(Qureg qureg, int target, qreal angle);


/// @notyetdoced
void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle);


/// @notyetdoced
void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle);


/// @notyetdoced
void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle);


/// @notyetdoced
void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle);


/// @notyetdoced
void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle);


/// @notyetdoced
void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle);


/// @notyetdoced
/// @see 
/// - applyRotateX()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);


/// @notyetdoced
/// @see 
/// - applyRotateY()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);


/// @notyetdoced
/// @see 
/// - applyRotateZ()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);


/** @notyetdoced
 *
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle and  @f$ \vec{n} = ( @f$ @p axisX, @p axisY, @p axisZ @f$ ) @f$,
 * with corresponding unit vector @f$ \bar{n} @f$. 
 * Further, let @f$ \vec{\sigma} = (\hat{\sigma}_x, \hat{\sigma}_y, \hat{\sigma}_z) @f$ denote a vector of the Pauli matrices.
 * 
 * This function effects unitary
 * @f[
      \hat{R}_{\bar{n}}(\theta) 
        = 
        \exp \left(
          - \iu \frac{\theta}{2} 
            \bar{n} \cdot \vec{\sigma}
        \right)
 * @f]
 * upon the target qubit. Explicitly,
 * @f[
      \hat{R}_{\bar{n}}(\theta) 
        \equiv 
        \begin{pmatrix}
        \cos\left( \frac{\theta}{2} \right) - \iu \, \bar{n}_z \sin\left( \frac{\theta}{2} \right)
          &
        - \, (\bar{n}_y + \bar{n}_x \, \iu ) \sin\left( \frac{\theta}{2} \right)
          \\
        (\bar{n}_y - \bar{n}_x \, \iu ) \sin\left( \frac{\theta}{2} \right)
          &
        \cos\left( \frac{\theta}{2} \right) + \iu \, \bar{n}_z \sin\left( \frac{\theta}{2} \right)
        \end{pmatrix}
 * @f]
 * where 
 * @f[
      \bar{n}_i 
        = 
      \frac{\vec{n}_i}{\| \vec{n} \|_2}
        =
      \frac{\vec{n}_i}{ \sqrt{ {\vec{n}_x}^2 + {\vec{n}_y}^2 + {\vec{n}_z}^2 } }.
 * @f]
 *
 * @equivalences
 * - Assuming @f$ \| \vec{n} \|_2 \ne 0 @f$, this function is agnostic to the normalisation
 *   of the axis vector.
 *   ```
     applyRotateAroundAxis(qureg, target, angle, x,  y,  z);
     applyRotateAroundAxis(qureg, target, angle, 5*x,5*y,5*z); // equivalent
 *   ```
 * - This function is entirely equivalent to preparing @f$ \hat{R}_{\bar{n}}(\theta) @f$
 *   as a CompMatr1 and effecting it upon the state via applyCompMatr1().
 * - This function is both more accurate and efficient than equivalently instantiating a 
 *   three-term PauliStrSum @f$ \hat{H} = \bar{n} \cdot \vec{\sigma}@f$ and effecting
 *   @f$ \exp \left(\iu \alpha \hat{H} \right) @f$ via applyTrotterizedPauliStrSumGadget() 
 *   with @f$ \alpha = - \theta/2 @f$ and very many repetitions.
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyRotateAroundAxis(Qureg qureg, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/// @notyetdoced
void applyControlledRotateAroundAxis(Qureg qureg, int ctrl, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/// @notyetdoced
void applyMultiControlledRotateAroundAxis(Qureg qureg, int* ctrls, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/// @notyetdoced
/// @see 
/// - applyRotateAroundAxis()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledRotateAroundAxis(Qureg qureg, int* ctrls, int* states, int numCtrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledRotateX()
void applyMultiControlledRotateX(Qureg qureg, std::vector<int> controls, int target, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledRotateY()
void applyMultiControlledRotateY(Qureg qureg, std::vector<int> controls, int target, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledRotateZ()
void applyMultiControlledRotateZ(Qureg qureg, std::vector<int> controls, int target, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledRotateX()
void applyMultiStateControlledRotateX(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledRotateY()
void applyMultiStateControlledRotateY(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledRotateZ()
void applyMultiStateControlledRotateZ(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledRotateAroundAxis()
void applyMultiControlledRotateAroundAxis(Qureg qureg, std::vector<int> ctrls, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledRotateAroundAxis()
void applyMultiStateControlledRotateAroundAxis(Qureg qureg, std::vector<int> ctrls, std::vector<int> states, int targ, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_pauligadget Pauli gadgets
 * @brief Functions for applying many-qubit rotations around arbitrary PauliStr.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - multiplyCompMatr1()
/// - applyPauliGadget()
void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \hat{\sigma} = @f$ @p str and @f$ \theta = @f$ @p angle. 
 * 
 * This function effects unitary
 * @f[
      R_{\hat{\sigma}}(\theta) = \exp \left( - \iu \, \frac{\theta}{2} \, \hat{\sigma} \right),
 * @f]
 * which affects only the qubits for which @f$ \hat{\sigma} @f$ is not the identity
 * Pauli. As such, this effects a multi-qubit rotation around an arbitrary Pauli string.
 * 
 * @equivalences
 * - Because @f$ R_{\hat{\sigma}}(\theta) @f$ satisfies
 *   @f[
        R_{\hat{\sigma}}(\theta) \equiv 
          \cos\left( \frac{\theta}{2} \right) \, \id 
          - \iu  \sin\left( \frac{\theta}{2} \right) \, \hat{\sigma},
 *   @f]
 *   this function is equivalent to (but much faster than) effecting @f$ \hat{\sigma} @f$
 *   upon a clone which is subsequently superposed.
 *   ```
     // prepare |temp> = str |qureg>
     Qureg temp = createCloneQureg(qureg);
     applyPauliStr(temp, str);

     // set |qureg> = cos(theta/2) |qureg> - i sin(theta/2) str |qureg>
     setQuregToSuperposition(cos(theta/2), qureg, - 1.0i * sin(theta/2), temp, 0, temp);
 *   ```
 * - When @p str contains only @f$ \hat{Z} @f$ or @f$ \id @f$ Paulis, this function will
 *   automatically invoke applyPhaseGadget() which leverages an optimised implementation.
 * - When @p str contains only @f$ \id @f$ Paulis, this function merely effects a change
 *   of global phase upon statevectors of @f$ -\theta/2 @f$, leaving density matrices
 *   unchanged.
 *   ```
     qcomp factor = cexp(- theta / 2 * 1.i);
     setQuregToSuperposition(factor, qureg, 0,qureg,0,qureg);
 *   ```
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(10);
    qreal theta = 3.14;
    
    // verbosely
    int numPaulis = 4;
    char* paulis = "XYIZ";
    int targets[] = {0,1,5,7};
    PauliStr str = getPauliStr(paulis, targets, numPaulis);
    applyPauliGadget(qureg, str, angle);

    // concisely
    applyPauliGadget(qureg, getInlinePauliStr("XYZ",{0,1,7}), theta);
 * ```
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


/// @notyetdoced
void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle);


/// @notyetdoced
void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle);


/// @notyetdoced
/// @see
/// - applyPauliGadget()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledPauliGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStr str, qreal angle);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledPauliGadget()
void applyMultiControlledPauliGadget(Qureg qureg, std::vector<int> controls, PauliStr str, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledPauliGadget()
void applyMultiStateControlledPauliGadget(Qureg qureg, std::vector<int> controls, std::vector<int> states, PauliStr str, qreal angle);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_phasegadget Phase gates
 * @brief Functions for applying many-qubit rotations around Pauli Z axis, and phase flips and shifts.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - multiplyCompMatr1()
/// - applyPhaseGadget
void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


/** @notyetdoced
 * 
 * @formulae
 * 
 * Let @f$ \vec{t} = @f$ @p targets and @f$ \theta = @f$ @p angle.
 * 
 * This function effects diagonal unitary
 * @f[
      R_{\hat{Z}}(\theta) = \exp \left( - \iu \, \frac{\theta}{2} \, \bigotimes_{t \,\in\, \vec{t}} \hat{Z}_t \right).
 * @f]
 *
 * @equivalences
 * - This function is equivalent to calling applyPauliGadget() with a PauliStr containing only @f$ \hat{Z} @f$ and @f$ \id @f$.
 *   This latter function will actually automatically invoke applyPhaseGadget() which has an optimised implementation.
 * - This function is equivalent to, albeit much faster than, preparing a DiagMatr with @f$ \pm 1 @f$ elements (depending upon
 *   the parity of the targeted set bits) and effecting it with applyDiagMatr().
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


/// @notyetdoced
void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle);


/// @notyetdoced
void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle);


/// @notyetdoced
/// @see
/// - applyPhaseGadget()
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledPhaseGadget(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, qreal angle);


/** @notyetdoced
 * 
 * This function is a mere alias of applyPauliZ(), meaningfully differing only for many targets.
 */
void applyPhaseFlip(Qureg qureg, int target);


/** @notyetdoced
 * 
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle. This function effects diagonal unitary
 * 
 * @f[
      \hat{U}(\theta) = \begin{pmatrix} 1 & 0 \\ 0 & e^{\iu \theta} \end{pmatrix}
 * @f]
 * upon the @p target qubit.
 * 
 * @equivalences
 * - This function is equivalent to, albeit much faster than, a Z-axis rotation with
 *   an adjustment to the global phase (which is redundant upon density matrices).
 *   @f[
 *      \hat{U}(\theta) \equiv \hat{R}_z(\theta) \cdot e^{\iu \frac{\theta}{2}} \hat{\id}
 *   @f]
 *   ```
     applyRotateZ(qureg, target, angle);
     applyPauliGadget(qureg, getPauliStr("I"), angle); // global phase
 *   ```
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyPhaseShift(Qureg qureg, int target, qreal angle);


/** @notyetdoced
 * 
 * Applies a two-qubit phase flip upon qubits @p target1 and @p target2 of @p qureg.
 * 
 * @formulae
 * 
 * This function flips the sign of all computational basis states for which
 * the targeted qubits are in state @f$ \ket{1}\ket{1} @f$. This is equivalent
 * to the diagonal unitary
 * 
 * @f[
      \hat{U}(\theta) = \begin{pmatrix} 1 \\ & 1 \\ & & 1 \\ & & & -1 \end{pmatrix},
 * @f]
 * effected upon the target qubits.
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
 * @equivalences
 * - The target qubits are interchangeable, ergo
 *   ```
     applyTwoQubitPhaseFlip(qureg, target1, target2);
     applyTwoQubitPhaseFlip(qureg, target2, target1); // equivalent
 *   ```
 * - This function is entirely equivalent to a controlled Pauli-Z unitary (or a hypothetical
 *   controlled variant of applyPhaseFlip()) with either target qubit substituted for the control qubit.
 *   ```
     applyControlledPauliZ(qureg, target1, target2);
 *   ```
 * - This function is faster and more accurate than, but otherwise equivalent to, a two-qubit phase shift
 *   with angle @f$ = \pi @f$.
 *   ```
     applyTwoQubitPhaseShift(qureg, target1, target2, 3.141592653); // approx equiv
 *   ```
 */
void applyTwoQubitPhaseFlip(Qureg qureg, int target1, int target2);


/** @notyetdoced
 * 
 * Applies a two-qubit phase shift upon qubits @p target1 and @p target2 of @p qureg.
 * 
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle.
 * This function multiplies factor @f$ e^{\iu \theta} @f$ upon all computational basis states 
 * for which the targeted qubits are in state @f$ \ket{1}\ket{1} @f$. This is equivalent
 * to the diagonal unitary
 * 
 * @f[
      \hat{U}(\theta) = \begin{pmatrix} 1 \\ & 1 \\ & & 1 \\ & & & e^{\iu \theta} \end{pmatrix},
 * @f]
 * effected upon the target qubits.
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
 * @equivalences
 * - The target qubits are interchangeable, ergo
 *   ```
     applyTwoQubitPhaseShift(qureg, target1, target2, angle);
     applyTwoQubitPhaseShift(qureg, target2, target1, angle); // equivalent
 *   ```
 * - This function is equivalent to a controlled variant of applyPhaseShift(), treating
 *   either target qubit as the control qubit.
 * - This function generalises applyTwoQubitPhaseFlip() to arbitrary changes in phase.
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyTwoQubitPhaseShift(Qureg qureg, int target1, int target2, qreal angle);


/** @notyetdoced
 * 
 * @formulae
 * 
 * This function flips the sign of all computational basis states for which
 * the targeted qubits are all in state @f$ \ket{1} @f$. This is equivalent
 * to the diagonal unitary
 * @f[
      \hat{U}(\theta) = \begin{pmatrix} 1 \\  & \ddots \\ & & 1 \\ & & & -1 \end{pmatrix},
 * @f]
 * effected upon the target qubits.
 * 
 * @equivalences
 * - The ordering of @p targets has no affect on the effected operation.
 * - This function is entirely equivalent to a multi-controlled Pauli-Z unitary (or a hypothetical
 *   many-controlled variant of applyPhaseFlip()) with all but one arbitrary target qubit becoming
 *   control qubits.
 *   ```
     applyMultiControlledPauliZ(qureg, targets, numTargets-1, targets[0]);
 *   ```
 * - This function is faster and more accurate than, but otherwise equivalent to, a multi-qubit phase shift
 *   with angle @f$ = \pi @f$.
 *   ```
     applyMultiQubitPhaseShift(qureg, targets, numTargets, 3.141592653); // approx equiv
 *   ```
 */
void applyMultiQubitPhaseFlip(Qureg qureg, int* targets, int numTargets);


/** @notyetdoced
 *
 * @formulae
 * 
 * Let @f$ \theta = @f$ @p angle.
 * This function multiplies factor @f$ e^{\iu \theta} @f$ upon all computational basis states 
 * for which all targeted qubits are in state @f$ \ket{1} @f$. This is equivalent
 * to the diagonal unitary
 * @f[
      \hat{U}(\theta) = \begin{pmatrix} 1 \\  & \ddots \\ & & 1 \\ & & & e^{\iu \theta} \end{pmatrix},
 * @f]
 * effected upon the target qubits.
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
 * @equivalences
 * - The ordering of @p targets has no affect on the effected operation.
 * - This function is equivalent to a multi-controlled variant of applyPhaseShift(), treating all
 *   but one arbitrary target qubit as control qubits.
 * - This function generalises applyMultiQubitPhaseFlip() to arbitrary changes in phase.
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyMultiQubitPhaseShift(Qureg qureg, int* targets, int numTargets, qreal angle);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see multiplyPhaseGadget()
void multiplyPhaseGadget(Qureg qureg, std::vector<int> targets, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyPhaseGadget()
void applyPhaseGadget(Qureg qureg, std::vector<int> targets, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyControlledPhaseGadget()
void applyControlledPhaseGadget(Qureg qureg, int control, std::vector<int> targets, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledPhaseGadget()
void applyMultiControlledPhaseGadget(Qureg qureg, std::vector<int> controls, std::vector<int> targets, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledPhaseGadget()
void applyMultiStateControlledPhaseGadget(Qureg qureg, std::vector<int> controls, std::vector<int> states, std::vector<int> targets, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitPhaseFlip()
void applyMultiQubitPhaseFlip(Qureg qureg, std::vector<int> targets);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitPhaseShift()
void applyMultiQubitPhaseShift(Qureg qureg, std::vector<int> targets, qreal angle);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_paulistrsum PauliStrSum
 * @brief Functions for applying, exponentiating or Trotterising a weigthed sum of Pauli tensors.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
/// @see multiplyCompMatr1()
void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);


/** @notyetdoced
 * @notyettested
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p sum and @f$ \theta = @f$ @p angle. This function approximates the action of
 * @f[
      \exp \left(\iu \, \theta \, \hat{H} \right)
 * @f]
 * via a Trotter-Suzuki decomposition of the specified @p order and number of repetitions (@p reps).
 * 
 * 
 * To be precise, let @f$ r = @f$ @p reps and assume @p sum is composed of
 * @f$ T @f$-many terms of the form
 * @f[
      \hat{H} = \sum\limits_j^T c_j \, \hat{\sigma}_j
 * @f]
 * where @f$ c_j @f$ is the (necessarily real) coefficient of the @f$ j @f$-th PauliStr @f$ \hat{\sigma}_j @f$.
 * 
 * - When @p order=1, this function performs first-order Trotterisation, whereby
 *   @f[
       \exp(\iu \, \theta \, \hat{H} )
          \approx 
        \prod\limits^{r} 
        \prod\limits_{j=1}^{T} 
        \exp \left( \iu \, \frac{\theta \, c_j}{r} \, \hat\sigma_j \right).
 *   @f]
 * - When @p order=2, this function performs the lowest order "symmetrized" Suzuki decomposition, whereby 
 *   @f[
       \exp(\iu \, \theta \, \hat{H} )
          \approx 
        \prod\limits^{r} \left[
             \prod\limits_{j=1}^{T} \exp \left( \iu \frac{\theta \, c_j}{2 \, r}  \hat\sigma_j \right)
              \prod\limits_{j=T}^{1} \exp \left( \iu \frac{\theta \, c_j}{2 \, r}  \hat\sigma_j \right)
         \right].
 *   @f]
 * - Greater, even values of @p order (denoted by symbol @f$ n @f$) invoke higher-order symmetrized decompositions 
 *   @f$ S[\theta,n,r] @f$. Letting @f$ p = \left( 4 - 4^{1/(n-1)} \right)^{-1} @f$, these satisfy
 *   @f{align*}
        S[\theta, n, 1] &= 
            \left( \prod\limits^2 S[p \, \theta, n-2, 1] \right)
            S[ (1-4p)\,\theta, n-2, 1]
            \left( \prod\limits^2 S[p \, \theta, n-2, 1] \right),
        \\
        S[\theta, n, r] &= 
            \prod\limits^{r} S\left[\frac{\theta}{r}, n, 1\right].
 *   @f}
 * 
 * > These formulations are taken from 'Finding Exponential Product Formulas
 * > of Higher Orders', Naomichi Hatano and Masuo Suzuki (2005) (<a href="https://arxiv.org/abs/math-ph/0506007">arXiv</a>).
 */
void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);


// end de-mangler
#ifdef __cplusplus
}
#endif


/** @} */



/** 
 * @defgroup op_nots Many-not gates
 * @brief Functions for effecting many-qubit NOT gates
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see multiplyCompMatr1()
void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


/// @notyetdoced
void applyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


/// @notyetdoced
void applyControlledMultiQubitNot(Qureg qureg, int control, int* targets, int numTargets);


/// @notyetdoced
void applyMultiControlledMultiQubitNot(Qureg qureg, int* controls, int numControls, int* targets, int numTargets);


/// @notyetdoced
/// @see
/// - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledMultiQubitNot(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see multiplyMultiQubitNot()
void multiplyMultiQubitNot(Qureg qureg, std::vector<int> targets);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitNot()
void applyMultiQubitNot(Qureg qureg, std::vector<int> targets);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyControlledMultiQubitNot()
void applyControlledMultiQubitNot(Qureg qureg, int control, std::vector<int> targets);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledMultiQubitNot()
void applyMultiControlledMultiQubitNot(Qureg qureg, std::vector<int> controls, std::vector<int> targets);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledMultiQubitNot()
void applyMultiStateControlledMultiQubitNot(Qureg qureg, std::vector<int> controls, std::vector<int> states, std::vector<int> targets);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_measurement Measurements
 * @brief Functions for effecting destructive measurements.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
int applyQubitMeasurement(Qureg qureg, int target);


/// @notyetdoced
/// @notyetvalidated
int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability);


/// @notyetdoced
/// @notyetvalidated
qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome);


/// @notyetdoced
/// @notyetvalidated
qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);


/// @notyetdoced
/// @notyetvalidated
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);


/// @notyetdoced
/// @notyetvalidated
qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitMeasurementAndGetProb()
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, std::vector<int> qubits, qreal* probability);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyForcedMultiQubitMeasurement()
qreal applyForcedMultiQubitMeasurement(Qureg qureg, std::vector<int> qubits, std::vector<int> outcomes);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_projectors Projectors
 * @brief Functions for effecting projectors which break the state normalisation.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
void applyQubitProjector(Qureg qureg, int target, int outcome);


/// @notyetdoced
/// @notyetvalidated
void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitProjector()
void applyMultiQubitProjector(Qureg qureg, std::vector<int> qubits, std::vector<int> outcomes);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_qft QFT
 * @brief Functions for applying the Quantum Fourier Transform.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets);


/// @notyetdoced
/// @notyetvalidated
void applyFullQuantumFourierTransform(Qureg qureg);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyQuantumFourierTransform()
void applyQuantumFourierTransform(Qureg qureg, std::vector<int> targets);


#endif // __cplusplus

/** @} */



#endif // OPERATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
