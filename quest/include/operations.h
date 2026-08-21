/** @file
 * API signatures for effecting mostly physical and/or trace
 * preserving operators, such as unitaries, gates and 
 * measurements, upon Quregs which are instantiated as both 
 * statevectors or density matrices. This excludes Trotterised
 * gadgets and evolutions (exposed instead in trotterisation.h),
 * functions to pre- or post-multiply operators upon density
 * matrices (multiplication.h) and decoherence channels
 * (decoherence.h).
 * 
 * @author Tyson Jones
 * @author Diogo Pratas Maia (non-unitary Pauli gadget)
 * 
 * @defgroup operations Operations
 * @ingroup api
 * @brief Functions for effecting standard operators upon Quregs.
 * @{
 */

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"

#include <stdbool.h>

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
 * 
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
 * 
 * - Unitarity of @f$ \hat{U} = @f$ @p matrix requires that 
 *   @f$ \hat{U} \hat{U}^\dagger = \id @f$. Validation will check that @p matrix is
 *   approximately unitary via
 *   @f[ 
        \max\limits_{ij} \Big|\left(\hat{U} \hat{U}^\dagger - \id\right)_{ij}\Big|^2 \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setQuESTValidationEpsilon().
 * 
 * @myexample
 * 
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
 * - leftapplyCompMatr1()
 * - rightapplyCompMatr1()
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
 * @formulae
 * 
 * Let @f$ \hat{U} = @f$ @p matrix, @f$ t = @f$ @p target, @f$ c = @f$ @p control,
 * and let @f$\hat{O}_q@f$ denote an operator upon the @f$q@f$-th qubit.
 * This function effects operator
 * @f[
    C_c[\hat{U}_t] = \ketbra{0}{0}_c \otimes \id_t + \ketbra{1}{1}_c \otimes \hat{U}_t,
 * @f]
 * where @f$\hat{U}@f$ is effected upon basis states for which qubit @f$c@f$ has value `1`.
 * For illustration, when @p control=0 and @p target=1, this function would effect
 * @f[
    C_1[\hat{U}_0] \equiv 
    \begin{pmatrix} 
      1 \\ & 1 \\ & & u_{00} & u_{01} \\ & & u_{10} & u_{11}
    \end{pmatrix}.
 * @f]
 *
 * This operation can be performed upon statevectors and density matrices.
 *
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function effects
 *   @f[ 
        \svpsi \rightarrow C_c[\hat{U}_t] \, \svpsi.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function effects
 *   @f[ 
        \dmrho \rightarrow C_c[\hat{U}_t] \, \dmrho \, {C_c[\hat{U}_t]}^\dagger.
 *   @f]
 *
 * The amplitudes which _are_ modified, are done so in an identical fashion as in applyCompMatr1().
 *
 * @constraints
 * 
 * - Unitarity of @f$ \hat{U} = @f$ @p matrix requires that 
 *   @f$ \hat{U} \hat{U}^\dagger = \id @f$. Validation will check that @p matrix is
 *   approximately unitary via
 *   @f[ 
        \max\limits_{ij} \Big|\left(\hat{U} \hat{U}^\dagger - \id\right)_{ij}\Big|^2 \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setQuESTValidationEpsilon().
 *
 * @equivalences
 * 
 * - This function is faster than, but mathematically equivalent to, initialising a two-qubit
 *   matrix (CompMatr2) to the @f$C_1[\hat{U}_0]@f$ matrix above, and calling applyCompMatr2(),
 *   passing @p control as the most significant target.
 * ```
     CompMatr2 m = getInlineCompMatr2({
         {1,0,0,0}, 
         {0,1,0,0}, 
         {0,0,u00,u01}, 
         {0,0,u10,u11}});
     
     applyCompMatr2(qureg, target, control, m);
 * ```
 *
 * @myexample
 * 
 * ```
    Qureg qureg = createQureg(5);

    CompMatr1 matrix = getInlineCompMatr1({
        {-1i/sqrt(2), 1i/sqrt(2)},
        {(1i-1)/2,    (1i-1)/2}
    });

    // C_0[U_2]
    applyControlledCompMatr1(qureg, 0, 2, matrix); 
 * ```

 * @param[in,out] qureg   the state to modify.
 * @param[in]     control the index of the control qubit.
 * @param[in]     target  the index of the target qubit.
 * @param[in]     matrix  the Z-basis unitary matrix to effect.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix is not approximately unitary.
 * - if @p control or @p target are an invalid qubit index.
 * - if @p control and @p target are equal.
 * @see
 * - applyCompMatr1()
 * - applyMultiControlledCompMatr1()
 * - applyMultiStateControlledCompMatr1()
 * @author Tyson Jones
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
 * @formulae
 * 
 * Let @f$ \vec{c} = @f$ @p controls, @f$ t = @f$ @p target, and @f$ \hat{U} = @f$ @p matrix.
 * This functions effects operator
 * 
 * @f[
    C_{\vec{c}}[\hat{U}_t]
 * @f]
 *
 * which is equivalent to applying @f$ \hat{U}_t @f$ upon only the computational basis states for which 
 * all control qubits are in the @f$ \ket{1} @f$ state.
 *
 * Precisely, let @f$n = 2^{|\vec{c}|}-1@f$. Then
 * @f[
    C_{\vec{c}}[\hat{U}_t] = \sum\limits_{i=0}^{n-1} \ketbra{i}{i}_{\vec{c}} \otimes \hat{\id}_t
      + \ketbra{n}{n}_{\vec{c}} \otimes \hat{U}_t
 * @f]
 *
 * The amplitudes which _are_ modified, are done so in an identical fashion as in applyCompMatr1().
 *
 * @param[in,out] qureg       the state to modify.
 * @param[in]     controls    a list of control qubits.
 * @param[in]     numControls the length of @p controls.
 * @param[in]     target      the target qubit.
 * @param[in]     matrix      the Z-basis unitary matrix to effect.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix is not approximately unitary.
 * - if @p target or any element of @p controls are an invalid qubit index.
 * - if @p controls contains duplicates, or includes @p target.
 * - if @p numControls is negative.
 * @see
 * - applyCompMatr1()
 * - applyMultiStateControlledCompMatr1()
 * @author Tyson Jones
 */
void applyMultiControlledCompMatr1(Qureg qureg, int* controls, int numControls, int target, CompMatr1 matrix);


/** Applies an arbitrarily-controlled one-qubit dense unitary @p matrix to the specified 
 * @p target qubit of @p qureg, conditioned upon the @p controls being in the corresponding @p states.
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
 * @formulae
 * 
 * Let @f$ \vec{c} = @f$ @p controls, @f$ t = @f$ @p target, @f$ \hat{U} = @f$ @p matrix and
 * @f$n = 2^{|\vec{c}|}-1@f$. Let @f$ \ket{s}_{\vec{c}} @f$ be the computational substate formed by
 * the qubits in @p controls being in the corresponding @p states.
 * 
 * This function applies the operator 
 * 
 * @f[
      \sum\limits_{i=0, i \ne s}^{n} \ketbra{i}{i}_{\vec{c}} \otimes \hat{\id}_t
      + \ketbra{s}{s}_{\vec{c}} \otimes \hat{U}_t
 * @f]
 *
 * The amplitudes which _are_ modified, are done so in an identical fashion as in applyCompMatr1().
 * 
 * @equivalences
 * 
 * - This function is faster than, but mathematically equivalent to, applying a Pauli @c X
 *   upon every zero-controlled qubit, applying the matrix with all one-controls, then undoing
 *   the flipped qubits.
 *   ```cpp
      for (int i=0; i<numControls; i++)
          if (states[i] == 0)
              applyPauliX(qureg, controls[i]);

      applyMultiControlledCompMatr1(qureg, controls, numControls, target, matrix);

      for (int i=0; i<numControls; i++)
          if (states[i] == 0)
              applyPauliX(qureg, controls[i]);
 *   ```
 *
 * @param[in,out] qureg       the state to modify.
 * @param[in]     controls    a list of control qubits.
 * @param[in]     states      a list of corresponding qubit states (each, @c 0 or @c 1).
 * @param[in]     numControls the length of @p controls and @p states.
 * @param[in]     target      the target qubit.
 * @param[in]     matrix      the Z-basis unitary matrix to effect.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix is not approximately unitary.
 * - if @p target or any element of @p controls are an invalid qubit index.
 * - if @p controls contains duplicates, or includes @p target.
 * - if @p numControls is negative.
 * - if @p states contains any element besides @c 0 or @c 1.
 * @see
 * - applyCompMatr1()
 * @author Tyson Jones
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
 * - leftapplyCompMatr2()
 * - rightapplyCompMatr2()
 * @author Tyson Jones
 */
void applyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix);


/** @notyetdoced
 * 
 * Applies a singly-controlled two-qubit dense unitary @p matrix to qubits 
 * @p target1 and @p target2 (treated as increasing significance) of @p qureg.
 * 
 * > - See applyCompMatr2() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about the @p control qubit.
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
 * @author Tyson Jones
 */
void applyControlledCompMatr2(Qureg qureg, int control, int target1, int target2, CompMatr2 matrix);


/** @notyetdoced
 * 
 * Applies a multiply-controlled two-qubit dense unitary @p matrix to qubits 
 * @p target1 and @p target2 (treated as increasing significance) of @p qureg.
 * 
 * > - See applyCompMatr2() for information about @p target1, @p target2 and @p matrix.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
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
void applyMultiControlledCompMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, CompMatr2 matrix);


/** @notyetdoced
 * 
 * Applies an arbitrarily-controlled two-qubit dense unitary @p matrix to qubits 
 * @p target1 and @p target2 of @p qureg,
 * conditioned upon the @p controls being in the corresponding @p states.
 * 
 * > - See applyCompMatr2() for information about @p target1, @p target2 and @p matrix.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
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
 * @author Tyson Jones
 */
void applyMultiStateControlledCompMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, CompMatr2 matrix);


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
void applyMultiControlledCompMatr2(Qureg qureg, std::vector<int> controls, int target1, int target2, CompMatr2 matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledCompMatr2()
void applyMultiStateControlledCompMatr2(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target1, int target2, CompMatr2 matrix);


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
 * Applies an arbitrarily-sized dense unitary @p matrix to the
 * @p targets (treated as increasing significance) of @p qureg.
 * 
 * @formulae
 * 
 * Let @f$ M = @f$ @p matrix.
 * The qubits within @p targets are treated to be ordered least to most significant with respect
 * to @f$ M @f$. That is, if @f$ M @f$ was hypothetically separable single-qubit matrices
 * @f[
      M \equiv \dots \otimes C \otimes B \otimes A
 * @f]
 * then this function would effect
 * @f[
      \hat{M}_{\text{targets}} \equiv A_{\text{targets}[0]} \cdot B_{\text{targets}[1]} \cdot C_{\text{targets}[2]} \cdot \dots
 * @f]
 *
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order.
 *
 * @see
 * - applyCompMatr1()
 * - leftapplyCompMatr()
 * - rightapplyCompMatr()
 * @author Tyson Jones
 */
void applyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix);


/** @notyetdoced
 * 
 * Applies a singly-controlled arbitrarily-sized dense unitary @p matrix to the
 * @p targets (treated as increasing significance) of @p qureg.
 * 
 * > - See applyCompMatr() for information about @p targets and @p matrix.
 * > - See applyControlledCompMatr1() for information about the @p control qubit.
 * 
 * @author Tyson Jones
 */
void applyControlledCompMatr(Qureg qureg, int control, int* targets, int numTargets, CompMatr matrix);


/** @notyetdoced
 * 
 * Applies a multiply-controlled arbitrarily-sized dense unitary @p matrix upon the
 * @p targets of @p qureg.
 * 
 * > - See applyCompMatr() for information about @p targets and @p matrix.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledCompMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, CompMatr matrix);


/** @notyetdoced
 * 
 * Applies an arbitrarily-controlled arbitrarily-sized dense unitary @p matrix upon the
 * @p targets of @p qureg,
 * conditioned upon the @p controls being in the corresponding @p states.
 * 
 * > - See applyCompMatr() for information about @p targets and @p matrix.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledCompMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, CompMatr matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyCompMatr()
void applyCompMatr(Qureg qureg, std::vector<int> targets, CompMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyControlledCompMatr()
void applyControlledCompMatr(Qureg qureg, int control, std::vector<int> targets, CompMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledCompMatr()
void applyMultiControlledCompMatr(Qureg qureg, std::vector<int> controls, std::vector<int> targets, CompMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledCompMatr()
void applyMultiStateControlledCompMatr(Qureg qureg, std::vector<int> controls, std::vector<int> states, std::vector<int> targets, CompMatr matrix);


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


/** @notyetdoced
 * 
 * Applies a one-qubit diagonal unitary @p matrix to the @p target qubit of @p qureg.
 * 
 * @see 
 * - applyCompMatr1()
 * - leftapplyCompMatr2()
 * - rightapplyCompMatr2()
 */
void applyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix);


/** @notyetdoced
 * 
 * Applies a singly-controlled one-qubit diagonal unitary @p matrix to the
 * @p target qubit of @p qureg.
 * 
 * > - See applyDiagMatr1() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about the @p control qubit.
 * 
 * @author Tyson Jones
 */
void applyControlledDiagMatr1(Qureg qureg, int control, int target, DiagMatr1 matrix);


/** @notyetdoced
 *
 * Applies a multiply-controlled one-qubit diagonal unitary @p matrix upon the
 * @p target qubit of @p qureg.
 * 
 * > - See applyDiagMatr1() for information about @p target and @p matrix.
 * > - See applyMultiControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledDiagMatr1(Qureg qureg, int* controls, int numControls, int target, DiagMatr1 matrix);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled one-qubit diagonal unitary @p matrix upon the
 * @p target qubit of @p qureg,
 * conditioned upon the @p controls being in the corresponding @p states.
 * 
 * > - See applyDiagMatr1() for information about @p target and @p matrix.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledDiagMatr1(Qureg qureg, int* controls, int* states, int numControls, int target, DiagMatr1 matrix);


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
void applyMultiControlledDiagMatr1(Qureg qureg, std::vector<int> controls, int target, DiagMatr1 matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledDiagMatr1()
void applyMultiStateControlledDiagMatr1(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, DiagMatr1 matrix);


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


/** @notyetdoced
 * 
 * Applies a two-qubit diagonal unitary @p matrix to qubits
 * @p target1 and @p target2 of @p qureg.
 *
 * @author Tyson Jones
 */
void applyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix);


/** @notyetdoced
 * 
 * Applies a singly-controlled two-qubit diagonal unitary @p matrix to qubits
 * @p target1 and @p target2 of @p qureg.
 * 
 * > - See applyDiagMatr2() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about the @p control qubit.
 * 
 * @author Tyson Jones
 */
void applyControlledDiagMatr2(Qureg qureg, int control, int target1, int target2, DiagMatr2 matrix);


/** @notyetdoced
 *
 * Applies a multiply-controlled two-qubit diagonal unitary @p matrix upon
 * qubits @p target1 and @p target2 of @p qureg.
 * 
 * > - See applyDiagMatr2() for information about @p target1, @p target2 and @p matrix.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledDiagMatr2(Qureg qureg, int* controls, int numControls, int target1, int target2, DiagMatr2 matrix);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled two-qubit diagonal unitary @p matrix upon
 * qubits @p target1 and @p target2 of @p qureg,
 * conditioned upon the @p controls being in the corresponding @p states.
 * 
 * > - See applyDiagMatr2() for information about @p target1, @p target2 and @p matrix.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledDiagMatr2(Qureg qureg, int* controls, int* states, int numControls, int target1, int target2, DiagMatr2 matrix);


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
void applyMultiControlledDiagMatr2(Qureg qureg, std::vector<int> controls, int target1, int target2, DiagMatr2 matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledDiagMatr2()
void applyMultiStateControlledDiagMatr2(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target1, int target2, DiagMatr2 matrix);


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


/** @notyetdoced
 *
 * Applies an arbitrarily-sized diagonal unitary @p matrix upon the @p targets of @p qureg.
 * 
 * > [!TIP]
 * > To efficiently apply a diagonal matrix upon _all_ targets of @p qureg, 
 * > use applyFullStateDiagMatr().
 * 
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order.
 * 
 * @see
 * - applyDiagMatrPower()
 * - applyControlledDiagMatr()
 * - applyFullStateDiagMatr()
 * @author Tyson Jones
 */
void applyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/** @notyetdoced
 * 
 * Applies a singly-controlled arbitrarily-sized diagonal unitary @p matrix
 * upon the @p targets of @p qureg.
 * 
 * > - See applyDiagMatr() for information about the base operation, @p targets and @p matrix.
 * > - See applyControlledCompMatr1() for information about the @p control qubit.
 * 
 * @author Tyson Jones
 */
void applyControlledDiagMatr(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix);


/** @notyetdoced
 *
 * Applies a multiply-controlled arbitrarily-sized diagonal unitary @p matrix upon
 * the @p targets of @p qureg.
 * 
 * > - See applyDiagMatr() for information about @p targets and @p matrix.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledDiagMatr(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled arbitrarily-sized diagonal unitary @p matrix upon
 * the @p targets of @p qureg,
 * conditioned upon the @p controls being in the corresponding @p states.
 * 
 * > - See applyDiagMatr() for information about @p targets and @p matrix.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledDiagMatr(Qureg qureg, int* controls, int* states, int numControls, int* targets, int numTargets, DiagMatr matrix);


/** @notyetdoced
 * 
 * Applies an arbitrarily-sized diagonal unitary @p matrix, raised to the power @p exponent,
 * upon the @p targets of @p qureg.
 * 
 * > [!TIP]
 * > To efficiently apply an exponentiated diagonal matrix upon _all_ targets of @p qureg, 
 * > use applyFullStateDiagMatrPower().
 * 
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order.
 *
 * @formulae
 * 
 * This function is equivalent to applyDiagMatr() except that @p matrix is raised to the given @p exponent.
 * 
 * @see
 * - applyControlledDiagMatrPower()
 * - applyFullStateDiagMatrPower()
 */
void applyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/** @notyetdoced
 * 
 * Applies a singly-controlled arbitrarily-sized diagonal unitary @p matrix, 
 * raised to the power @p exponent, upon the @p targets of @p qureg.
 * 
 * > - See applyDiagMatr() for information about @p targets and @p matrix.
 * > - See applyDiagMatrPower() for information about @p exponent.
 * > - See applyControlledCompMatr1() for information about the @p control qubit.
 * 
 * @author Tyson Jones
 */
void applyControlledDiagMatrPower(Qureg qureg, int control, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/** @notyetdoced
 *
 * Applies a multiply-controlled arbitrarily-sized diagonal unitary @p matrix,
 * raised to the power @p exponent, 
 * upon the @p targets of @p qureg.
 * 
 * > - See applyDiagMatr() for information about @p targets and @p matrix.
 * > - See applyDiagMatrPower() for information about @p exponent.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledDiagMatrPower(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled arbitrarily-sized diagonal unitary @p matrix,
 * raised to the power @p exponent, 
 * upon the @p targets of @p qureg,
 * conditioned upon the @p controls being in the corresponding @p states.
 * 
 * > - See applyDiagMatr() for information about @p targets and @p matrix.
 * > - See applyDiagMatrPower() for information about @p exponent.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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


/** @notyetdoced
 * 
 * Applies a singly-controlled S gate on the @p target qubit of @p qureg.
 * 
 * > - See applyS() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledS(Qureg qureg, int control, int target);


/** @notyetdoced
 *
 * Applies a multiply-controlled S gate on the @p target qubit of @p qureg.
 * 
 * > - See applyS() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledS(Qureg qureg, int* controls, int numControls, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled S gate on the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyS() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledS(Qureg qureg, int* controls, int* states, int numControls, int target);


/// @notyetdoced
void applyT(Qureg qureg, int target);


/** @notyetdoced
 * 
 * Applies a singly-controlled T gate on the @p target qubit of @p qureg.
 * 
 * > - See applyT() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledT(Qureg qureg, int control, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled T gate on the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyT() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledT(Qureg qureg, int* controls, int numControls, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled T gate on the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyT() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledT(Qureg qureg, int* controls, int* states, int numControls, int target);


/// @notyetdoced
void applyHadamard(Qureg qureg, int target);


/** @notyetdoced
 * 
 * Applies a singly-controlled Hadamard gate on the @p target qubit of @p qureg.
 * 
 * > - See applyHadamard() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledHadamard(Qureg qureg, int control, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled Hadamard gate on the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyHadamard() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledHadamard(Qureg qureg, int* controls, int numControls, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled Hadamard gate on the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyHadamard() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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


/** @notyetdoced
 * 
 * Applies a singly-controlled SWAP gate upon @p qubit1 and @p qubit2 of @p qureg.
 * 
 * > - See applySwap() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledSwap(Qureg qureg, int control, int qubit1, int qubit2);


/** @notyetdoced
 *
 * Applies a multiply-controlled SWAP gate on @p qubit1 and @p qubit2 of @p qureg.
 * 
 * > - See applySwap() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled SWAP gate on @p qubit1 and @p qubit2 of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applySwap() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledSwap(Qureg qureg, int* controls, int* states, int numControls, int qubit1, int qubit2);


/// @notyetdoced
void applySqrtSwap(Qureg qureg, int qubit1, int qubit2);


/** @notyetdoced
 * 
 * Applies a singly-controlled square-root-of-SWAP gate upon @p qubit1 and @p qubit2 of @p qureg.
 * 
 * > - See applySqrtSwap() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledSqrtSwap(Qureg qureg, int control, int qubit1, int qubit2);


/** @notyetdoced
 *
 * Applies a multiply-controlled square-root-of-SWAP gate on @p qubit1 and @p qubit2 of @p qureg.
 * 
 * > - See applySqrtSwap() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledSqrtSwap(Qureg qureg, int* controls, int numControls, int qubit1, int qubit2);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled square-root-of-SWAP gate on @p qubit1 and @p qubit2 of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applySqrtSwap() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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
void applyMultiStateControlledSqrtSwap(Qureg qureg, std::vector<int> controls, std::vector<int> states, int qubit1, int qubit2);


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
void applyPauliX(Qureg qureg, int target);


/// @notyetdoced
void applyPauliY(Qureg qureg, int target);


/// @notyetdoced
void applyPauliZ(Qureg qureg, int target);


/** @notyetdoced
 * 
 * Applies a singly-controlled Pauli @c X operator (or NOT gate) upon the @p target qubit of @p qureg.
 * 
 * > - See applyPauliX() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledPauliX(Qureg qureg, int control, int target);


/** @notyetdoced
 * 
 * Applies a singly-controlled Pauli @c Y operator upon the @p target qubit of @p qureg.
 * 
 * > - See applyPauliY() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledPauliY(Qureg qureg, int control, int target);


/** @notyetdoced
 * 
 * Applies a singly-controlled Pauli @c Z operator upon the @p target qubit of @p qureg.
 * 
 * > - See applyPauliZ() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledPauliZ(Qureg qureg, int control, int target);


/** @notyetdoced
 *
 * Applies a multiply-controlled Pauli @c X operator (or NOT gate) upon the @p target qubit of @p qureg.
 * 
 * > - See applyPauliX() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledPauliX(Qureg qureg, int* controls, int numControls, int target);


/** @notyetdoced
 *
 * Applies a multiply-controlled Pauli @c Y operator (or NOT gate) upon the @p target qubit of @p qureg.
 * 
 * > - See applyPauliY() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledPauliY(Qureg qureg, int* controls, int numControls, int target);


/** @notyetdoced
 *
 * Applies a multiply-controlled Pauli @c Z operator (or NOT gate) upon the @p target qubit of @p qureg.
 * 
 * > - See applyPauliZ() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledPauliZ(Qureg qureg, int* controls, int numControls, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled Pauli @c X operator (or NOT gate) upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyPauliX() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledPauliX(Qureg qureg, int* controls, int* states, int numControls, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled Pauli @c Y operator upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyPauliY() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledPauliY(Qureg qureg, int* controls, int* states, int numControls, int target);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled Pauli @c Z operator upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyPauliZ() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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
void applyPauliStr(Qureg qureg, PauliStr str);


/** @notyetdoced
 *
 * Applies a singly-controlled tensor product of Pauli operators @p str upon @p qureg.
 * 
 * > - See applyPauliStr() for information about @p str.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledPauliStr(Qureg qureg, int control, PauliStr str);


/** @notyetdoced
 *
 * Applies a multiply-controlled tensor product of Pauli operators @p str upon @p qureg.
 * 
 * > - See applyPauliStr() for information about @p str.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledPauliStr(Qureg qureg, int* controls, int numControls, PauliStr str);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled tensor product of Pauli operators @p str upon @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyPauliStr() for information about @p str.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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
 * upon the @p target qubit, where @f$ \hat{\sigma}_x @f$ is the Pauli @c X matrix.
 *
 * @equivalences
 * 
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
 * upon the @p target qubit, where @f$ \hat{\sigma}_y @f$ is the Pauli @c Y matrix.
 *
 * @equivalences
 * 
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
 * upon the @p target qubit, where @f$ \hat{\sigma}_z @f$ is the Pauli @c Z matrix.
 *
 * @equivalences
 * 
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


/** @notyetdoced
 *
 * Applies a singly-controlled one-qubit rotation of @p angle around the @c X axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateX() for information about the base operation, and @p angle.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledRotateX(Qureg qureg, int control, int target, qreal angle);


/** @notyetdoced
 *
 * Applies a singly-controlled one-qubit rotation of @p angle around the @c Y axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateY() for information about the base operation, and @p angle.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledRotateY(Qureg qureg, int control, int target, qreal angle);


/** @notyetdoced
 *
 * Applies a singly-controlled one-qubit rotation of @p angle around the @c Z axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateZ() for information about the base operation, and @p angle.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledRotateZ(Qureg qureg, int control, int target, qreal angle);


/** @notyetdoced
 *
 * Applies a multiply-controlled one-qubit rotation of @p angle around the @c X axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateX() for information about the base operation, and @p angle.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledRotateX(Qureg qureg, int* controls, int numControls, int target, qreal angle);


/** @notyetdoced
 *
 * Applies a multiply-controlled one-qubit rotation of @p angle around the @c Y axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateY() for information about the base operation, and @p angle.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledRotateY(Qureg qureg, int* controls, int numControls, int target, qreal angle);


/** @notyetdoced
 *
 * Applies a multiply-controlled one-qubit rotation of @p angle around the @c Z axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateZ() for information about the base operation, and @p angle.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledRotateZ(Qureg qureg, int* controls, int numControls, int target, qreal angle);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled one-qubit rotation of @p angle around the @c X axis,
 * upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyRotateX() for information about the base operation, and @p angle.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledRotateX(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled one-qubit rotation of @p angle around the @c Y axis,
 * upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyRotateY() for information about the base operation, and @p angle.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledRotateY(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled one-qubit rotation of @p angle around the @c Z axis,
 * upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyRotateZ() for information about the base operation, and @p angle.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledRotateZ(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle);


/** @notyetdoced
 * 
 * Rotates the @p target qubit of @p qureg by @p angle around an arbitrary axis specified by
 * vector @p axisX, @p axisY, @p axisZ.
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
 * 
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


/** @notyetdoced
 *
 * Applies a singly-controlled one-qubit rotation of @p angle around an arbitrary axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateAroundAxis() for information about the base operation, @p angle, @p axisX, @p axisY and @p axisZ.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledRotateAroundAxis(Qureg qureg, int control, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/** @notyetdoced
 *
 * Applies a multiply-controlled one-qubit rotation of @p angle around an arbitrary axis,
 * upon the @p target qubit of @p qureg.
 * 
 * > - See applyRotateAroundAxis() for information about the base operation, @p angle, @p axisX, @p axisY and @p axisZ.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledRotateAroundAxis(Qureg qureg, int* controls, int numControls, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled one-qubit rotation of @p angle around an arbitrary axis,
 * upon the @p target qubit of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyRotateAroundAxis() for information about the base operation, @p angle, @p axisX, @p axisY and @p axisZ.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
void applyMultiStateControlledRotateAroundAxis(Qureg qureg, int* controls, int* states, int numControls, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


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
void applyMultiControlledRotateAroundAxis(Qureg qureg, std::vector<int> controls, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledRotateAroundAxis()
void applyMultiStateControlledRotateAroundAxis(Qureg qureg, std::vector<int> controls, std::vector<int> states, int target, qreal angle, qreal axisX, qreal axisY, qreal axisZ);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup op_pauligadget PauliStr gadgets
 * @brief Functions for applying many-qubit rotations around arbitrary PauliStr.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 * 
 * Applies a many-qubit rotation upon @p qureg, generated by tensor product of Pauli operators @p str.
 * 
 * @formulae
 * 
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
 * 
 * - Because @f$ R_{\hat{\sigma}}(\theta) @f$ satisfies
 *   @f[
        R_{\hat{\sigma}}(\theta) \equiv 
          \cos\left( \frac{\theta}{2} \right) \, \id 
          - \iu  \sin\left( \frac{\theta}{2} \right) \, \hat{\sigma},
 *   @f]
 *   this function is equivalent to (but much faster than) effecting @f$ \hat{\sigma} @f$
 *   upon a clone which is subsequently combined.
 *   ```
     // prepare |temp> = str |qureg>
     Qureg temp = createCloneQureg(qureg);
     applyPauliStr(temp, str);

     // set |qureg> = cos(theta/2) |qureg> - i sin(theta/2) str |qureg>
     qcomp coeffs[] = {cos(theta/2), -1i * sin(theta/2)};
     Qureg quregs[] = {qureg, temp};
     setQuregToWeightedSum(qureg, coeffs, quregs, 2);
 *   ```
 * - When @p str contains only @f$ \hat{Z} @f$ or @f$ \id @f$ Paulis, this function will
 *   automatically invoke applyPhaseGadget() which leverages an optimised implementation.
 * - When @p str contains only @f$ \id @f$ Paulis, this function merely effects a change
 *   of global phase upon statevectors of @f$ -\theta/2 @f$, leaving density matrices
 *   unchanged.
 *   ```
     qcomp factor = cexp(- theta / 2 * 1.i);
     setQuregToWeightedSum(qureg, &factor, &qureg, 1);
 *   ```
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
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
 *
 * @see
 *  - applyNonUnitaryPauliGadget()
 */
void applyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


/** @notyetdoced
 * 
 * This function generalises applyPauliGadget() to accept a complex angle.
 */
void applyNonUnitaryPauliGadget(Qureg qureg, PauliStr str, qcomp angle);


/** @notyetdoced
 *
 * Applies a singly-controlled many-qubit rotation upon @p qureg, 
 * generated by tensor product of Pauli operators @p str.
 * 
 * > - See applyPauliGadget() for information about the base operation, @p angle, and @p str.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledPauliGadget(Qureg qureg, int control, PauliStr str, qreal angle);


/** @notyetdoced
 *
 * Applies a multiply-controlled many-qubit rotation upon @p qureg, 
 * generated by tensor product of Pauli operators @p str.
 * 
 * > - See applyPauliGadget() for information about the base operation, @p angle, and @p str.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledPauliGadget(Qureg qureg, int* controls, int numControls, PauliStr str, qreal angle);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled many-qubit rotation upon @p qureg, generated by tensor product of Pauli operators @p str,
 * and conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyPauliGadget() for information about the base operation, @p angle, and @p str.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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
 * @brief Functions for applying many-qubit rotations around the Pauli @c Z axis, and phase flips and shifts.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 * 
 * Applies a many-qubit @c Z rotation upon @p qureg, generated by a tensor product of Pauli @c Z operators
 * upon @p targets.
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
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order,
 * > though the effect of this function is incidentally unaffected by the ordering of @p targets.
 *
 * @equivalences
 * 
 * - This function is equivalent to calling applyPauliGadget() with a PauliStr containing only @f$ \hat{Z} @f$ and @f$ \id @f$.
 *   This latter function will actually automatically invoke applyPhaseGadget() which has an optimised implementation.
 * - This function is equivalent to, albeit much faster than, preparing a DiagMatr with @f$ \pm 1 @f$ elements (depending upon
 *   the parity of the targeted set bits) and effecting it with applyDiagMatr().
 * - Passing @p angle=0 is equivalent to effecting the identity, leaving the state unchanged.
 */
void applyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


/** @notyetdoced
 *
 * Applies a singly-controlled many-qubit @c Z rotation upon @p qureg, 
 * generated by a tensor product of Pauli @c Z operators upon @p targets.
 * 
 * > - See applyPhaseGadget() for information about the base operation, @p angle, and @p targets.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledPhaseGadget(Qureg qureg, int control, int* targets, int numTargets, qreal angle);


/** @notyetdoced
 *
 * Applies a multiply-controlled many-qubit @c Z rotation upon @p qureg, 
 * generated by a tensor product of Pauli @c Z operators upon @p targets.
 * 
 * > - See applyPhaseGadget() for information about the base operation, @p angle, and @p targets.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledPhaseGadget(Qureg qureg, int* controls, int numControls, int* targets, int numTargets, qreal angle);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled many-qubit @c Z rotation upon @p qureg, generated by a tensor product of Pauli @c Z operators
 * upon @p targets, and conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyPhaseGadget() for information about the base operation, @p angle, and @p targets.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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
 * 
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
 * 
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
 * 
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
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order,
 * > though the effect of this function is incidentally unaffected by the ordering of @p targets.
 * 
 * @equivalences
 * 
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
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order,
 * > though the effect of this function is incidentally unaffected by the ordering of @p targets.
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
 * 
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
 * @defgroup op_nots Many-not gates
 * @brief Functions for effecting many-qubit NOT gates
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 *
 * Applies a many-qubit NOT gate (or tensor product of Pauli @c X operators) upon the @p targets of @p qureg.
 *
 * > [!TIP]
 * > This function is sometimes more efficient when @p targets are specified in increasing order,
 * > though the effect of this function is incidentally unaffected by the ordering of @p targets.
 * 
 * @author Tyson Jones
 */
void applyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


/** @notyetdoced
 *
 * Applies a singly-controlled many-qubit NOT gate (or tensor product of Pauli @c X operators) 
 * upon the @p targets of @p qureg.
 * 
 * > - See applyMultiQubitNot() for information about the base operation.
 * > - See applyControlledCompMatr1() for information about @p control.
 * 
 * @author Tyson Jones
 */
void applyControlledMultiQubitNot(Qureg qureg, int control, int* targets, int numTargets);


/** @notyetdoced
 *
 * Applies a multiply-controlled many-qubit NOT gate (or tensor product of Pauli @c X operators)
 * upon the @p targets of @p qureg.
 * 
 * > - See applyMultiQubitNot() for information about the base operation.
 * > - See applyMultiControlledCompMatr1() for information about @p controls.
 * 
 * @author Tyson Jones
 */
void applyMultiControlledMultiQubitNot(Qureg qureg, int* controls, int numControls, int* targets, int numTargets);


/** @notyetdoced
 *
 * Applies an arbitrarily-controlled many-qubit NOT gate (or tensor product of Pauli @c X operators) upon the @p targets of @p qureg,
 * conditioned upon @p controls being in the corresponding @p states.
 * 
 * > - See applyMultiQubitNot() for information about the base operation.
 * > - See applyMultiStateControlledCompMatr1() for information about @p controls and @p states.
 * 
 * @author Tyson Jones
 */
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
int applyQubitMeasurement(Qureg qureg, int target);


/// @notyetdoced
int applyQubitMeasurementAndGetProb(Qureg qureg, int target, qreal* probability);


/// @notyetdoced
qreal applyForcedQubitMeasurement(Qureg qureg, int target, int outcome);


/// @notyetdoced
qindex applyMultiQubitMeasurement(Qureg qureg, int* qubits, int numQubits);


/// @notyetdoced
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, int* qubits, int numQubits, qreal* probability);


/// @notyetdoced
qreal applyForcedMultiQubitMeasurement(Qureg qureg, int* qubits, int* outcomes, int numQubits);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitMeasurement()
qindex applyMultiQubitMeasurement(Qureg qureg, std::vector<int> qubits);


/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiQubitMeasurementAndGetProb()
qindex applyMultiQubitMeasurementAndGetProb(Qureg qureg, std::vector<int> qubits, qreal* probability);


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
void applyQubitProjector(Qureg qureg, int target, int outcome);


/// @notyetdoced
void applyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


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


/** 
 * Applies the Quantum Fourier Transform upon the specified @p targets of @p qureg.
 * Alternatively, applies the Inverse Quantum Fourier Transform according to @p inverse.
 * 
 * @formulae
 * 
 * Letting @f$ N @f$ = @p numTargets, the @f$ N @f$ qubit Quantum Fourier Transform maps each
 * computational basis state of the targeted qubits, @f$ \ket{j} @f$, according to
 * @f[ 
        \ket{j} \rightarrow \frac{1}{\sqrt{2^N}} \sum_{k=0}^{2^N-1} e^{2 \pi i j k / 2^N} \ket{k}.
 * @f]
 * Similarly the Inverse Quantum Fourier Transform maps each basis state like
 * @f[ 
        \ket{j} \rightarrow \frac{1}{\sqrt{2^N}} \sum_{k=0}^{2^N-1} e^{-2 \pi i j k / 2^N} \ket{k}.
 * @f]
 *
 * @param[in,out] qureg      the state to modify.
 * @param[in]     targets    the indices of the target qubits.
 * @param[in]     numTargets the length of list @p targets
 * @param[in]     inverse    whether to apply the inverse QFT or forward QFT
 * @throws @validationerror
 * - if @p qureg is uninitialised.
*  - if @p targets are invalid qubit indices.
*  - if @p targets are not unique.
 * - if @p numTargets < 1.
 * @see
 * - applyFullQuantumFourierTransform()
 * @author Vasco Ferreira
 */
void applyQuantumFourierTransform(Qureg qureg, int* targets, int numTargets, bool inverse);


/** 
 * Applies the Quantum Fourier Transform upon all qubits in @p qureg. Alternatively,
 * applies the Inverse Quantum Fourier Transform according to @p inverse.
 * 
 * @formulae
 * 
 * The Quantum Fourier Transform maps each computational basis state @f$ \ket{j} @f$
 * in an @f$ N @f$ qubit @p qureg according to
 * @f[ 
        \ket{j} \rightarrow \frac{1}{\sqrt{2^N}} \sum_{k=0}^{2^N-1} e^{2 \pi i j k / 2^N} \ket{k}.
 * @f]
 * Similarly the Inverse Quantum Fourier Transform maps each basis state like
 * @f[ 
        \ket{j} \rightarrow \frac{1}{\sqrt{2^N}} \sum_{k=0}^{2^N-1} e^{-2 \pi i j k / 2^N} \ket{k}.
 * @f]
 *
 * @equivalences
 *
 * - This function wraps applyQuantumFourierTransform(), passing all qubits in the @p qureg as targets.
 *
 * @param[in,out] qureg      the state to modify.
 * @param[in]     inverse    whether to apply the inverse QFT or forward QFT
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @see
 * - applyQuantumFourierTransform()
 * @author Vasco Ferreira
 */
void applyFullQuantumFourierTransform(Qureg qureg, bool inverse);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyetdoced
/// @cppvectoroverload
/// @see applyQuantumFourierTransform()
void applyQuantumFourierTransform(Qureg qureg, std::vector<int> targets, bool inverse);


#endif // __cplusplus

/** @} */



#endif // OPERATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
