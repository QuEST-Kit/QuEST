/** @file
 * API signatures for directly pre- and post-multiplying 
 * operators upon density matrices, likely constituting 
 * non-physical operations which break state normalisation.
 * 
 * @author Tyson Jones
 * 
 * @defgroup multiplication Multiplication
 * @ingroup api
 * @brief Functions for directly pre- or post-multiplying operators 
 *        upon density matrices.
 * @{
 */

#ifndef MULTIPLICATION_H
#define MULTIPLICATION_H

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
 * @defgroup mult_compmatr1 CompMatr1
 * @brief Functions for pre- or post-multiplying general one-qubit dense matrices
 *        (as CompMatr1) upon density matrices.
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

    leftapplyCompMatr1(qureg, 2, matrix); 
 * ```
 *
 * @param[in,out] qureg  the state to modify.
 * @param[in]     target the index of the target qubit.
 * @param[in]     matrix the Z-basis matrix to multiply upon the left.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p target is an invalid qubit index.
 * @see
 * - getCompMatr1()
 * - getInlineCompMatr1()
 * - applyCompMatr1()
 * - rightapplyCompMatr1()
 * - applyQubitProjector()
 * - leftapplyCompMatr()
 * @author Tyson Jones
 */
void leftapplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix);


/** Multiplies a general one-qubit dense @p matrix upon the specified @p target 
 * qubit of the density matrix @p qureg, from the right-hand side.
 *  
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ \hat{M} = @f$ @p matrix and @f$ t = @f$ @p target, 
 * and notate @f$\hat{M}_t@f$ as per applyCompMatr1(). Unlike applyCompMatr1() however,
 * this function only ever right-multiplies @p matrix upon @p qureg.
 * 
 * Explicitly
 *   @f[ 
        \dmrho \rightarrow \dmrho \, \hat{M}_t
 *   @f]
 * where @f$ \hat{M} @f$ is not conjugated nor transposed, and there are no additional 
 * constraints like unitarity.
 * 
 * In general, this function will break the normalisation of @p qureg and result in a
 * non-physical state, and is useful for preparing sub-expressions of formulae like
 * the Linbladian.
 *
 * @myexample
 * ```
    Qureg qureg = createDensityQureg(5);

    CompMatr1 matrix = getInlineCompMatr1({
        {0.1, 0.2},
        {0.3i, 0.4i}
    });

    rightapplyCompMatr1(qureg, 2, matrix); 
 * ```
 *
 * @param[in,out] qureg  the state to modify.
 * @param[in]     target the index of the target qubit.
 * @param[in]     matrix the Z-basis matrix to post-multiply.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p qureg is not a density matrix.
 * - if @p target is an invalid qubit index.
 * @see
 * - getCompMatr1()
 * - getInlineCompMatr1()
 * - applyCompMatr1()
 * - leftapplyCompMatr1()
 * - leftapplyCompMatr()
 * @author Tyson Jones
 */
void rightapplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_compmatr2 CompMatr2
 * @brief Functions for pre- or post-multiplying general two-qubit dense matrices
 *        (as CompMatr2) upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see
/// - applyCompMatr2()
/// - leftapplyCompMatr1()
void leftapplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);


/// @notyetdoced
/// @see
/// - rightapplyCompMatr1()
void rightapplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */




/** 
 * @defgroup mult_compmatr CompMatr
 * @brief Functions for pre- or post-multiplying general many-target dense matrices
 *        (as CompMatr) upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 * 
 * @see
 * - applyCompMatr()
 * - leftapplyCompMatr1()
 */
void leftapplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix);


/// @notyetdoced
/// @see
/// - rightapplyCompMatr1()
void rightapplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see leftapplyCompMatr()
void leftapplyCompMatr(Qureg qureg, std::vector<int> targets, CompMatr matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see rightapplyCompMatr()
void rightapplyCompMatr(Qureg qureg, std::vector<int> targets, CompMatr matr);


#endif 

/** @} */



/** 
 * @defgroup mult_diagmatr1 DiagMatr1
 * @brief Functions for pre- or post-multiplying general single-qubit diagonal 
 *        matrices (as DiagMatr1) upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see leftapplyCompMatr1()
void leftapplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);


/// @notyetdoced
/// @see rightapplyCompMatr1()
void rightapplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_diagmatr2 DiagMatr2
 * @brief Functions for pre- or post-multiplying general two-qubit diagonal 
 *        matrices (as DiagMatr2) upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see leftapplyCompMatr1()
void leftapplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);


/// @notyetdoced
/// @see rightapplyCompMatr1()
void rightapplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_diagmatr DiagMatr
 * @brief Functions for pre- or post-multiplying general any-target diagonal 
 *        matrices (as DiagMatr), or powers thereof, upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see leftapplyCompMatr1()
void leftapplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see rightapplyCompMatr1()
void rightapplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see
/// - leftapplyCompMatr1()
/// - applyDiagMatrPower()
void leftapplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @see 
/// - rightapplyCompMatr1()
/// - applyDiagMatrPower()
void rightapplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see leftapplyDiagMatr()
void leftapplyDiagMatr(Qureg qureg, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see rightapplyDiagMatr()
void rightapplyDiagMatr(Qureg qureg, std::vector<int> targets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see leftapplyDiagMatrPower()
void leftapplyDiagMatrPower(Qureg qureg, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see rightapplyDiagMatrPower()
void rightapplyDiagMatrPower(Qureg qureg, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


#endif 

/** @} */



/** 
 * @defgroup mult_fullstatediagmatr FullStateDiagMatr
 * @brief Functions for pre- or post-multiplying general full-state diagonal 
 *        matrices (FullStateDiagMatr), or powers thereof, upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - leftapplyCompMatr1()
void leftapplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - rightapplyCompMatr1()
/// - applyFullStateDiagMatr()
void rightapplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - leftapplyCompMatr1()
/// - applyFullStateDiagMatr()
void leftapplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - rightapplyCompMatr1()
/// - applyFullStateDiagMatr()
void rightapplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup multi_swap Swap
 * @brief Functions for pre- or post-multiplying the two-qubit SWAP
 *        gate upon density matrices
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applySwap()
void leftapplySwap(Qureg qureg, int qubit1, int qubit2);


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applySwap()
void rightapplySwap(Qureg qureg, int qubit1, int qubit2);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_pauli Pauli
 * @brief Functions for pre- or post-multiplying the individual one-qubit 
 *        Pauli operators upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyPauliX()
void leftapplyPauliX(Qureg qureg, int target);


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyPauliY()
void leftapplyPauliY(Qureg qureg, int target);


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyPauliZ()
void leftapplyPauliZ(Qureg qureg, int target);


/// @notyetdoced
/// @see 
/// - rightapplyCompMatr1()
/// - applyPauliX()
void rightapplyPauliX(Qureg qureg, int target);


/// @notyetdoced
/// @see 
/// - rightapplyCompMatr1()
/// - applyPauliY()
void rightapplyPauliY(Qureg qureg, int target);


/// @notyetdoced
/// @see 
/// - rightapplyCompMatr1()
/// - applyPauliZ()
void rightapplyPauliZ(Qureg qureg, int target);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_paulistr PauliStr
 * @brief Functions for pre- or post-multiplying a tensor product of 
 *       Pauli operators (as a PauliStr) upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyPauliStr()
void leftapplyPauliStr(Qureg qureg, PauliStr str);


/// @notyetdoced
/// @see 
/// - rightapplyCompMatr1()
/// - applyPauliStr()
void rightapplyPauliStr(Qureg qureg, PauliStr str);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_pauligadget Pauli gadgets
 * @brief Functions for pre- or post-multiplying many-qubit rotations around 
 *        arbitrary PauliStr upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyPauliGadget()
void leftapplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


/// @notyetdoced
/// @see 
/// - rightapplyCompMatr1()
/// - applyPauliGadget()
void rightapplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



/** 
 * @defgroup mult_phasegadget Phase gates
 * @brief Functions for pre- or post-multiplying many-qubit rotations around 
 *        the Pauli Z axis upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyPhaseGadget()
void leftapplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


/// @notyetdoced
/// @see
/// - rightapplyCompMatr1()
/// - applyPhaseGadget()
void rightapplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see leftapplyPhaseGadget()
void leftapplyPhaseGadget(Qureg qureg, std::vector<int> targets, qreal angle);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see rightapplyPhaseGadget()
void rightapplyPhaseGadget(Qureg qureg, std::vector<int> targets, qreal angle);


#endif

/** @} */



/** 
 * @defgroup mult_nots Many-not gates
 * @brief Functions for pre- or post-multiplying many-qubit NOT gates 
 *        upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @see 
/// - leftapplyCompMatr1()
/// - applyMultiQubitNot()
void leftapplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - rightapplyCompMatr1()
/// - applyMultiQubitNot()
void rightapplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see leftapplyMultiQubitNot()
void leftapplyMultiQubitNot(Qureg qureg, std::vector<int> targets);


/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see rightapplyMultiQubitNot()
void rightapplyMultiQubitNot(Qureg qureg, std::vector<int> targets);


#endif

/** @} */



/** 
 * @defgroup mult_projectors Projectors
 * @brief Functions for pre- or post-multiplying projectors upon density matrices.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - leftapplyCompMatr1()
/// - applyQubitProjector()
void leftapplyQubitProjector(Qureg qureg, int qubit, int outcome);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - leftapplyCompMatr1()
/// - applyMultiQubitProjector()
void leftapplyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - rightapplyCompMatr1()
/// - applyQubitProjector()
void rightapplyQubitProjector(Qureg qureg, int qubit, int outcome);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - rightapplyCompMatr1()
/// - applyMultiQubitProjector()
void rightapplyMultiQubitProjector(Qureg qureg, int* qubits, int* outcomes, int numQubits);


// end de-mangler
#ifdef __cplusplus
}
#endif



/** 
 * @defgroup mult_paulistrsum PauliStrSum
 * @brief Functions for pre- or post-multiplying weighted sums of Pauli 
 *        tensors upon a density matrix.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/// @notyetdoced
/// @notyetvalidated
/// @see leftapplyCompMatr1()
void leftapplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);


/// @notyetdoced
/// @notyetvalidated
/// @see leftapplyCompMatr1()
void rightapplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



#endif // MULTIPLICATION_H

/** @} */ // (end file-wide doxygen defgroup)
