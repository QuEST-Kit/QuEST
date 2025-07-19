/** @file
 * API signatures for directly pre- and post-multiplying 
 * operators upon density matrices, likely constituting 
 * non-physical operations which break state normalisation.
 * 
 * @author Tyson Jones
 * 
 * @defgroup multiplication Multiplication
 * @ingroup api
 * @brief Functions for directly multiplying operators upon 
 *        density matrices.
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
 * - postMultiplyCompMatr1()
 * - applyQubitProjector()
 * - multiplyCompMatr()
 * @author Tyson Jones
 */
void multiplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix);


/** @notyettested
 * 
 * Multiplies a general one-qubit dense @p matrix upon the specified @p target 
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

    postMultiplyCompMatr1(qureg, 2, matrix); 
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
 * - multiplyCompMatr1()
 * - multiplyCompMatr()
 * @author Tyson Jones
 */
void postMultiplyCompMatr1(Qureg qureg, int target, CompMatr1 matrix);


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
/// - multiplyCompMatr1()
void multiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matr);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
/// @see
/// - postMultiplyCompMatr1
void postMultiplyCompMatr2(Qureg qureg, int target1, int target2, CompMatr2 matrix);


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
 * - multiplyCompMatr1()
 */
void multiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
/// @see
/// - postMultiplyCompMatr1
void postMultiplyCompMatr(Qureg qureg, int* targets, int numTargets, CompMatr matrix);


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
/// @see postMultiplyCompMatr()
void postMultiplyCompMatr(Qureg qureg, std::vector<int> targets, CompMatr matr);


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
/// @see multiplyCompMatr1()
void multiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
void postMultiplyDiagMatr1(Qureg qureg, int target, DiagMatr1 matrix);


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
/// @see multiplyCompMatr1()
void multiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matr);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
void postMultiplyDiagMatr2(Qureg qureg, int target1, int target2, DiagMatr2 matrix);


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
/// @see multiplyCompMatr1()
void multiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
void postMultiplyDiagMatr(Qureg qureg, int* targets, int numTargets, DiagMatr matrix);


/// @notyetdoced
/// @see
/// - multiplyCompMatr1()
/// - applyDiagMatrPower()
void multiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
void postMultiplyDiagMatrPower(Qureg qureg, int* targets, int numTargets, DiagMatr matrix, qcomp exponent);


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
/// @see postMultiplyDiagMatr()
void postMultiplyDiagMatr(Qureg qureg, std::vector<int> targets, DiagMatr matrix);


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
/// @see postMultiplyDiagMatrPower()
void postMultiplyDiagMatrPower(Qureg qureg, std::vector<int> targets, DiagMatr matrix, qcomp exponent);


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
/// - multiplyCompMatr1
void multiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix);


/// @notyetdoced
/// @notyetvalidated
/// @see
/// - multiplyCompMatr1
/// - applyDiagMatrPower
void multiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


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
/// @see multiplyCompMatr1()
void multiplySwap(Qureg qureg, int qubit1, int qubit2);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplySwap(Qureg qureg, int qubit1, int qubit2);


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
/// @notyettested
/// @see multiplyCompMatr1()
void multiplyPauliX(Qureg qureg, int target);


/// @notyetdoced
/// @notyettested
/// @see multiplyCompMatr1()
void multiplyPauliY(Qureg qureg, int target);


/// @notyetdoced
/// @notyettested
/// @see multiplyCompMatr1()
void multiplyPauliZ(Qureg qureg, int target);


/// @notyetdoced
/// @notyettested
/// @see postMultiplyCompMatr1()
void postMultiplyPauliX(Qureg qureg, int target);


/// @notyetdoced
/// @notyettested
/// @see postMultiplyCompMatr1()
void postMultiplyPauliY(Qureg qureg, int target);


/// @notyetdoced
/// @notyettested
/// @see postMultiplyCompMatr1()
void postMultiplyPauliZ(Qureg qureg, int target);


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
/// @see multiplyCompMatr1()
void multiplyPauliStr(Qureg qureg, PauliStr str);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyPauliStr(Qureg qureg, PauliStr str);


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
/// - multiplyCompMatr1()
/// - applyPauliGadget()
void multiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyPauliGadget(Qureg qureg, PauliStr str, qreal angle);


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
/// - multiplyCompMatr1()
/// - applyPhaseGadget
void multiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyPhaseGadget(Qureg qureg, int* targets, int numTargets, qreal angle);


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
/// @see postMultiplyPhaseGadget()
void postMultiplyPhaseGadget(Qureg qureg, std::vector<int> targets, qreal angle);


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
/// @see multiplyCompMatr1()
void multiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyMultiQubitNot(Qureg qureg, int* targets, int numTargets);


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
/// @see postMultiplyMultiQubitNot()
void postMultiplyMultiQubitNot(Qureg qureg, std::vector<int> targets);


#endif

/** @} */



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
/// @see multiplyCompMatr1()
void multiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);


/// @notyetdoced
/// @notyettested
/// @notyetvalidated
void postMultiplyPauliStrSum(Qureg qureg, PauliStrSum sum, Qureg workspace);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



#endif // MULTIPLICATION_H

/** @} */ // (end file-wide doxygen defgroup)
