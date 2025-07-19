/** @file
 * API signatures for effecting Trotterised operators which
 * approximate the action of exponentials of PauliStrSum
 * 
 * @author Tyson Jones
 * 
 * @defgroup trotterisation Trotterisation
 * @ingroup api
 * @brief Functions for Trottersing operations upon Quregs.
 * @{
 */

#ifndef TROTTERISATION_H
#define TROTTERISATION_H

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#ifdef __cplusplus
    #include <vector>
#endif



/** 
 * @defgroup trotter_paulistrsum PauliStrSum gadgets
 * @brief Functions for using Trotterisation to approximate the action of 
 *        exponentials of weighted sums of Pauli tensors upon Quregs.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyettested
 * 
 * Effects (an approximation to) the exponential of @p sum, weighted by @p angle, upon @p qureg,
 * via the symmetrized Trotter-Suzuki decomposition (<a href="https://arxiv.org/abs/math-ph/0506007">arXiv</a>).
 * Increasing @p reps (the number of Trotter repetitions) or @p order (an even, positive integer or one) 
 * improves the accuracy of the approximation (reducing the "Trotter error" due to non-commuting 
 * terms of @p sum), though increases the runtime linearly and exponentially respectively.
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p sum and @f$ \theta = @f$ @p angle. This function approximates the action of
 * @f[
      \exp \left(\iu \, \theta \, \hat{H} \right)
 * @f]
 * via a Trotter-Suzuki decomposition of the specified @p order and number of repetitions (@p reps).
 * Simulation is exact, regardless of @p order or @p reps, only when all terms in @p sum commute.
 * 
 * @important
 *   Note that @f$ \theta @f$ lacks the @f$ -\frac{1}{2} @f$ prefactor present in other functions like
 *   applyPauliGadget().
 * 
 * To be precise, let @f$ r = @f$ @p reps and assume @p sum is composed of
 * @f$ T @f$-many terms of the form
 * @f[
      \hat{H} = \sum\limits_j^T c_j \, \hat{\sigma}_j
 * @f]
 * where @f$ c_j @f$ is the coefficient of the @f$ j @f$-th PauliStr @f$ \hat{\sigma}_j @f$.
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
 * 
 * @equivalences
 * 
 * - Time evolution of duration @f$ t @f$ under a time-independent Hamiltonian @p sum = @f$ \hat{H} @f$, as
 *   per the unitary time evolution operator
 *   @f[
        \hat{U}(t) = \exp(- \iu \, t  \,\hat{H} \, / \, \hbar) 
 *   @f]
 *   is approximated via @f$ \theta = - t / \hbar @f$.
 *   ```
     qreal time = 3.14;
     qreal angle = - time / hbar;
     applyTrotterizedPauliStrSumGadget(qureg, sum, angle, order, reps);
 *   ```
 * - This function is equivalent to applyNonUnitaryTrotterizedPauliStrSumGadget() when passing
 *   a @p qcomp instance with a zero imaginary component as the @p angle parameter. This latter 
 *   function is useful for generalising dynamical simulation to imaginary-time evolution.
 * 
 * @constraints
 * - Unitarity of the prescribed exponential(s) requires that @p sum is Hermitian, ergo containing
 *   only real coefficients. Validation will check that @p sum is approximately Hermitian, permitting
 *   coefficients with imaginary components smaller (in magnitude) than epsilon.
 *   @f[ 
        \max\limits_{i} \Big|c_i| \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setValidationEpsilon().
 *   Otherwise, use applyNonUnitaryTrotterizedPauliStrSumGadget() to permit non-Hermitian @p sum
 *   and ergo effect a non-unitary exponential(s). 
 * - The @p angle parameter is necessarily real despite the validation epsilon, but can be relaxed
 *   to an arbitrary complex scalar using applyNonUnitaryTrotterizedPauliStrSumGadget().
 * - This function only ever effects @f$ \exp \left(\iu \, \theta \, \hat{H} \right) @f$ exactly
 *   when all PauliStr in @p sum = @f$ \hat{H} @f$ commute. 
 * 
 * @param[in,out] qureg  the state to modify.
 * @param[in]     sum    a weighted sum of Pauli strings to approximately exponentiate.
 * @param[in]     angle  an effective prefactor of @p sum in the exponent.
 * @param[in]     order  the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...)
 * @param[in]     reps   the number of Trotter repetitions
 * 
 * @throws @validationerror
 * - if @p qureg or @p sum are uninitialised.
 * - if @p sum is not approximately Hermitian.
 * - if @p sum contains non-identities on qubits beyond the size of @p qureg.
 * - if @p order is not 1 nor a positive, @b even integer.
 * - if @p reps is not a positive integer.
 * 
 * @see
 *  - applyPauliGadget()
 *  - applyNonUnitaryTrotterizedPauliStrSumGadget()
 * 
 * @author Tyson Jones
 */
void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyetdoced
/// @notyettested
/// @see
///  - applyTrotterizedPauliStrSumGadget()
///  - applyControlledCompMatr1()
void applyControlledTrotterizedPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyetdoced
/// @notyettested
/// @see
///  - applyTrotterizedPauliStrSumGadget()
///  - applyMultiControlledCompMatr1()
void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyetdoced
/// @notyettested
/// @see
///  - applyTrotterizedPauliStrSumGadget()
///  - applyMultiStateControlledCompMatr1()
void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps);


/** @notyettested
 * 
 * A generalisation of applyTrotterizedPauliStrSumGadget() which accepts a complex angle and permits
 * @p sum to be non-Hermitian, thereby effecting a potentially non-unitary and non-CPTP operation.
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p sum and @f$ \theta = @f$ @p angle. This function approximates the action of
 * @f[
      \exp \left(\iu \, \theta \, \hat{H} \right)
 * @f]
 * via a Trotter-Suzuki decomposition of the specified @p order and number of repetitions (@p reps). 
 * 
 * See applyTrotterizedPauliStrSumGadget() for more information about the decomposition.
 *
 * @equivalences
 * 
 * - When @p angle is set to @f$ \theta = \iu \, \tau @f$ and @p sum = @f$ \hat{H} @f$ is Hermitian,
 *   this function (approximately) evolves @p qureg in imaginary-time. That is, letting 
 *   @f$ \hat{U}(t) = \exp(-\iu \, t \, \hat{H}) @f$ be the normalised unitary evolution operator, this 
 *   function effects the imaginary-time operator
     @f[
        \hat{V}(\tau) = \hat{U}(t=-\iu \tau) = \exp(- \tau \hat{H}).
 *   @f]
 *   This operation drives the system toward the (unnormalised) groundstate.
 *   Let @f$ \{ \ket{\phi_i} \} @f$ and @f$ \{ \ket{\lambda_i} \} @f$ be the eigenstates and respective
 *   eigenvalues of @f$ \hat{H} @f$, which are real due to Hermiticity.
 *   @f[
         \hat{H} = \sum \limits_i \lambda_i \ket{\phi_i}\bra{\phi_i},
         \;\;\;\;\; \lambda_i \in \mathbb{R}.
 *   @f]
 *   
 *   - When @p qureg is a statevector @f$ \svpsi @f$ and can ergo be expressed in the basis of 
 *     @f$ \{ \ket{\phi_i} \} @f$ as @f$ \svpsi = \sum_i \alpha_i \ket{\phi_i} @f$, 
 *     this function approximates
 *     @f[
          \svpsi \, \rightarrow  \, \hat{V}(\tau) \svpsi =
          \sum\limits_i \alpha_i \exp(- \tau \, \lambda_i) \ket{\phi_i}.
 *     @f]
 *   - When @p qureg is a density matrix and is ergo expressible as
 *     @f$ \dmrho = \sum\limits_{ij} \alpha_{ij} \ket{\phi_i}\bra{\phi_j} @f$, this function effects
 *     @f[
          \dmrho \, \rightarrow \, \hat{V}(\tau) \dmrho \hat{V}(\tau)^\dagger =
          \sum\limits_{ij} \alpha_{ij} \exp(-\tau (\lambda_i + \lambda_j)) \ket{\phi_i}\bra{\phi_j}.
 *     @f]
 *
 *   As @f$ \tau \rightarrow \infty @f$, the resulting unnormalised state approaches statevector
 *   @f$ \svpsi \rightarrow \alpha_0 \exp(-\tau \lambda_0) \ket{\phi_0} @f$ or density matrix
 *   @f$ \dmrho \rightarrow \alpha_{0,0} \exp(-2 \tau \lambda_0) \ket{\phi_0}\bra{\phi_0} @f$,
 *   where @f$ \lambda_0 @f$ is the minimum eigenvalue and @f$ \ket{\phi_0} @f$ is the groundstate.
 *   Assuming the initial overlap @f$ \alpha_0 @f$ is not zero (or exponentially tiny), 
 *   subsequent renormalisation via setQuregToRenormalized() produces the pure 
 *   ground-state @f$ \ket{\phi_0} @f$.
 *
 *   ```
     // pray for a non-zero initial overlap
     initRandomPureState(qureg); // works even for density matrices

     // minimize then renormalise
     qreal tau = 10; // impatient infinity
     int order = 4;
     int reps = 100;
     applyNonUnitaryTrotterizedPauliStrSumGadget(qureg, hamil, tau * 1i, order, reps);
     setQuregToRenormalized(qureg);

     // ground-state (phi_0)
     reportQureg(qureg);

     // lowest lying eigenvalue (lambda_0)
     qreal expec = calcExpecPauliStrSum(qureg, hamil);
     reportScalar("expec", expec);
 *   ```
 *
 *   Note degenerate eigenvalues will yield a pure superposition of the corresponding eigenstates, with 
 *   coefficients informed by the initial, relative populations.
 * 
 * - When @p angle is real and @p sum is Hermitian (has approximately real coefficients), this
 *   function is equivalent to applyTrotterizedPauliStrSumGadget()
 * 
 * @constraints
 * - This function only ever effects @f$ \exp \left(\iu \, \theta \, \hat{H} \right) @f$ exactly
 *   when all PauliStr in @p sum = @f$ \hat{H} @f$ commute. 
 * 
 * @param[in,out] qureg  the state to modify.
 * @param[in]     sum    a weighted sum of Pauli strings to approximately exponentiate.
 * @param[in]     angle  an effective prefactor of @p sum in the exponent.
 * @param[in]     order  the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...)
 * @param[in]     reps   the number of Trotter repetitions
 * 
 * @throws @validationerror
 * - if @p qureg or @p sum are uninitialised.
 * - if @p sum contains non-identities on qubits beyond the size of @p qureg.
 * - if @p order is not 1 nor a positive, @b even integer.
 * - if @p reps is not a positive integer.
 * 
 * @author Tyson Jones
 */
void applyNonUnitaryTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiControlledTrotterizedPauliStrSumGadget()
void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, std::vector<int> controls, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyMultiStateControlledTrotterizedPauliStrSumGadget()
void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, std::vector<int> controls, std::vector<int> states, PauliStrSum sum, qreal angle, int order, int reps);


#endif // __cplusplus

/** @} */



#endif // TROTTERISATION_H

/** @} */ // (end file-wide doxygen defgroup)
