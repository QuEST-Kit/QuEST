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
 * Effects an approximation to the exponential of @p sum, weighted by @p angle times @f$ i @f$, upon @p qureg,
 * via the symmetrized Trotter-Suzuki decomposition (<a href="https://arxiv.org/abs/math-ph/0506007">arXiv</a>).
 * Increasing @p reps (the number of Trotter repetitions) or @p order (an even, positive integer or one) 
 * improves the accuracy of the approximation by reducing the "Trotter error" due to non-commuting 
 * terms of @p sum, though increases the runtime linearly and exponentially respectively.
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p sum and @f$ \theta = @f$ @p angle @f$ \in \mathbb{R} @f$. This function approximates 
 * the action of
 * @f[
      \exp \left(\iu \, \theta \, \hat{H} \right)
 * @f]
 * via a Trotter-Suzuki decomposition of the specified @p order and number of repetitions (@p reps).
 * Simulation is exact, regardless of @p order or @p reps, only when all terms in @p sum commute.
 * 
 * @important
 *   Observe that @f$ \theta @f$ lacks the @f$ -\frac{1}{2} @f$ prefactor present in other functions like
 *   applyPauliGadget().
 * 
 * To be precise, let @f$ r = @f$ @p reps and assume @p sum is composed of
 * @f$ T @f$-many terms of the form
 * @f[
      \hat{H} = \sum\limits_j^T c_j \, \hat{\sigma}_j
 * @f]
 * where @f$ c_j @f$ is the coefficient of the @f$ j @f$-th PauliStr @f$ \hat{\sigma}_j @f$.
 * 
 * - When @p order=1, this function performs first-order Trotterisation, where the terms of @p sum
 *   are effected in a repeated, arbitrary but fixed order.
 *   @f[
       \exp(\iu \, \theta \, \hat{H} )
          \approx 
        \prod\limits^{r} 
        \prod\limits_{j=1}^{T} 
        \exp \left( \iu \, \frac{\theta \, c_j}{r} \, \hat\sigma_j \right).
 *   @f]
 *
 * - When @p order=2, this function performs the lowest order "symmetrized" Suzuki decomposition, whereby
 *   each repetition effects the terms of @p sum forward then in reverse.
 *   @f[
       \exp(\iu \, \theta \, \hat{H} )
          \approx 
        \prod\limits^{r} \left[
             \prod\limits_{j=1}^{T} \exp \left( \iu \frac{\theta \, c_j}{2 \, r}  \hat\sigma_j \right)
              \prod\limits_{j=T}^{1} \exp \left( \iu \frac{\theta \, c_j}{2 \, r}  \hat\sigma_j \right)
         \right].
 *   @f]
 *
 * - Greater, even values of @p order (denoted by symbol @f$ n @f$) invoke higher-order symmetrized decompositions 
 *   @f$ S[\theta,n,r] @f$. These see the lower order Trotter circuits repeated twice forward, then reversed, then 
 *   twice forward again, recursively. To be precise, letting @f$ p = \left( 4 - 4^{1/(n-1)} \right)^{-1} @f$, these
 *   satisfy
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
 * - By passing @f$ \theta = - \Delta t / \hbar @f$, this function approximates unitary time evolution of a closed 
 *   system under the time-independent Hamiltonian @p sum = @f$ \hat{H} @f$ over a duration of @f$ \Delta t @f$, as
 *   described by propagator
 *   @f[
        \hat{U}(\Delta t) = \exp(- \iu \, \Delta t  \,\hat{H} \, / \, \hbar),
 *   @f]
 *   as utilised by the function applyTrotterizedUnitaryTimeEvolution().
 * 
 * - This function is equivalent to applyTrotterizedNonUnitaryPauliStrSumGadget() when passing
 *   a @p qcomp instance with a zero imaginary component as the @p angle parameter. This latter 
 *   function is useful for generalising dynamical simulation to imaginary-time evolution.
 * 
 * @constraints
 * 
 * - Unitarity of the prescribed exponential(s) requires that @p sum is Hermitian, ergo containing
 *   only real coefficients. Validation will check that @p sum is approximately Hermitian, permitting
 *   coefficients with imaginary components smaller (in magnitude) than epsilon.
 *   @f[ 
        \max\limits_{i} |c_i| \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setValidationEpsilon().
 *   Otherwise, use applyTrotterizedNonUnitaryPauliStrSumGadget() to permit non-Hermitian @p sum
 *   and ergo effect a non-unitary exponential(s). 
 * 
 * - The @p angle parameter is necessarily real to retain unitarity, but can be relaxed to an arbitrary 
 *   complex scalar (i.e. a @p qcomp) using applyTrotterizedNonUnitaryPauliStrSumGadget(). This permits
 *   cancelling the complex unit @f$ i @f$ to effect non-unitary @f$ \exp(\theta \, \hat{H}) @f$ as
 *   is useful for imaginary-time evolution.
 * 
 * - This function only ever effects @f$ \exp \left(\iu \, \theta \, \hat{H} \right) @f$ exactly
 *   when all PauliStr in @p sum = @f$ \hat{H} @f$ commute, or @p reps @f$ \rightarrow \infty @f$.
 * 
 * @param[in,out] qureg  the state to modify.
 * @param[in]     sum    a weighted sum of Pauli strings to approximately exponentiate.
 * @param[in]     angle  the prefactor of @p sum times @f$ i @f$ in the exponent.
 * @param[in]     order  the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...).
 * @param[in]     reps   the number of Trotter repetitions.
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
 *  - applyTrotterizedNonUnitaryPauliStrSumGadget()
 *  - applyTrotterizedUnitaryTimeEvolution()
 * 
 * @author Tyson Jones
 */
void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyetdoced
/// @notyettested
/// @see
///  - applyTrotterizedPauliStrSumGadget()
///  - applyControlledCompMatr1()
void applyTrotterizedControlledPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyetdoced
/// @notyettested
/// @see
///  - applyTrotterizedPauliStrSumGadget()
///  - applyMultiControlledCompMatr1()
void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyetdoced
/// @notyettested
/// @see
///  - applyTrotterizedPauliStrSumGadget()
///  - applyMultiStateControlledCompMatr1()
void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps);


/** @notyettested
 * 
 * A generalisation of applyTrotterizedPauliStrSumGadget() which accepts a complex @p angle and permits
 * @p sum to be non-Hermitian, thereby effecting a potentially non-unitary and non-CPTP operation.
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p sum and @f$ \theta = @f$ @p angle @f$ \in \mathbb{C} @f$. This function 
 * approximates the action of
 * @f[
      \exp \left(\iu \, \theta \, \hat{H} \right)
 * @f]
 * via a Trotter-Suzuki decomposition of the specified @p order and number of repetitions (@p reps). 
 * 
 * > See applyTrotterizedPauliStrSumGadget() for more information about the decomposition.
 *
 * @equivalences
 * 
 * - When @p angle is set to @f$ \theta = \iu \, \Delta \tau @f$ and @p sum = @f$ \hat{H} @f$ is Hermitian,
 *  this function (approximately) evolves @p qureg in imaginary-time for duration @f$ \Delta \tau @f$,
 *  effecting non-unitary propagator
    @f[
        \exp(- \Delta \tau \hat{H})
 *  @f]
 *  as utilised by applyTrotterizedImaginaryTimeEvolution().
 * 
 * - When @p angle is real and @p sum is Hermitian (i.e. has approximately real coefficients), the effected
 *   operation is unitary and this function becomes equivalent to applyTrotterizedPauliStrSumGadget().
 * 
 * @constraints
 * 
 * - This function only ever effects @f$ \exp \left(\iu \, \theta \, \hat{H} \right) @f$ exactly
 *   when all PauliStr in @p sum = @f$ \hat{H} @f$ commute. 
 * 
 * @param[in,out] qureg  the state to modify.
 * @param[in]     sum    a weighted sum of Pauli strings to approximately exponentiate.
 * @param[in]     angle  an effective prefactor of @p sum in the exponent.
 * @param[in]     order  the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...).
 * @param[in]     reps   the number of Trotter repetitions.
 * 
 * @throws @validationerror
 * - if @p qureg or @p sum are uninitialised.
 * - if @p sum contains non-identities on qubits beyond the size of @p qureg.
 * - if @p order is not 1 nor a positive, @b even integer.
 * - if @p reps is not a positive integer.
 * 
 * @author Tyson Jones
 */
void applyTrotterizedNonUnitaryPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps);


// end de-mangler
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyTrotterizedMultiControlledPauliStrSumGadget()
void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, std::vector<int> controls, PauliStrSum sum, qreal angle, int order, int reps);


/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cppvectoroverload
/// @see applyTrotterizedMultiStateControlledPauliStrSumGadget()
void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, std::vector<int> controls, std::vector<int> states, PauliStrSum sum, qreal angle, int order, int reps);


#endif // __cplusplus

/** @} */



/** 
 * @defgroup trotter_timeevol Time evolution
 * @brief Functions for approximate dynamical simulation.
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** @notyettested
 * 
 * Unitarily time evolves @p qureg for the duration @p time under the time-independent Hamiltonian @p hamil, 
 * as approximated by symmetrized Trotterisation of the specified @p order and number of cycles @p reps. 
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p hamil and @f$ t = @f$ @p time @f$ \in \mathbb{R} @f$. This function approximates 
 * the action of the unitary-time evolution operator/propagator
 * @f[
      \hat{U}(t) = \exp \left(- \iu \, t \, \hat{H} \right),
 * @f]
 * as solves the time-independent Schrödinger equation. When @p qureg is a statevector @f$ \svpsi @f$, the 
 * resulting state approximates
 * @f[
      \approx U(t) \svpsi
 * @f]
 * while when @p qureg is a density matrix @f$ \dmrho @f$, the result approximates
 * @f[
      \approx U(t) \, \dmrho \, U(t)^\dagger.
 * @f]
 *
 * > See applyTrotterizedPauliStrSumGadget() for information about the Trotter method.
 * 
 * @equivalences
 * 
 * - This function merely wraps applyTrotterizedPauliStrSumGadget() which effects @f$ \exp(\iu \theta \hat{H}) @f$,
 *   passing @f$ \theta = - t @f$.
 * 
 * @constraints
 * 
 * - Unitarity requires that @p hamil is Hermitian and ergo contains only real coefficients. Validation will check that 
 *   @p hamil is approximately Hermitian, permitting coefficients with imaginary components smaller (in magnitude) than 
 *   epsilon.
 *   @f[ 
        \max\limits_{i} |c_i| \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setValidationEpsilon(). The imaginary components
 *   of the Hamiltonian _are_ considered during simulation.
 * 
 * - The @p time parameter is necessarily real to retain unitarity. It can be substituted for a strictly imaginary
 *   scalar to perform imaginary-time evolution (as per Wick rotation @f$ t \rightarrow - \iu \tau @f$) via 
 *   applyTrotterizedImaginaryTimeEvolution(), or generalised to an arbitrary complex number through direct use of 
 *   applyTrotterizedNonUnitaryPauliStrSumGadget().
 * 
 * - The simulated system is _closed_ with dynamics described fully by the Hamiltonian @p hamil. Open or otherwise noisy
 *   system dynamics can be simulated with applyTrotterizedNoisyTimeEvolution().
 * 
 * - Simulation is exact such that the effected operation is precisely @f$ \exp(-\iu t \hat{H}) @f$ only when 
 *   @p reps @f$ \rightarrow \infty @f$ or all terms in @p hamil commute with one another. Conveniently, Trotter error
 *   does _not_ break normalisation of the state since the approximating circuit remains unitary.
 * 
 * @myexample
 * 
 *   ```
     Qureg qureg = createDensityQureg(10);
     PauliStrSum hamil =  createInlinePauliStrSum(R"(
         1   ZZI
         2   IZZ
         3   ZIZ
         1.5 XII
         2.5 IXI
         3.5 IIX
     )");

     qreal time = 0.8 * hbar;
     int order = 4;
     int reps = 20;
     applyTrotterizedUnitaryTimeEvolution(qureg, hamil, time, order, reps);
 *   ```
 *
 * @see
 *  - applyTrotterizedImaginaryTimeEvolution()
 *  - applyTrotterizedNoisyTimeEvolution()
 *  - applyTrotterizedNonUnitaryPauliStrSumGadget()
 * 
 * @param[in,out] qureg  the state to modify.
 * @param[in]     hamil  the Hamiltonian as a a weighted sum of Pauli strings.
 * @param[in]     time   the duration over which to simulate evolution.
 * @param[in]     order  the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...).
 * @param[in]     reps   the number of Trotter repetitions.
 * 
 * @throws @validationerror
 * - if @p qureg or @p hamil are uninitialised.
 * - if @p hamil contains non-identities on qubits beyond the size of @p qureg.
 * - if @p hamil is not approximately Hermitian.
 * - if @p order is not 1 nor a positive, @b even integer.
 * - if @p reps is not a positive integer.
 * 
 * @author Tyson Jones
 */
void applyTrotterizedUnitaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal time, int order, int reps);


/** @notyettested
 * 
 * Simulates imaginary-time evolution of @p qureg for the duration @p tau under the time-independent 
 * Hamiltonian @p hamil, as approximated by symmetrized Trotterisation of the specified @p order and
 * number of cycles @p reps. 
 * 
 * > [!IMPORTANT]
 * > This is a non-physical operation and breaks the normalisation of state which can be restored
 * > via setQuregToRenormalized().
 * 
 * @formulae 
 * 
 * Let @f$ \hat{H} = @f$ @p hamil and @f$ \tau = @f$ @p tau @f$ \in \mathbb{R} @f$. This function 
 * approximates the action of the non-unitary imaginary-time propagator
 * @f[
      \hat{V}(\tau) = \exp \left(- \tau \, \hat{H} \right),
 * @f]
 * as prescribed by Wick rotating (substituting time @f$ t @f$ for @f$ t \rightarrow -\iu \tau @f$)
 * the time-independent Schrödinger equation. When @p qureg is a statevector @f$ \svpsi @f$, the 
 * resulting state approximates
 * @f[
      \approx V(\tau) \svpsi
 * @f]
 * while when @p qureg is a density matrix @f$ \dmrho @f$, the result approximates
 * @f[
      \approx V(\tau) \, \dmrho \, V(\tau)^\dagger.
 * @f]
 *
 * > See applyTrotterizedPauliStrSumGadget() for information about the Trotter method.
 * 
 * @par Utility
 * 
 * Imaginary-time evolution drives the system toward the (unnormalised) groundstate of the Hamiltonian.
 * Let @f$ \{ \ket{\phi_i} \} @f$ and @f$ \{ \ket{\lambda_i} \} @f$ be the eigenstates and respective
 * eigenvalues of @f$ \hat{H} @f$, which are real due to Hermiticity.
 * @f[
    \hat{H} = \sum \limits_i \lambda_i \ket{\phi_i}\bra{\phi_i},
    \;\;\;\;\; \lambda_i \in \mathbb{R}.
 * @f]
 *
 * - When @p qureg is a statevector @f$ \svpsi @f$ and can ergo be expressed in the basis of 
 *   @f$ \{ \ket{\phi_i} \} @f$ as @f$ \svpsi = \sum_i \alpha_i \ket{\phi_i} @f$, 
 *   this function approximates
 *   @f[
        \svpsi \, \rightarrow  \, \hat{V}(\tau) \svpsi =
        \sum\limits_i \alpha_i \exp(- \tau \, \lambda_i) \ket{\phi_i}.
 *   @f]
 * - When @p qureg is a density matrix and is ergo expressible as
 *   @f$ \dmrho = \sum\limits_{ij} \alpha_{ij} \ket{\phi_i}\bra{\phi_j} @f$, this function effects
 *   @f[
        \dmrho \, \rightarrow \, \hat{V}(\tau) \dmrho \hat{V}(\tau)^\dagger =
        \sum\limits_{ij} \alpha_{ij} \exp(-\tau (\lambda_i + \lambda_j)) \ket{\phi_i}\bra{\phi_j}.
 *   @f]
 *
 * As @f$ \tau \rightarrow \infty @f$, the resulting unnormalised state approaches statevector
 * @f$ \svpsi \rightarrow \alpha_0 \exp(-\tau \lambda_0) \ket{\phi_0} @f$ or density matrix
 * @f$ \dmrho \rightarrow \alpha_{0,0} \exp(-2 \tau \lambda_0) \ket{\phi_0}\bra{\phi_0} @f$,
 * where @f$ \lambda_0 @f$ is the minimum eigenvalue and @f$ \ket{\phi_0} @f$ is the groundstate.
 * Assuming the initial overlap @f$ \alpha_0 @f$ is not zero (or exponentially tiny), 
 * subsequent renormalisation via setQuregToRenormalized() produces the pure 
 * ground-state @f$ \ket{\phi_0} @f$ or @f$ \ket{\phi_0}\bra{\phi_0} @f$.
 * 
 * Note degenerate minimum eigenvalues will yield a pure superposition of the corresponding 
 * eigenstates, with coefficients informed by the initial, relative populations.
 * 
 * @equivalences
 * 
 * - This function merely wraps applyTrotterizedNonUnitaryPauliStrSumGadget() which effects @f$ \exp(\iu \theta \hat{H}) @f$,
 *   passing @f$ \theta = \tau \iu @f$.
 * 
 * @constraints
 * 
 * - While the process of imaginary-time evolution is non-unitary (and non-physical), Hermiticity of @p hamil is still
 *   assumed, requiring it contains only real coefficients. Validation will check that @p hamil is _approximately_ Hermitian, 
 *   permitting coefficients with imaginary components smaller (in magnitude) than epsilon.
 *   @f[ 
        \max\limits_{i} |c_i| \le \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setValidationEpsilon(). Beware however that 
 *   imaginary-time evolution under a non-Hermitian Hamiltonian will _not_ necessarily approach the lowest lying eigenstate
 *   (the eigenvalues may be non-real) so is likely of limited utility.
 * 
 * - The @p tau parameter is necessarily real such that evolution approaches the groundstate (modulo renormalisation).
 *   It can generalised to an arbitrary complex number through direct use of applyTrotterizedNonUnitaryPauliStrSumGadget().
 * 
 * - Simulation is exact such that the effected operation is precisely @f$ \exp(-\tau \hat{H}) @f$ only when 
 *   @p reps @f$ \rightarrow \infty @f$ or all terms in @p hamil commute with one another.
 * 
 * @myexample
 *
 * ```
   // pray for a non-zero initial overlap
   initRandomPureState(qureg); // works even for density matrices

   // minimize then renormalise
   qreal tau = 10; // impatient infinity
   int order = 4;
   int reps = 100;
   applyTrotterizedImaginaryTimeEvolution(qureg, hamil, tau, order, reps);
   setQuregToRenormalized(qureg);

   // ground-state (phi_0)
   reportQureg(qureg);

   // lowest lying eigenvalue (lambda_0)
   qreal expec = calcExpecPauliStrSum(qureg, hamil);
   reportScalar("expec", expec);
 * ```
 *
 * @see
 *  - applyTrotterizedUnitaryTimeEvolution()
 *  - applyTrotterizedNonUnitaryPauliStrSumGadget()
 * 
 * @param[in,out] qureg  the state to modify.
 * @param[in]     hamil  the Hamiltonian as a a weighted sum of Pauli strings.
 * @param[in]     tau    the duration over which to simulate imaginary-time evolution.
 * @param[in]     order  the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...).
 * @param[in]     reps   the number of Trotter repetitions.
 * 
 * @throws @validationerror
 * - if @p qureg or @p hamil are uninitialised.
 * - if @p hamil contains non-identities on qubits beyond the size of @p qureg.
 * - if @p hamil is not approximately Hermitian.
 * - if @p order is not 1 nor a positive, @b even integer.
 * - if @p reps is not a positive integer.
 * 
 * @author Tyson Jones
 */
void applyTrotterizedImaginaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal tau, int order, int reps);


/** @notyettested
 * 
 * Simulates open dynamics of @p qureg as per the Lindblad master equation, under the time-independent
 * Hamiltonian @p hamil and jump operators @p jumps with corresponding damping rates @p damps, with 
 * evolution approximated by symmetrized Trotterisation of the specified @p order and number of cycles
 * @p reps.
 * 
 * @formulae 
 * 
 * Let @f$ \rho = @f$ @p qureg, @f$ \hat{H} = @f$ @p hamil, @f$ t = @f$ @p time, and denote the @f$ i @f$-th
 * element of @p damps and @p jumps as @f$ \gamma_i @f$ and @f$ \hat{J}_i @f$ respectively. The Lindblad
 * master equation prescribes that @f$ \rho @f$ time-evolves according to
 * @f[
     \frac{\mathrm{d}}{\mathrm{d}t} \rho = -\iu [\hat{H}, \rho] + \sum\limits_i \gamma_i \left(
          \hat{J}_i \rho \hat{J}_i^\dagger - \frac{1}{2} \left\{ \hat{J}_i^\dagger \hat{J}_i, \rho \right\}
     \right).
 * @f]
 * This function works by building a superoperator of the right-hand-side which acts upon the space of
 * linearised @f$\rho@f$,
 * @f[
     \boldsymbol{L} = -\iu \left( \hat{\id} \otimes \hat{H} - \hat{H}^* \otimes \hat{\id} \right) +
          \sum\limits_i \gamma_i \left(
               \hat{J}_i^* \otimes \hat{J}_i - \frac{1}{2} \hat{\id} \otimes (\hat{J}^\dagger J_i)
               - \frac{1}{2} (\hat{J}^\dagger J_i)^* \otimes \hat{\id}
          \right),
 * @f]
 * as a non-Hermitian weighted sum of Pauli strings (a PauliStrSum). The superoperator @f$ \boldsymbol{L} @f$
 * informs a superpropagator which exactly solves evolution as:
 * @f[
     \ket{\rho(t)} = \exp\left( t \boldsymbol{L} \right) \ket{\rho(0)}.
 * @f]
 * This function approximates the superpropagator @f$ \exp\left( t \boldsymbol{L} \right) @f$ using a higher-order 
 * symmetrized Suzuki-Trotter decomposition, as informed by parameters @p order and @p reps.
 * 
 * > See applyTrotterizedPauliStrSumGadget() for information about the Trotter method.
 * 
 * @par Utility
 * 
 * This function simulates time evolution of an open system, where the jump operators model interactions with
 * the environment. This can capture sophisticated decoherence processes of the quantum state which are untenable
 * to model as discrete operations with functions like mixKrausMap(). This function also proves useful for
 * preparing realistic, physical input states to quantum metrological circuits, or the general high-performance
 * simulation of digital time evolution of condensed matter systems.
 *
 * @equivalences
 * 
 * - When `numJumps = 0`, evolution is unitary and the Lindblad master equation simplifes to the Liouville–von Neumann 
 *   equation, which is equivalently (and more efficiently) simulated via applyTrotterizedUnitaryTimeEvolution().
 * 
 * @constraints
 * 
 * - Each damping rate in @p damps is expected to be a zero or positive number, in order for evolution to be trace 
 *   preserving. Validation will assert that each damping rate @f$ \gamma_i @f$ satisfies
 *   @f[
          \min\limits_{i} \gamma_i \ge - \valeps
 *   @f]
 *   where the validation epsilon @f$ \valeps @f$ can be adjusted with setValidationEpsilon(). Non-trace-preserving,
 *   negative damping rates can be simulated by disabling numerical validation via `setValidationEpsilon(0)`.
 * 
 * - The @p time parameter is necessarily real, and cannot be generalised to imaginary or complex like in other
 *   functions. Generalisation is trivially numerically possible, but has no established physical meaning and so
 *   is not exposed in the API. Please open an issue on Github for advice on complex-time simulation.
 * 
 * - Simulation is exact only when @p reps @f$ \rightarrow \infty @f$ or all terms in the superoperator 
 *   @f$ \boldsymbol{L} @f$ incidentally commute with one another, and otherwise incorporates Trotter error.
 *   Unlike for unitary evolution, Trotter error _does_ break normalisation of the state and so this function
 *   is generally non-trace-preserving. In theory, normalisation can be restored with setQuregToRenormalized()
 *   though noticable norm-breaking indicates evolution was inaccurate, and should instead be repeated with 
 *   increased @p order or @p reps parameters.
 * 
 * - The function instantiates superoperator @f$ \boldsymbol{L} @f$ above as a temporary PauliStrSum, incurring a 
 *   memory and time overhead which grows quadratically with the number of terms in @p hamil, plus quadratically
 *   with the number in each jump operator. These overheads may prove prohibitively costly for PauliStrSum
 *   containing very many terms.
 * 
 * @myexample
 *
 * ```
    // |+><+|
    Qureg qureg = createDensityQureg(3);
    initPlusState(qureg);

    PauliStrSum hamil = createInlinePauliStrSum(R"(
        1  IIX
        2  IYI
        3  ZZZ
    )");

    // |0><0|
    PauliStrSum jump1 = createInlinePauliStrSum(R"(
        0.5  I
        0.5  Z
    )");

    // |1><0|
    PauliStrSum jump2 = createInlinePauliStrSum(R"(
         0.5  X
        -0.5i Y
    )");

    // "noisiness"
    qreal damps[] = {.3, .4};
    PauliStrSum jumps[] = {jump1, jump2};
    int numJumps = 2;

    reportScalar("initial energy", calcExpecPauliStrSum(qureg, hamil));

    // time and accuracy
    qreal time = 0.5;
    int order = 4;
    int reps = 100;
    applyTrotterizedNoisyTimeEvolution(qureg, hamil, damps, jumps, numJumps, time, order, reps);

    reportScalar("final energy", calcExpecPauliStrSum(qureg, hamil));
 * ```
 * 
 * @see
 *  - applyTrotterizedUnitaryTimeEvolution()
 *  - applyTrotterizedImaginaryTimeEvolution()
 * 
 * @param[in,out] qureg     the density-matrix state to evolve and modify.
 * @param[in]     hamil     the Hamiltonian of the qubit system (excludes any environment).
 * @param[in]     damps     the damping rates of each jump operator in @p jumps.
 * @param[in]     jumps     the jump operators specified as PauliStrSum.
 * @param[in]     numJumps  the length of list @p jumps (and @p damps).
 * @param[in]     time      the duration through which to evolve the state.
 * @param[in]     order     the order of the Trotter-Suzuki decomposition (e.g. @p 1, @p 2, @p 4, ...).
 * @param[in]     reps      the number of Trotter repetitions.
 * 
 * @throws @validationerror
 * - if @p qureg, @p hamil or any element of @p jumps are uninitialised.
 * - if @p qureg is not a density matrix.
 * - if @p hamil or any element of @p jumps contains non-identities on qubits beyond the size of @p qureg.
 * - if @p hamil is not approximately Hermitian.
 * - if @p numJumps is negative.
 * - if any element of @p damps is not approximately positive.
 * - if the total number of Lindbladian superoperator terms overflows the `qindex` type.
 * - if all Lindbladian superoperator terms cannot simultaneously fit into CPU memory.
 * - if memory allocation of the Lindbladian superoperator terms unexpectedly fails.
 * - if @p order is not 1 nor a positive, @b even integer.
 * - if @p reps is not a positive integer.
 * 
 * @author Tyson Jones
 */
void applyTrotterizedNoisyTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal* damps, PauliStrSum* jumps, int numJumps, qreal time, int order, int reps);


// end de-mangler
#ifdef __cplusplus
}
#endif

/** @} */



#endif // TROTTERISATION_H

/** @} */ // (end file-wide doxygen defgroup)
