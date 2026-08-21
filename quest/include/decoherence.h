/** @file
 * API signatures for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 * 
 * @author Tyson Jones
 * 
 * @defgroup decoherence Decoherence
 * @ingroup api
 * @brief Functions for effecting decoherence channels upon density matrices.
 * @{
 */

#ifndef DECOHERENCE_H
#define DECOHERENCE_H

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/channels.h"



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


/** Applies a one-qubit dephasing channel upon the density matrix @p qureg,
 * where @p prob is the probability of a @f$ \hat{Z} @f$ error upon the 
 * @p target qubit.
 * 
 * This is also known as a phase-flip channel.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob and @f$ t = @f$ @p target. 
 * 
 * This function effects
 * @f[ 
        \dmrho 
            \;\rightarrow\;
        (1 - p) \, \dmrho 
            \,+\, 
        p \, \hat{Z}_t \,\dmrho\, \hat{Z}_t.
 * @f]
 *
 * This is a physically valid operation (is completely positive and trace preserving)
 * when @f$ 0 \le p \le 1 @f$, and is a meaningful noise channel (i.e. induces mixing)
 * for @f$ 0 < p \le 1/2 @f$.
 *
 * @constraints
 *
 * - Parameter @p prob must be a valid probability and ergo satisfy @f$ 0 \le p \le 1 @f$
 *   unless validation is disabled via setQuESTValidationOff(), which permits quasi-probability
 *   channels with, for example, negative probabilities.
 * - The maximum permitted error probability is @f$ p = 1/2 @f$ (unless validation is disabled)
 *   at which the qubit becomes completely "dephased", and the off-diagonal "coherences" become zero.
 * - With validation disabled, the channel remains CPTP for @f$ 1/2 \le p \le 1 @f$.
 * 
 * @equivalences
 * 
 * This function is equivalent to (but much faster than):
 * - mixPaulis() with a zero probability for the @f$\hat{X}@f$ and @f$\hat{Y}@f$ components.
 *   ```
    mixPaulis(qureg, target, 0, 0, prob);
 *   ```
 * - mixKrausMap() with (scaled) @f$\hat{\id}@f$ and @f$\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-prob);
    qreal b = sqrt(prob);

    KrausMap map = createInlineKrausMap(1, 2, {
        {{a,0},{0, a}}, // a * I
        {{b,0},{0,-b}}  // b * Z
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 * - mixQureg() with a duplicated Qureg modified under applyPauliZ().
 *   ```
    Qureg clone = createCloneQureg(qureg);
    applyPauliZ(clone);
    mixQureg(qureg, other, prob);
 *   ```
 *
 * @notyetvalidated
 * 
 * @param[in,out] qureg  the density matrix to modify.
 * @param[in]     target the index of the target qubit.
 * @param[in]     prob   the probability of any error.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix (unless @p prob is 0).
 * - if @p target is an invalid qubit index.
 * - if @p prob < 0, or @p prob > 1/2.
 * @see
 * - mixTwoQubitDephasing()
 * - mixDepolarising()
 * - mixPaulis()
 * @author Tyson Jones
 */
void mixDephasing(Qureg qureg, int target, qreal prob);


/** Applies a two-qubit dephasing channel upon the density matrix @p qureg,
 * where @p prob is the probability of a @f$ \hat{Z} @f$ error upon either
 * or both of qubits @p target1 and @p target2, which are interchangeable.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob, @f$ t_1 = @f$ @p target1 and @f$ t_2 = @f$ @p target2.
 * 
 * This function effects
 * @f[ 
        \dmrho 
            \;\rightarrow\;
            (1 - p) \, \dmrho 
                \,+\,
            \frac{p}{3} \left(
            \hat{Z}_{t_1} \dmrho \hat{Z}_{t_1}
                \,+\,
            \hat{Z}_{t_2} \dmrho \hat{Z}_{t_2}
                \,+\,
            \hat{Z}_{t_1}  \hat{Z}_{t_2} \dmrho \hat{Z}_{t_1} \hat{Z}_{t_2}
        \right).
 * @f]
 *
 * This is a physically valid operation (is completely positive and trace preserving)
 * when @f$ 0 \le p \le 1 @f$, and is a meaningful noise channel (i.e. induces mixing)
 * for @f$ 0 < p \le 3/4 @f$.
 *
 * @constraints
 *
 * - Parameter @p prob must be a valid probability and ergo satisfy @f$ 0 \le p \le 1 @f$
 *   unless validation is disabled via setQuESTValidationOff(), which permits quasi-probability
 *   channels with, for example, negative probabilities.
 * - The maximum permitted error probability is @f$ p = 3/4 @f$ (unless validation is disabled),
 *   at which the channel is "maximum strength" and the off-diagonal "coherences" of the two qubits
 *   become zero.
 * - With validation disabled, the channel remains CPTP for @f$ 3/4 \le p \le 1 @f$.
 * 
 * @equivalences
 * 
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with (scaled) @f$\hat{\id}\otimes\hat{\id}@f$, @f$\hat{\id}\otimes\hat{Z}@f$, 
 *   @f$\hat{Z}\otimes\hat{\id}@f$ and @f$\hat{Z}\otimes\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-prob);
    qreal b = sqrt(prob/3);

    KrausMap map = createInlineKrausMap(2, 4, {
        {{a,0,0,0},{0, a,0,0},{0,0, a,0},{0,0,0, a}}, // a * II
        {{b,0,0,0},{0,-b,0,0},{0,0, b,0},{0,0,0,-b}}, // b * IZ
        {{b,0,0,0},{0, b,0,0},{0,0,-b,0},{0,0,0,-b}}, // b * ZI
        {{b,0,0,0},{0,-b,0,0},{0,0,-b,0},{0,0,0, b}}  // b * ZZ
    });
    
    int targets[] = {target1, target2};
    mixKrausMap(qureg, targets, 2, map);
 *   ```
 *
 * @notyetvalidated
 * 
 * @param[in,out] qureg   the density matrix to modify.
 * @param[in]     target1 the index of the first target qubit.
 * @param[in]     target2 the index of the second target qubit.
 * @param[in]     prob    the probability of any error.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix (unless @p prob is 0).
 * - if @p target1 or @p target2 is an invalid qubit index, or are the same.
 * - if @p prob < 0, or @p prob > 3/4.
 * @see
 * - createDensityQureg()
 * - mixDephasing()
 * - mixDepolarising()
 * - mixPaulis()
 * @author Tyson Jones
 */
void mixTwoQubitDephasing(Qureg qureg, int target1, int target2, qreal prob);


/** Applies a one-qubit homogeneous depolarising channel upon the density matrix @p qureg,
 * where @p prob is the probability of any error upon the @p target qubit.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob and @f$ t = @f$ @p target. 
 * 
 * This function effects
 * @f[ 
        \dmrho \;\rightarrow\; 
            (1 - p) \, \dmrho \,+\, \frac{p}{3} \left( 
            \hat{X}_t \dmrho \hat{X}_t \,+\,
            \hat{Y}_t \dmrho \hat{Y}_t \,+\,
            \hat{Z}_t \dmrho \hat{Z}_t
        \right).
 * @f]
 *
 * This is a physically valid operation (is completely positive and trace preserving)
 * when @f$ 0 \le p \le 1 @f$, and is a meaningful noise channel (i.e. induces mixing)
 * for @f$ 0 < p \le 3/4 @f$. 
 *
 * @constraints
 *
 * - Parameter @p prob must be a valid probability and ergo satisfy @f$ 0 \le p \le 1 @f$
 *   unless validation is disabled via setQuESTValidationOff(), which permits quasi-probability
 *   channels with, for example, negative probabilities.
 * - The maximum permitted error probability is @f$ p = 3/4 @f$ (unless validation is disabled),
 *   at which the channel is maximum strength, and the qubit enters the maximally mixed state.
 * - With validation disabled, the channel remains CPTP for @f$ 3/4 \le p \le 1 @f$.
 * 
 * @equivalences
 * 
 * This function is equivalent to (but much faster than):
 * - mixPaulis() with a uniform probability.
 *   ```
    mixPaulis(qureg, target, prob/3, prob/3, prob/3);
 *   ```
 * - mixKrausMap() with (scaled) @f$\hat{\id}@f$, @f$\hat{X}@f$, @f$\hat{Y}@f$ and @f$\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-prob);
    qreal b = sqrt(prob/3);

    KrausMap map = createInlineKrausMap(1, 4, {
        {{a,0},{0, a}}, // a * I
        {{0,b},{b, 0}}, // b * X
        {{b,0},{0,-b}}, // b * Z
        {{0,-1i*b},{1i*b,0}}, // b * Y
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 *
 * @notyetvalidated
 * 
 * @param[in,out] qureg   the density matrix to modify.
 * @param[in]     target  the index of the target qubit.
 * @param[in]     prob    the probability of any error.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix (unless @p prob is 0).
 * - if @p target is an invalid qubit index.
 * - if @p prob < 0, or @p prob > 3/4.
 * @see
 * - mixTwoQubitDepolarising()
 * - mixDephasing()
 * - mixPaulis()
 * @author Tyson Jones
 */
void mixDepolarising(Qureg qureg, int target, qreal prob);


/** Applies a two-qubit homogeneous depolarising channel upon the density matrix @p qureg,
 * where @p prob is the probability of any error upon either
 * or both of qubits @p target1 and @p target2, which are interchangeable.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob, @f$ t_1 = @f$ @p target1 and @f$ t_2 = @f$ @p target2.
 * 
 * This function effects:
 * @f[
        \dmrho \; \rightarrow \;
            (1 - p) \dmrho 
                +
            \frac{p}{15} \left(
                \sum_{\hat{\sigma} \in \{\hat{\id},\hat{X},\hat{Y},\hat{Z}\}}
                \sum_{\hat{\sigma}' \in \{\hat{\id},\hat{X},\hat{Y},\hat{Z}\}}
                \hat{\sigma}_{t_1} \hat{\sigma}_{t_2}'
                \; \dmrho \;
                \hat{\sigma}_{t_1} \hat{\sigma}_{t_2}'
            \right)
            - \frac{p}{15} \hat{\id}_{t_1} \hat{\id}_{t_2} \dmrho \hat{\id}_{t_1} \hat{\id}_{t_2},
 * @f]
 *
 * or verbosely:
 * 
 * @f[
    \dmrho \; \rightarrow \;
    (1 - p) \, \rho + \frac{p}{15} \; 
    \left( 
    \begin{gathered}
        \hat{X}_{t_1} \, \rho \, \hat{X}_{t_1} + 
        \hat{Y}_{t_1} \, \rho \, \hat{Y}_{t_1} + 
        \hat{Z}_{t_1} \, \rho \, \hat{Z}_{t_1} + 
            \\
        \hat{X}_{t_2} \, \rho \, \hat{X}_{t_2} + 
        \hat{Y}_{t_2} \, \rho \, \hat{Y}_{t_2} + 
        \hat{Z}_{t_2} \, \rho \, \hat{Z}_{t_2} + 
            \\
        \hat{X}_{t_1} \hat{X}_{t_2} \, \rho \, \hat{X}_{t_1} \hat{X}_{t_2} + 
        \hat{Y}_{t_1} \hat{Y}_{t_2} \, \rho \, \hat{Y}_{t_1} \hat{Y}_{t_2} + 
        \hat{Z}_{t_1} \hat{Z}_{t_2} \, \rho \, \hat{Z}_{t_1} \hat{Z}_{t_2} +
            \\
        \hat{X}_{t_1} \hat{Y}_{t_2} \, \rho \, \hat{X}_{t_1} \hat{Y}_{t_2} + 
        \hat{Y}_{t_1} \hat{Z}_{t_2} \, \rho \, \hat{Y}_{t_1} \hat{Z}_{t_2} + 
        \hat{Z}_{t_1} \hat{X}_{t_2} \, \rho \, \hat{Z}_{t_1} \hat{X}_{t_2} +
            \\
        \hat{X}_{t_1} \hat{Z}_{t_2} \, \rho \, \hat{X}_{t_1} \hat{Z}_{t_2} + 
        \hat{Y}_{t_1} \hat{X}_{t_2} \, \rho \, \hat{Y}_{t_1} \hat{X}_{t_2} + 
        \hat{Z}_{t_1} \hat{Y}_{t_2} \, \rho \, \hat{Z}_{t_1} \hat{Y}_{t_2}
    \end{gathered}
    \right).
 * @f]
 *
 * This is a physically valid operation (is completely positive and trace preserving)
 * when @f$ 0 \le p \le 1 @f$, and is a meaningful noise channel (i.e. induces mixing)
 * for @f$ 0 < p \le 15/16 @f$.
 * 
 * @constraints
 *
 * - Parameter @p prob must be a valid probability and ergo satisfy @f$ 0 \le p \le 1 @f$
 *   unless validation is disabled via setQuESTValidationOff(), which permits quasi-probability
 *   channels with, for example, negative probabilities.
 * - The maximum permitted error probability is @f$ p = 15/16 @f$ (unless validation is disabled),
 *   at which the channel is maximum strength, and the target qubits become maximally mixed.
 * - With validation disabled, the channel remains CPTP for @f$ 15/16 \le p \le 1 @f$.
 *
 * @equivalences
 * 
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with Kraus operators containing every possible tensor product 
 *   of two Pauli matrices, all scaled by @f$ (p/15)^{1/2} @f$, _except_ for
 *   @f$ \hat{\id} \otimes \hat{\id} @f$ which is scaled by @f$ (1-16p/15)^{1/2} @f$.
 *
 * @notyetvalidated
 *
 * @param[in,out] qureg   the density matrix to modify.
 * @param[in]     target1 the index of the first target qubit.
 * @param[in]     target2 the index of the second target qubit.
 * @param[in]     prob    the probability of any error.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix (unless @p prob is 0).
 * - if @p target1 or @p target2 is an invalid qubit index, or are the same.
 * - if @p prob < 0, or @p prob > 15/16.
 * @see
 * - createDensityQureg()
 * - mixTwoQubitDephasing()
 * @author Tyson Jones
 */
void mixTwoQubitDepolarising(Qureg qureg, int target1, int target2, qreal prob);


/** Applies a one-qubit amplitude damping channel upon the density matrix @p qureg,
 * where @p prob is the probability of the @p target qubit relaxing to the zero state.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob and @f$ t = @f$ @p target.
 * 
 * This function effects
 * @f[ 
        \dmrho \; \rightarrow \; 
        \hat{K}_t^{(1)} \dmrho \, {\hat{K}_t^{(1)}}^\dagger 
            \,+\,
        \hat{K}_t^{(2)} \dmrho \, {\hat{K}_t^{(2)}}^\dagger
 * @f]
 * where @f$ \hat{K}^{(1)} @f$ and @f$ \hat{K}^{(2)} @f$ are one-qubit Kraus operators
 * @f[
    \hat{K}^{(1)} = \begin{pmatrix} 1 & 0 \\ 0 & \sqrt{1-p} \end{pmatrix},
    \;\;
    \hat{K}^{(2)} = \begin{pmatrix} 0 & \sqrt{p} \\ 0 & 0 \end{pmatrix}.
 * @f]
 *
 * This is a physically valid operation (the channel is completely positive and trace
 * preserving) for @f$ 0 \le p \le 1 @f$. Note however that it may actually reduce
 * mixing and increase purity, depending on @f$ p @f$ and the qubit's initial state.
 * 
 * @constraints
 *
 * - Parameter @p prob must be a valid probability and ergo satisfy @f$ 0 \le p \le 1 @f$.
 *   Beware that disabling validation with setQuESTValidationOff() and passing @p prob
 *   outside this domain will result in mathematically erroneous results; amplitudes which
 *   are erroneously zero, or @c NaN, as output by @c sqrt().
 *
 * @equivalences
 * 
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with the above Kraus operators.
 *   ```
    KrausMap map = createInlineKrausMap(1, 2, {
        {{1,0},{0,sqrt(1-prob)}}, // K1
        {{0,sqrt(p)},{0,0}}       // K2
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 *
 * @notyetvalidated
 * 
 * @param[in,out] qureg   the density matrix to modify.
 * @param[in]     target  the index of the target qubit.
 * @param[in]     prob    the probability of relaxing to zero.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix (unless @p prob is 0).
 * - if @p target is an invalid qubit index.
 * - if @p prob < 0 or @p prob > 1.
 * @see
 * - mixKrausMap()
 * @author Tyson Jones
 */
void mixDamping(Qureg qureg, int target, qreal prob);


/** Applies a one-qubit inhomogeneous Pauli channel upon the density matrix @p qureg.
 * 
 * This is a generalisation of mixDepolarising(), permitting inhomogeneous error probabilities.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ t = @f$ @p target, and
 * @f$ p_x = @f$ @p probX, @f$ p_y = @f$ @p probY, @f$ p_z = @f$ @p probZ.
 * 
 * This function effects
 * @f[ 
        \dmrho \;\rightarrow\;
             (1 - p) \, \dmrho 
             \,+\,
            p_x \, \hat{X}_t \dmrho \hat{X}_t
             \,+\,
            p_y \, \hat{Y}_t \dmrho \hat{Y}_t
             \,+\,
            p_z \, \hat{Z}_t \dmrho \hat{Z}_t.
 * @f]
 *
 * This operation is physically valid (completely positive and trace preserving) when each
 * probability is valid (@f$ 0 \le p_i \le 1 @f$), and together satisfy
 * @f[
 *  p_x + p_y + p_z \le 1.
 * @f]
 * The operation is a meaningful noise channel (decreases purity) when the probabilities are
 * below that which induces maximal mixing; when the probability of no error is greater than
 * (or equal to) the probability of any error.
 * @f[
 *   1 - (p_x + p_y + p_z) \ge \max(p_x, p_y, p_z).
 * @f]
 * 
 * @constraints
 *
 * - Each of @p probX, @p probY, and @p probZ must be a valid probability, i.e. @f$ 0 \le p_i \le 1 @f$,
 *   and the probability of no error (one minus their sum) must also be valid. This particular validation
 *   is insensitive to the validation epsilon as controlled with setQuESTValidationEpsilon(), but can instead
 *   be relaxed with setQuESTValidationOff(), to effect channels which are not completely-positive and
 *   trace-preserving, such as quasi-probability channels.
 * - The channel strength must not exceed that which induces maximal mixing (unless validation is disabled),
 *   whereby the probability of any particular error equals that of no error, as discussed above.
 * 
 * @equivalences
 * 
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with (scaled) @f$\hat{\id}@f$, @f$\hat{X}@f$, @f$\hat{Y}@f$ and @f$\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-probX-probY-probZ);
    qreal b = sqrt(probX);
    qreal c = sqrt(probY);
    qreal d = sqrt(probZ);

    KrausMap map = createInlineKrausMap(1, 4, {
        {{a,0},{0, a}}, // a * I
        {{0,b},{b, 0}}, // b * X
        {{d,0},{0,-d}}, // d * Z
        {{0,-1i*c},{1i*c,0}}, // c * Y
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 *
 * @notyetvalidated
 * 
 * @param[in,out] qureg  the density matrix to modify.
 * @param[in]     target the index of the target qubit.
 * @param[in]     probX  the probability of an X operator upon @p target.
 * @param[in]     probY  the probability of an Y operator upon @p target.
 * @param[in]     probZ  the probability of an Z operator upon @p target.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix (unless @p probX = @p probY = @p probZ = 0).
 * - if @p target is an invalid qubit index.
 * - if any probability is invalid (below zero or above one).
 * - if the probability of any error exceeds that of no error.
 * @see
 * - mixDephasing()
 * - mixDepolarising()
 * - mixKrausMap
 * @author Tyson Jones
 */
void mixPaulis(Qureg qureg, int target, qreal probX, qreal probY, qreal probZ);


/** Modifies the density matrix @p qureg to the mixture of itself and the density matrix or
 * statevector @p other.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho_1 = @f$ @p qureg and @f$ p = @f$ @p prob.
 * 
 * - When @p other is a density matrix @f$ \dmrho_2 @f$, this function effects
 *   @f[ 
        \dmrho_1 \;\rightarrow \;
            (1 - p) \, \dmrho_1 
                \,+\,
            p \, \dmrho_2.
 *   @f]
 * - When @p other is a statevector @f$ \ket{\psi_2} @f$, this function effects
 *   @f[ 
        \dmrho_1 \;\rightarrow \;
            (1 - p) \, \dmrho_1 
                \,+\,
            p \, \ketbra{\psi_2}{\psi_2}.
 *   @f]
 *
 * @constraints
 *
 * - Parameter @p prob must be a valid probability, satisfying @f$ 0 \le p \le 1 @f$,
 *   though can be relaxed to any real scalar by disabling validation with
 *   setQuESTValidationOff().
 * - @p qureg and @p other must contain the same number of qubits.
 * - If @p other is a density matrix, it must be identically distributed to @p qureg 
 *   (although the parallelisation backends, like multithreading and GPU acceleration, are
 *    permitted to differ).
 * - If @p other is a statevector and @p qureg is not distributed, neither too must @p other.
 * 
 * @notyetvalidated
 * 
 * @param[in,out] qureg the density matrix to modify.
 * @param[in]     other the density matrix or statevector to mix into @p qureg.
 * @param[in]     prob  the coefficient of @p other in the mixture.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix.
 * - if @p qureg and @p other contain a different number of qubits.
 * - if @p prob is not a valid probability.
 * - if @p other is a density matrix which is differently distributed to @p qureg.
 * - if @p other is a distributed statevector, but @p qureg is not distributed.
 * @see
 * - mixKrausMap()
 * @author Tyson Jones
 */
void mixQureg(Qureg qureg, Qureg other, qreal prob);


/** Applies a general, any-size channel described as a Kraus map upon the density matrix @p qureg.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg, @f$ \vec{t} = @f$ @p targets and @f$ \hat{K}^{(i)} @f$
 * denote the @f$i@f$-th Kraus operator in @p map.
 * 
 * This function effects
 * @f[ 
        \dmrho \; \rightarrow \; 
        \sum\limits_i
        \hat{K}_{\vec{t}}^{(i)} \dmrho \, {\hat{K}_{\vec{t}}^{(i)}}^\dagger.
 * @f]
 *
 * The channel is completely positive and trace preserving (CPTP) when the Kraus operators satisfy
 * @f[
        \sum\limits_i {\hat{K}_{\vec{t}}^{(i)}}^\dagger \hat{K}_{\vec{t}}^{(i)} = \mathbb{1}.
 * @f]
 *
 * @constraints
 * 
 * - The number of channel targets @p numTargets must agree with the size of the Kraus map.
 * - The channel must be approximately CPTP, such that difference between
 *   @f$ \sum\limits_i {\hat{K}_{\vec{t}}^{(i)}}^\dagger \hat{K}_{\vec{t}}^{(i)} @f$ and 
 *   @f$ \mathbb{1} @f$ has no element of absolute value greater than the validation
 *   epsilon @f$ \valeps @f$. This can be adjusted with setQuESTValidationEpsilon(), and
 *   relaxed entirely by setting @f$ \valeps = 0 @f$.
 * - When @p qureg is distributed, each node must contain at least @c pow(2,2*numTargets)
 *   many amplitudes, to ensure sufficient communication buffers are allocated.
 * 
 * @equivalences
 * 
 * This function calls mixSuperOp(), passing the corresponding superoperator of @p map, which has the form
 * @f[
        \hat{S} = \sum\limits_i {\hat{K}_{\vec{t}}^{(i)}}^* \otimes \hat{K}_{\vec{t}}^{(i)}.
 * @f]
 *
 * @notyetvalidated
 * 
 * @param[in,out] qureg      the density matrix to modify.
 * @param[in]     targets    the list of target qubit indices.
 * @param[in]     numTargets the length of @p targets
 * @param[in]     map        a compatible-sized KrausMap.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix.
 * - if @p targets contains a duplicate or invalid qubit index.
 * - if @p numTargets is less than one, or exceeds the size of @p qureg.
 * - if @p map is not initialised.
 * - if @p map contains a different number of qubits than @p numTargets.
 * - if @p map is not CPTP.
 * - if @p qureg is distributed and @p numTargets exceeds the number of
 *   qubits in @p qureg minus half log-2 of the number of processes.
 * @see
 * - createKrausMap()
 * - [createInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__create.html#gae9c49a6443896ef590ff1e4cfaa4912b)
 * - setKrausMap()
 * @author Tyson Jones
 */
void mixKrausMap(Qureg qureg, int* targets, int numTargets, KrausMap map);


/** Applies a superoperator upon the linearised density matrix @p qureg, where @p targets
 * span the ket space.
 * 
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg contain @f$N@f$ qubits, with amplitudes @f$ \alpha_{ij} @f$.
 * @f[
        \dmrho = \sum\limits_i^{2^N} \sum\limits_j^{2^N} \alpha_{ij} \ket{i}\bra{j}.
 * @f] 
 * Internally, this matrix of dimension @f$ 2^N \times 2^N @f$ is stored in a vectorised
 * form @f$ \ket{\rho} @f$ of dimension @f$ 2^{2N} \times 1 @f$, which concatenates the columns
 * of @f$ \dmrho @f$.
 * @f[  
        \begin{aligned}
        \ket{\rho} &= \sum\limits_i^{2^N} \sum\limits_j^{2^N} \alpha_{ij} \ket{j} \ket{i} \\
                   &= \sum\limits_k^{2^{2N}} \beta_k \ket{k}
        \end{aligned}
 * @f]
 * This resembles an unnormalised statevector of twice as many qubits as @f$ \dmrho @f$, whereby
 * operators are left- and right-multiplied as
 * @f[
        \ket{ \hat{A} \, \rho \, \hat{B} } = \hat{B}^T\otimes \hat{A} \ket{\rho}.
 * @f]
 * 
 * Let @f$ \vec{t} = @f$ @p targets, @f$ n = @f$ @p numTargets, and let @f$ \hat{S} = @f$ superop.
 * The @f$n@f$-qubit superoperator @f$\hat{S}@f$ is a @f$ 2^{2n} \times 2^{2n} @f$ complex matrix
 * which operates upon both the ket and bra partitions of the linearised density matrix @f$\ket{\rho}@f$.
 * 
 * The targets @f$\vec{t}@f$ are treated as the ket qubits, specified in order of increasing
 * significance, where the first qubit corresponds to the rightmost partition of the matrix
 * @f$\hat{S}@f$. Concretely, let @f$\vec{t}+N@f$ notate the list of indices obtained by adding @f$N@f$
 * to every element of @f$\vec{t}@f$, and let @f$(\vec{t} \cup \vec{t}+N)@f$ the result of concatenating
 * this new list with @f$\vec{t}@f$. Then, this function effects
 * @f[
        \ket{\rho} \rightarrow \hat{S}_{(\vec{t} \cup \vec{t}+N)} \ket{\rho},
 * @f]
 * which is mathematically identical to left-applying the @f$2n@f$-qubit matrix @f$\hat{S}@f$ upon a
 * @f$2N@f$-qubit statevector @f$\ket{\rho}@f$.
 * 
 * See mixKrausMap() for an example of the construction of @f$\hat{S}@f$.
 * 
 * @constraints
 * 
 * - There is no requirement nor validation that @p superop is CPTP, and so is permitted to
 *   break state normalisation and interpretability.
 * - The number of targets must agree with the number of qubits upon which @p superop acts.
 * - When @p qureg is distributed, each node must contain at least @c pow(2,2*numTargets)
 *   many amplitudes, to ensure sufficient communication buffers are allocated.
 * 
 * @equivalences
 * 
 * This function is equivalent to calling leftapplyCompMatr() upon @p qureg, having prepared
 * @p superop as a CompMatr, and passing a target list prepared as  @f$(\vec{t} \cup \vec{t}+N)@f$ above.
 * 
 * @notyetvalidated
 * 
 * @param[in,out] qureg      the density matrix to modify.
 * @param[in]     targets    the list of target qubit indices.
 * @param[in]     numTargets the length of @p targets
 * @param[in]     superop    a compatible-sized SuperOp.
 * @throws @validationerror
 * - if @p qureg is not initialised.
 * - if @p qureg is not a density matrix.
 * - if @p targets contains a duplicate or invalid qubit index.
 * - if @p numTargets is less than one, or exceeds the size of @p qureg.
 * - if @p superop is not initialised, or not sync'ed (e.g. via syncSuperOp()).
 * - if @p superop contains a different number of qubits than @p numTargets.
 * - if @p qureg is distributed and @p numTargets exceeds the number of
 *   qubits in @p qureg minus half log-2 of the number of processes.
 * @see
 * - createSuperOp()
 * - createInlineSuperOp()
 * - setSuperOp()
 * - mixKrausMap()
 * @author Tyson Jones
 */
void mixSuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop);


// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). These
 * are included in the file-wide doxygen group (no subgroups).
 */

#ifdef __cplusplus

#include <vector>

/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see mixKrausMap()
void mixKrausMap(Qureg qureg, std::vector<int> targets, KrausMap map);

/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see mixSuperOp()
void mixSuperOp(Qureg qureg, std::vector<int> targets, SuperOp superop);

#endif // __cplusplus



#endif // DECOHERENCE_H

/** @} */ // (end file-wide doxygen defgroup)
