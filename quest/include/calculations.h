/** @file
 * API signatures for calculating properties of quantum states,
 * such as probabilities, expectation values and partial traces.
 * 
 * @author Tyson Jones
 *
 * @defgroup calculations Calculations
 * @ingroup api
 * @brief Functions for calculating properties of quantum states without modifying them.
 * @{
 */

#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"



/*
 * These signatures are divided into three partitions; those which are
 * natively C and C++ compatible (first partition), then those which are
 * only exposed to C++ (second partition) because they return 'qcomp' 
 * which cannot cross the C++-to-C ABI, and then finally the C++-only
 * convenience overloads. The first partition defines the doc groups, and 
 * the latter partition functions are added into them.
 */



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup calc_expec Expectation values
 * @brief Functions for calculating expected values of Hermitian observables.
 * @{
 */


/** Calculates the expectation value of the given Pauli string observable @p str under the given 
 * state @p qureg without modifying it. 
 * 
 * @formulae
 * 
 * Let @f$ \pstr = @f$ @p str, which notates a tensor product of single-qubit Pauli operators.
 * - When @p qureg is a statevector @f$\svpsi@f$, this function returns
 *   @f[ 
    \brapsi \pstr \svpsi \in \mathbb{R}.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function returns the real component of
 *   @f[ 
    \tr{ \pstr \dmrho }
 *   @f]
 *   which is exact when @f$\dmrho@f$ is physical (specifically Hermitian).
 * 
 * @constraints
 * 
 * - Postcondition validation will check that the calculated expectation value is approximately
 *   real (i.e. the imaginary component is smaller in size than the validation epsilon), as admitted
 *   when @p qureg is correctly normalised. This behaviour can be adjusted using setValidationEpsilon(). 
 * - Regardless of the validation epsilon, the returned value is always real and the imaginary component
 *   is discarded. The full complex value can be obtained using calcExpecNonHermitianPauliStrSum().
 * 
 * @equivalences
 * 
 * - When @p str is general, this function is equivalent to calling calcExpecPauliStrSum() with a 
 *   PauliStrSum composed of only a single PauliStr term and a unity coefficient.
 * - When @p str @f$ = \id^\otimes @f$, the output is equivalent to that of calcTotalProb().
 * 
 * @myexample
 * ```
    Qureg qureg = createQureg(10);
    initRandomPureState(qureg);

    PauliStr str = getInlinePauliStr("XYZ", {0,2,3});

    qreal expec = calcExpecPauliStr(qureg, str);
    reportScalar("expec", expec);  
 * ```
 * 
 * @see
 * - calcExpecPauliStrSum()
 * - calcExpecFullStateDiagMatr()
 * @param[in] qureg the reference state.
 * @param[in] str   the observable operator.
 * @returns The real component of the expectation value.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p str contains a (non-identity) Pauli upon a higher-index qubit than exists in @p qureg.
 * - if the output (with unreturned imaginary component) is not approximately real.
 * @notyetvalidated
 * @author Tyson Jones
 */
qreal calcExpecPauliStr(Qureg qureg, PauliStr str);


/** Calculates the expectation value of the given Hermitian observable @p sum - a weighted sum of 
 * Pauli strings - under the given state @p qureg, without modifying it. 
 * 
 * @formulae
 * 
 * Let @f$ \hat{H} = @f$ @p sum.
 * - When @p qureg is a statevector @f$\svpsi@f$, this function returns
 *   @f[ 
    \brapsi \hat{H} \svpsi \in \mathbb{R}.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function returns the real component of
 *   @f[ 
     \tr{ \hat{H} \dmrho }
 *   @f]
 *   which is the exact expectation value when @f$\dmrho@f$ is physical (or at least, Hermitian).
 * 
 * @constraints
 * 
 * - Hermiticity of @p sum requires that every coefficient within is real. 
 *   Validation will check @p sum is _approximately_ Hermitian, i.e. that
 *   @f[ 
     |\im{c}| \le \valeps
 *   @f]
 *   for all @f$c \in @f$ `sum.coeffs`. Adjust @f$\valeps@f$ using setValidationEpsilon().
 *   The sub-epsilon imaginary components of the coefficients _are_ included in calculation.
 * - Postcondition validation will check that the calculated expectation value is approximately
 *   real (i.e. the imaginary component is smaller in size than the validation epsilon), as should be
 *   admitted when @p qureg is correctly normalised, and @p sum is Hermitian.
 * - The returned value is always real, and the imaginary component is neglected even when 
 *   Hermiticity validation is relaxed and/or @p qureg is an unnormalised density matrix. 
 *   The full complex value can be obtained using calcExpecNonHermitianPauliStrSum().
 * 
 * @equivalences
 * 
 * - This function is mathematically equivalent to (albeit faster than) calling calcExpecPauliStr() upon
 *   each constituent @p PauliStr within @p sum, weighting each by its corresponding coefficient, and
 *   summing the outputs.
 * - When @p sum contains only @f$\pauliz@f$ and @f$\id@f$ operators, its corresponding operator matrix
 *   is diagonal, and could be instead effected with calcExpecFullStateDiagMatr(). This may be faster when
 *   @p sum contains very-many terms and operates upon all qubits of the register.
 * 
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    PauliStrSum sum = createInlinePauliStrSum(R"(
        0.123 XXIXX
        1.234 XYZXZ
        -1E-2 IIIII
    )");

    qreal expec = calcExpecPauliStrSum(qureg, sum);
    reportScalar("expec", expec);  
 * ```
 *
 * @param[in] qureg the reference state.
 * @param[in] sum   the observable operator.
 * @returns The real component of the expectation value.
 * @throws @validationerror
 * - if @p qureg or @p sum are uninitialised.
 * - if any PauliStr in @p sum targets a higher-index qubit than exists in @p qureg.
 * - if @p sum is not approximately Hermitian.
 * - if the output (with unreturned imaginary component) is not approximately real.
* @notyetvalidated
 * @see
 * - calcExpecNonHermitianPauliStrSum()
 * - calcExpecFullStateDiagMatr()
 * @author Tyson Jones
 */
qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum);


/** Calculates the expectation value of the given Hermitian observable @p matr - a diagonal,
 * Hermitian matrix spanning the full Hilbert space - under the given state @p qureg, without 
 * modifying it. 
 * 
 * @formulae
 * 
 * Let @f$ \hat{D} = @f$ @p matr.
 * - When @p qureg is a statevector @f$\svpsi@f$, this function returns
 *   @f[ 
    \brapsi \hat{D} \svpsi \in \mathbb{R}.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function returns the real component of
 *   @f[ 
     \tr{ \hat{D} \dmrho }
 *   @f]
 *   which is the exact expectation value when @f$\dmrho@f$ is physical (or at least, Hermitian).
 * 
 * @constraints
 * 
 * - Hermiticity of @p matr requires that every element within is real. 
 *   Validation will check @p matr is _approximately_ Hermitian, i.e. that
 *   @f[ 
     |\im{c}| \le \valeps
 *   @f]
 *   for all @f$c \in @f$ `matr.cpuElems`. Adjust @f$\valeps@f$ using setValidationEpsilon().
 * - Postcondition validation will check that the calculated expectation value is approximately
 *   real (i.e. the imaginary component is smaller in size than the validation epsilon), as should be
 *   admitted when @p qureg is correctly normalised, and @p matr is Hermitian.
 * - The returned value is always real, and the imaginary component is neglected even when @p matr
 *   Hermiticity validation is relaxed and/or @p qureg is an unnormalised density matrix. 
 *   The full complex value can be obtained using calcExpecNonHermitianFullStateDiagMatr().
 * 
 * @equivalences
 * 
 * - This function is mathematically equivalent to (albeit much faster than) calling calcExpecPauliStrSum()
 *   with a PauliStrSum consisting of all permutations of @f$\hat{I}@f$ and @f$\hat{Z}@f$ Pauli operators
 *   with a precise, linear combination of coefficients.
 * 
 * @myexample
 * 
 * ```
    Qureg qureg = createQureg(5);
    initPlusState(qureg);

    FullStateDiagMatr matr = createFullStateDiagMatr(qureg.numQubits);

    // profanely inefficient per-element initialisation
    for (int n=0; n<matr.numElems; n++) {
        qcomp elem = getQcomp(n, 0);
        setFullStateDiagMatr(matr, n, &elem, 1);
    }

    // prints "expec: 15.5"
    qreal expec = calcExpecFullStateDiagMatr(qureg, matr);
    reportScalar("expec", expec);
 * ```
 *
 * @param[in] qureg the reference state.
 * @param[in] matr  the observable operator.
 * @returns The real component of the expectation value.
 * @throws @validationerror
 * - if @p qureg or @p matr are uninitialised.
 * - if @p matr does not match the dimension of @p qureg
 * - if @p matr is distributed but @p qureg is not
 * - if @p matr is not approximately Hermitian.
 * - if the output (with unreturned imaginary component) is not approximately real.
* @notyetvalidated
 * @see
 * - calcExpecFullStateDiagMatrPower()
 * - calcExpecNonHermitianFullStateDiagMatr()
 * - calcExpecPauliStrSum()
 * @author Tyson Jones
 */
qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);


/** Calculates the expectation value of the given Hermitian observable @p matrix - a diagonal,
 * Hermitian matrix spanning the full Hilbert space - when raised to the given @p exponent,
 * under the given state @p qureg, which is not modified.
 * 
 * @formulae
 * 
 * Let @f$ \hat{D} = @f$ @p matrix and @f$x = @f$ @p exponent.
 * - When @p qureg is a statevector @f$\svpsi@f$, this function returns
 *   @f[ 
    \brapsi \hat{D}^x \svpsi \in \mathbb{R}.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function returns the real component of
 *   @f[ 
     \tr{ \hat{D}^x \dmrho }
 *   @f]
 *   which is the exact expectation value when @f$\dmrho@f$ is physical (or at least, Hermitian).
 * 
 * @constraints
 * 
 * - Hermiticity of @p matrix itself requires that every element within is real. 
 *   Validation will check @p matrix is _approximately_ Hermitian, i.e. that
 *   @f[ 
     |\im{c}| \le \valeps
 *   @f]
 *   for all @f$c \in @f$ `matr.cpuElems`. Adjust @f$\valeps@f$ using setValidationEpsilon().
 * 
 *   > [!CAUTION]
 *   > Unlike other functions (including calcExpecFullStateDiagMatr()), this function will _NOT_
 *   > consult the imaginary components of the elements of @p matrix, since a non-complex exponentiation
 *   > function is used. That is, while validation permits the imaginary components to be small, they
 *   > will be internally treated as precisely zero. This is true even when Hermiticity validation
 *   > is disabled using setValidationOff(). To consult the imaginary components of @p matrix, use
 *   > calcExpecNonHermitianFullStateDiagMatrPower().
 * 
 * - Hermiticity of @p matrix when raised to @p exponent further requires that, when @p exponent is 
 *   a non-integer, @p matrix does not contain any negative elements which would otherwise produce 
 *   complex elements in @f$\hat{D}^x@f$. This validation is always strict (i.e. independent of 
 *   @f$\valeps@f$), and demands that
 *   @f[ 
     \min(\hat{D}) \ge 0 \text{ when } x \notin \mathbb{R}.
 *   @f]
 * - Numerical stability requires that if @p exponent is negative, @p matrix does not contain any
 *   zero elements which would otherwise create divergences in @f$\hat{D}^x@f$. Validation ergo
 *   checks that when @p exponent is (strictly) negative, @p matrix contains no elements within 
 *   distance @f$\valeps@f$ to zero (regardless of the magnitude of @p exponent). Adjust
 *   @f$\valeps@f$ using setValidationEpsilon().
 * - The passed @p exponent is always real, but can be relaxed to a general complex scalar via
 *   calcExpecNonHermitianFullStateDiagMatrPower().
 * - The returned value is always real, and the imaginary component is neglected even when 
 *   Hermiticity validation is relaxed and/or @p qureg is an unnormalised density matrix. 
 *   The full complex value can be obtained using calcExpecNonHermitianFullStateDiagMatrPower().
 * 
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    initPlusState(qureg);

    FullStateDiagMatr matrix = createFullStateDiagMatr(qureg.numQubits);

    // profanely inefficient per-element initialisation
    for (int n=0; n<matrix.numElems; n++) {
        qcomp elem = getQcomp(n+1, 0);
        setFullStateDiagMatr(matrix, n, &elem, 1);
    }

    // prints "expec: 0.044503"
    qreal exponent = -2.3;
    qreal expec = calcExpecFullStateDiagMatrPower(qureg, matrix, exponent);
    reportScalar("expec", expec);
 * ```
 * @param[in] qureg     the reference state.
 * @param[in] matrix    the observable operator.
 * @param[in] exponent  the exponent to which to raise @p matrix
 * @returns The real component of the expectation value of @p matrix raised to @p exponent.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix does not match the dimension of @p qureg
 * - if @p matrix is distributed but @p qureg is not
 * - if @p matrix is not approximately Hermitian.
 * - if @p exponent is (precisely) non-integer but @p matrix contains (precisely) negative elements.
 * - if @p exponent is (precisely) negative but @p matrix contains elements which are approximately zero. 
 * - if the output (with unreturned imaginary component) is not approximately real.
 * @notyetvalidated
 * @see
 * - calcExpecNonHermitianFullStateDiagMatrPower()
 * @author Tyson Jones
 */
qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qreal exponent);


/** @} */



/** 
 * @defgroup calc_prob Probabilities
 * @brief Functions for non-destructively calculating the probabilities of measurement outcomes.
 * @{
 */


/** Calculates the probability of the full computational basis state of the specified
 * @p index. This is the probability that, when measured in the @f$ \hat{Z} @f$ basis,
 * every qubit of @p qureg is consistent with the bits of @p index.
 * 
 * Indexing is little-endian and from zero, such that (for example) computational basis state 
 * @f$ \ket{0011} @f$ (where qubits at indices @f$0@f$ and @f$1@f$ are in the @f$\ket{1}@f$ state)
 * corresponds to @p index @f$ = 3 @f$. The maximum legal @p index of an @f$N@f$-qubit
 * register is @p index @f$ = 2^N-1 @f$.
 *
 * @formulae
 * 
 * Let @f$ i = @f$ @p index.
 * 
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function returns
 *   @f[ 
      P(i) = |\braket{i}{\psi}|^2 = |\psi_i|^2
 *   @f] 
 *   where @f$\psi_i@f$ is the @f$i@f$-th amplitude of @f$\svpsi@f$.
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function returns
 *   @f[ 
      P(i) = \re{ \tr{ \ketbra{i}{i} \dmrho } } = \re{ \bra{i} \dmrho \ket{i} } = \re{ \dmrho_{ii} }
 *   @f]
 *   where @f$ \dmrho_{ii} @f$ is the @f$i@f$-th diagonal element of @f$\dmrho@f$, and is
 *   real whenever @f$ \dmrho @f$ is valid (or at least, Hermitian).
 * 
 * When @p qureg is correctly normalised, these quantities are within @f$[0, 1]@f$, and satisfy
 * @f[
      \sum\limits_{i=0}^{2^N-1} P(i) = 1
 * @f]
 * where @f$N@f$ is the number of qubits in @p qureg.
 * 
 * @equivalences
 * 
 * - This function is equivalent to obtaining the corresponding @p qureg amplitude directly
 *   and evaluating the probability.
 *   ```
     // qureg is statevector
     qcomp amp = getQuregAmp(qureg, index);
     qreal prob = pow(abs(amp, 2));

     // qureg is a density matrix
     qcomp amp = getDensityQuregAmp(qureg, index, index);
     qreal prob = real(amp);
 *   ```
 * - This function is slightly faster than, but otherwise mathematically equivalent to, invoking
 *   calcProbOfMultiQubitOutcome() and passing explicitly the bits of @p index. I.e.
 *   ```
     int qubits[qureg.numQubits];
     int outcomes[qureg.numQubits];

     for (int q=0; q<qureg.numQubits; q++) {
         qubits[q] = q;
         outcomes[q] = (index >> q) & 1;
     }

     qreal prob = calcProbOfMultiQubitOutcome(qureg, qubits, outcomes, qureg.numQubits);
 *   ```
 *   Use of calcProbOfMultiQubitOutcome() may be more convenient if only the individual qubit 
 *   outcomes are known.
 * - This function is significantly faster than, but mathematically equivalent to, preparing
 *   a secondary Qureg in the basis state @p index and computing their overlap.
 *   ```
     Qureg alt = createCloneQureg(qureg);
     initClassicalState(alt, index);
     qcomp amp = calcInnerProduct(alt, qureg);
     qreal prob = pow(abs(amp), 2);
 *   ```
 * 
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    initPlusState(qureg);

    qreal prob = calcProbOfBasisState(qureg, 2);
    reportScalar("prob of |00010>", prob);
 * ```
 *
 * @param[in] qureg the reference state, which is unchanged.
 * @param[in] index the index of the queried basis state among the ordered set of all basis states.
 * @returns The probability of the basis state at @p index.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p index is less than zero or beyond (or equal to) the dimension of @p qureg.
* @notyetvalidated
 * @see
 * - calcProbOfQubitOutcome()
 * - calcProbOfMultiQubitOutcome()
 * - getQuregAmp()
 * - getDensityQuregAmp()
 * @author Tyson Jones
 */
qreal calcProbOfBasisState(Qureg qureg, qindex index);


/** Calculates the probability of the single qubit at index @p qubit being in the
 * given computational basis @p outcome (`0` or `1`).
 *
 * @formulae
 * 
 * Let @f$ q = @f$ @p qubit and @f$ x = @f$ @p outcome, and let @f$\ketbra{x}{x}_q@f$
 * notate a projector operating upon qubit @f$ q @f$. 
 * 
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function returns
 *   @f[
      P_q(x) = \tr{ \ketbra{x}{x}_q \, \ketbra{\psi}{\psi} }
         = \sum\limits_i |\psi_i|^2 \delta_{x,i_{[q]}}
 *   @f]
 *   where @f$\psi_i@f$ is the @f$i@f$-th amplitude of @f$\svpsi@f$, and @f$i_{[q]}@f$
 *   notates the @f$q@f$-th bit of @f$i@f$.
 * - When @p qureg is a density matrix @f$ \dmrho @f$, this function returns
 *   @f[
     P_q(x) = \tr{ \ketbra{x}{x}_q \, \dmrho }
         = \sum\limits_i \re{ \dmrho_{ii} } \delta_{x,i_{[q]}}
 *   @f]
 *   where @f$ \dmrho_{ii} @f$ is the @f$i@f$-th diagonal element of @f$\dmrho@f$. This 
 *   is real whenever @f$\dmrho@f$ is validly normalised (specifically, Hermitian).
 * 
 * When @p qureg is correctly normalised, these quantities are within @f$[0, 1]@f$, and
 * satisfy
 * @f[
     P_q(x=0) + P_q(x=1) = 1.
 * @f]
 *
 * @equivalences
 * 
 * - This function is a single-qubit convenience overload of calcProbOfMultiQubitOutcome(), 
 *   which itself has optimised implementations for few-qubit outcomes.
 *   ```
     calcProbOfMultiQubitOutcome(qureg, &qubit, &outcome, 1);
 *   ```
 * - This function is much faster than, but mathematically equivalent to, summing the probability
 *   of every computational basis state (e.g. via calcProbOfBasisState()) which is consistent
 *   with the given qubit outcome.
 *   ```
     qreal prob = 0;
     qindex dim = 1 << qureg.numQubits;
     for (qindex i=0; i<dim; i++)
         if (outcome == (i >> qubit) & 1)
            prob += calcProbOfBasisState(qureg, i);
 *   ```
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    
    int qubit = 2;
    int outcome = 1;
    qreal theta = 0.3;
    applyRotateX(qureg, qubit, theta);

    // prob = cos(theta/2)^2
    qreal prob = calcProbOfQubitOutcome(qureg, qubit, outcome);
 * ```
 *
 * @param[in] qureg   the reference state, which is unchanged.
 * @param[in] qubit   the target qubit to query.
 * @param[in] outcome the outcome of @p qubit to query (i.e. `0` oe `1`).
 * @returns The probability that the given qubit is in the given outcome.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p qubit is less than zero or beyond the number of qubits in @p qureg.
 * - if @p outcome is not `0` or `1`.
* @notyetvalidated
 * @see
 * - calcProbOfMultiQubitOutcome()
 * @author Tyson Jones
 */
qreal calcProbOfQubitOutcome(Qureg qureg, int qubit, int outcome);


/** Calculates the probability that the given list of @p qubits are simultaneously in the 
 * respective single-qubit states specified in @p outcomes.
 *
 * @formulae
 * 
 * Let @f$q_j@f$ and @f$x_j@f$ notate the @f$j@f$-th qubit in @p qubits and its respective
 * outcome in @p outcomes. 
 * 
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function returns
 *   @f[
         \tr{
            \bigotimes\limits_j \ketbra{x_j}{x_j}_{q_j} \; \ketbra{\psi}{\psi} 
         }
         =
         \sum\limits_i |\psi_i|^2 \prod\limits_j \delta_{x_j, \, i_{[q_j]}}
 *   @f]
 *   where @f$\psi_i@f$ is the @f$i@f$-th amplitude of @f$\svpsi@f$, and 
 *   @f$i_{[q]}@f$ notates the @f$q@f$-th bit of @f$i@f$.
 * - When @p qureg is a density matrix @f$ \dmrho @f$, this function returns
 *   @f[
         \tr{
            \bigotimes\limits_j \ketbra{x_j}{x_j}_{q_j} \; \dmrho
         }
         =
         \sum\limits_i \re{\dmrho_{ii}} \prod\limits_j \delta_{x_j, \, i_{[q_j]}}
 *   @f]
 *   where @f$ \dmrho_{ii} @f$ is the @f$i@f$-th diagonal element of @f$\dmrho@f$. This 
 *   is real whenever @f$\dmrho@f$ is validly normalised (specifically, Hermitian).
 *
 * When @p qureg is correctly normalised, these quantities are within @f$[0, 1]@f$, and their sum
 * across all possible values of @p outcomes equals one.
 *
 * @equivalences
 * 
 * - The output of this function is equal to that found by in-turn finding the probability of each
 *   qubit being in the specified outcome, then projecting @p qureg into it (i.e. forcing that 
 *   measurement outcome). That approach is however slower and modifies @p qureg, whereas this
 *   function leaves @p qureg unchanged.
 *   ```
     qreal prob = 1;
     for (int j=0; j<numQubits; j++)
         prob *= applyForcedQubitMeasurement(qureg, qubits[j], outcomes[j]);
 *   ```
 *
 * - This function is much faster than, but mathematically equivalent to, summing the probability
 *   of every computational basis state (e.g. via calcProbOfBasisState()) which is consistent
 *   with the given qubit outcomes.
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    initRandomPureState(qureg);

    int num = 3;
    int qubits[]   = {0, 3, 4};
    int outcomes[] = {1, 1, 0};

    qreal prob = calcProbOfMultiQubitOutcome(qureg, qubits, outcomes, num);
 * ```
 *
 * @param[in] qureg     the reference state, which is unchanged.
 * @param[in] qubits    a list of target qubits to query.
 * @param[in] outcomes  a list of corresponding qubit outcomes (each `0` or `1`).
 * @param[in] numQubits the length of list @p qubits (and @p outcomes).
 * @returns The probability that the given qubits are simultaneously in the specified outcomes.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p qubits contains any duplicates.
 * - if any element of @p qubits is less than zero or beyond the number of qubits in @p qureg.
 * - if any element of @p outcomes is not `0` or `1`.
 * - if @p numQubits is less than one or exceeds the number of qubits in @p qureg.
 * @throws @segfault
 * - if either of @p qubits or @p outcomes are not lists of length @p numQubits.
* @notyetvalidated
 * @see
 * - calcProbsOfAllMultiQubitOutcomes()
 * - calcProbOfBasisState()
 * @author Tyson Jones
 */
qreal calcProbOfMultiQubitOutcome(Qureg qureg, int* qubits, int* outcomes, int numQubits);


/** Populates @p outcomeProbs with the probabilities of the specified list of @p qubits
 * being in _all_ of their possible, simultaneous outcomes (of which there are `2^`
 * @p numQubits).
 * 
 * The list @p qubits is taken to be in order of _increasing_ significance, determining 
 * the ordering of the output @p outcomeProbs.
 * For example, if @p qubits @f$ = \{ 1, 3 \} @f$, then @p outcomeProbs will be populated
 * with _four_ values; the probabilities of qubits @f$(3,1)@f$ being in the respective
 * simultaneously outcomes @f$(0,0), \, (0,1), \, (1,0) @f$ and @f$(1,1)@f$. In contrast,
 * @p qubits @f$ = \{ 3, 1 \} @f$ would see the middle two outputs swapped.
 * 
 * @formulae
 * 
 * Let @f$ n = @f$ @p numQubits, and @f$ q_i @f$ be the @f$i@f$-th element of @p qubits,
 * such that @p qubits = @f$ \{ q_0, q_1, \dots, q_{n-1} \} @f$. 
 * Let @f$ P_{\ket{q_{n-1} \dots q_1 q_0}}(\ket{i}) @f$ denote the probability that the specified
 * substate is in the computational basis substate @f$\ket{i}@f$. Explicitly, that
 * qubit @f$q_j@f$ is in the outcome given by the @f$j@f$-th bit of @f$n@f$-digit integer 
 * @f$i@f$ (simultaneously for all @f$j@f$).
 * 
 * Then, this function sets
 * @f[
      \text{outcomeProbs}[i] = P_{\ket{q_{n-1} \dots q_1 q_0}}(\ket{i})
 * @f]
 * for all @f$i \in \{0, 1, \dots 2^n-1\} @f$.
 * 
 * Explicitly, expressing substate @f$\ket{i}@f$ in terms of its individual qubits;
 * @f[
      \begin{gathered}
      \text{outcomeProbs}[0] = P_{\ket{q_{n-1} \dots q_1 q_0}}( \ket{0\dots00} ) \\
      \text{outcomeProbs}[1] = P_{\ket{q_{n-1} \dots q_1 q_0}}( \ket{0\dots01} ) \\
      \text{outcomeProbs}[2] = P_{\ket{q_{n-1} \dots q_1 q_0}}( \ket{0\dots10} ) \\
      \text{outcomeProbs}[3] = P_{\ket{q_{n-1} \dots q_1 q_0}}( \ket{0\dots11} ) \\
      \vdots \\
      \text{outcomeProbs}[2^n-1] = P_{\ket{q_{n-1} \dots q_1 q_0}}( \ket{1\dots11} )
      \end{gathered}
 * @f]
 *
 * Each probability is that which would be output by calcProbOfMultiQubitOutcome() when
 * passed @p qubits and the bits of @f$ i @f$.
 *
 * When @p qureg is correctly normalised, all probabilities are within @f$[0, 1]@f$, and
 * the sum of all elements written to @p outcomeProbs equals one.
 * 
 * @equivalences
 * 
 * - This function is significantly faster than, but otherwise equivalent to, populating
 *   each element of @p outcomeProbs in-turn with the output of calcProbOfMultiQubitOutcome().
 *   ```
     qindex numOut = (1 << numQubits);

     for (qindex i=0; i<numOut; i++) {

         // set outcomes to the bits of i
         int outcomes[numQubits];
         for (int j=0; j<numQubits; j++)
            outcomes[j] = (i >> j) & 1;

         outcomeProbs[i] = calcProbOfMultiQubitOutcome(qureg, qubits, outcomes, numQubits);
     }
 *   ``` 
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    initRandomPureState(qureg);

    int num = 3;
    int qubits[] = {0, 3, 4};
    
    qreal probs[8];
    calcProbsOfAllMultiQubitOutcomes(probs, qureg, qubits, num);
 * ```
 * @param[out] outcomeProbs the array to which the output is written.
 * @param[in]  qureg        the reference state, which is unchanged.
 * @param[in]  qubits       a list of target qubits to query.
 * @param[in]  numQubits    the length of list @p qubits.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p qubits contains any duplicates.
 * - if any element of @p qubits is less than zero or beyond the number of qubits in @p qureg.
 * - if @p numQubits is less than one or exceeds the number of qubits in @p qureg.
 * @throws @segfault
 * - if @p outcomeProbs is not a pre-allocated list of length `2^` @p numQubits.
 * - if @p qubits is not a list of length @p numQubits.
* @notyetvalidated
 * @see
 * - calcProbOfMultiQubitOutcome()
 * - calcProbOfBasisState()
 * @author Tyson Jones
 */
void calcProbsOfAllMultiQubitOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits);


/** @} */



/** 
 * @defgroup calc_properties Properties
 * @brief Functions for calculating single-state properties like normalisation and purity.
 * @{
 */


/** Calculates the probability normalisation of the given @p qureg. This is the probability
 * of the @p qureg being in _any_ outcome state, which is expected to equal `1`.
 *
 * @formulae
 * 
 * Let @f$N@f$ be the number of qubits in @p qureg.
 * 
 * - When @p qureg is a statevector @f$ \svpsi @f$ with @f$i@f$-th amplitude @f$\psi_i@f$,
 *   this function returns
 *   @f[
         \sum\limits_{i=0}^{2^N-1} |\psi_i|^2.
 *   @f]
 * - When @p qureg is a density matrix @f$ \dmrho @f$ with @f$i@f$-th diagonal element
 *   @f$ \dmrho_{ii} @f$, this function returns
 *   @f[
         \sum\limits_{i=0}^{2^N-1} \re{ \rho_{ii} }
 *   @f]
 * 
 * @constraints
 * 
 * - As above, only the real components of the diagonal elements of a density matrix are consulted;
 *   these are the only amplitudes consulted by functions which calculate probabilities in the
 *   computational basis. As such, this function gives no indication of the general validity of density
 *   matrices, such as whether they are Hermitian, whether the diagonals are real, and whether the
 *   off-diagoanl elements are valid. 
 *
 * @equivalences
 *
 * - This function is faster than, but mathematically equivalent to, summing the outputs of other
 *   functions which calculate probabilitie across all possible outcomes.
 *   ```
     // choice is arbitrary
     int qubit = 0;

     qreal totalProb = (
         calcProbOfQubitOutcome(qureg, qubit, 0) + 
         calcProbOfQubitOutcome(qureg, qubit, 1));
 *   ```
 *
 * @myexample
 * ```
    Qureg qureg = createDensityQureg(5);
    initRandomMixedState(qureg, 1<<5);

    // differs from 1 by numerical error
    qreal totalProb = calcTotalProb(qureg);
 * ```
 *
 * @param[in] qureg the reference state, which is unchanged.
 * @returns The probability normalisation of @p qureg.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
* @notyetvalidated
 * @see
 * - calcPurity()
 * - calcProbsOfAllMultiQubitOutcomes()
 * @author Tyson Jones
 */
qreal calcTotalProb(Qureg qureg);


/** Calculates the purity of @p qureg, which is a measure of its mixedness.
 *
 * @formulae
 * 
 * Let @f$N@f$ be the number of qubits in @p qureg.
 * 
 * - When @p qureg is a density matrix @f$ \dmrho @f$ (as expected), this function returns
 *   @f[
         \tr{ \dmrho^2 } = \sum\limits_{i,j} \left| \dmrho_{ij} \right|^2
 *   @f]
 *   where @f$ \dmrho_{ij} @f$ is the @f$(i,j)@f$-th element of @f$ \dmrho @f$.
 *   
 *   A purity of `1` indicates that the matrix is _pure_ and can be expressed as
 *   @f[
         \dmrho \equiv \ketbra{\phi}{\phi}
 *   @f]
 *   where @f$ \ket{\phi} @f$ is some pure state expressible as a statevector.
 *   
 *   In contrast, a purity less than `1` indicates the matrix is _mixed_ and can be
 *   understood as a convex combination of multiple (at least _two_) pure states.
 *   That is,
 *   @f[
         \dmrho \equiv \sum\limits_n p_n \ketbra{\phi}{\phi}_n,
 *   @f]
 *   where @f$p_n \in [0,1]@f$ and sum to `1` whenever @f$\dmrho@f$ is a valid and correctly
 *   normalised density matrix. Mixedness can result, for example, from @ref decoherence.
 * 
 *   The minimum purity of an @f$N@f$-qubit density matrix is @f$ 1/2^N @f$, which is
 *   admitted only by the maximally-mixed state @f$ \dmrho = \hat{\id} / 2^N @f$.
 * 
 * - When @p qureg is a statevector @f$ \svpsi @f$, this function returns
 *   @f[
         \tr{ \ketbra{\psi}{\psi} \; \ketbra{\psi}{\psi} } 
            = \left( \sum\limits_i |\psi_i|^2 \right)^2
 *   @f]
 *   where @f$\psi_i@f$ is the @f$i@f$-th amplitude of @f$\svpsi@f$. This is always `1` for
 *   any valid statevector, and is otherwise equivalent to the output of calcTotalProb(), squared.
 * 
 * @constraints
 *
 * - The output of this function is only a reliable measure of purity when @p qureg is correctly 
 *   normalised. For example, an invalid density matrix can return a purity of `1`, such as the
 *   @f$N@f$-qubit maximally-mixed state scaled by factor @f$ 2^N @f$. Note that the function 
 *   calcTotalProb() alone _cannot_ be used to validate validity since it only consults diagonal 
 *   elements, whereas the purity is informed by all elements.
 *
 * @equivalences
 *
 * - When @p qureg is a valid density matrix (specifically, Hermitian), this function is faster
 *   than, but mathematically equivalent to, calling calcInnerProduct() and passing @p qureg twice.
 *   ```
     qcomp out = calcInnerProduct(qureg, qureg);
     qreal pur = real(out); // im=0
 *   ```
 * - When @p qureg is a statevector, this function returns the output of calcTotalProb(), squared.
 *
 * @myexample
 * ```
    Qureg qureg = createDensityQureg(5);
    initRandomPureState(qureg);

    // = 1
    qreal purity1 = calcPurity(qureg);
    reportScalar("purity1", purity1);

    mixTwoQubitDepolarising(qureg, 0, 1, 0.5);

    // < 1
    qreal purity2 = calcPurity(qureg);
    reportScalar("purity2", purity2);
 * ```
 *
 * @param[in] qureg the reference state, which is unchanged.
 * @returns The purity of @p qureg.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
* @notyetvalidated
 * @see
 * - calcFidelity()
 * - calcTotalProb()
 * @author Tyson Jones
 */
qreal calcPurity(Qureg qureg);


/** @} */



/** 
 * @defgroup calc_comparisons Comparisons
 * @brief Functions for comparing multiple quantum states.
 * @{
 */


/** Calculates the fidelity between @p qureg and @p other, where at least one is a
 * statevector.
 *
 * @formulae
 * 
 * - When both @p qureg and @p other are statevectors (respectively @f$\ket{\psi}@f$ and 
 *   @f$\ket{\phi}@f$), this function returns
 *   @f[
         \left| \braket{\phi}{\psi} \right|^2.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$ and @p other is a statevector @f$\svpsi@f$,
 *   this function returns
 *   @f[
         \bra{\psi} \dmrho \ket{\psi},
 *   @f]
 *   and similarly when @p qureg is a statevector and @p other is a density matrix.
 * 
 * @constraints
 *
 * - The output of this function is always real, which validation will check after computing the
 *   fidelity as a complex scalar. Specifically, validation will assert that the result has an
 *   absolute imaginary component less than the validation epsilon, which can be adjusted with
 *   setValidationEpsilon().
 * 
 * - This function does not yet support both @p qureg and @p other being density matrices, for
 *   which the fidelity calculation is more substantial.
 * 
 * - When @p qureg and @p other are _both_ statevectors, or _both_ density matrices, then _both_ or
 *   _neither_ must be GPU-accelerated. That is, their CPU vs GPU deployments must agree. They are
 *   permitted to differ in distribution however. Such considerations are only relevant when
 *   creating the registers using createCustomQureg(), since the automatic deployments of createQureg()
 *   and createDensityQureg() will always agree.
 * 
 * - When @p qureg and @p other dimensionally _differ_ (i.e. one is a statevector while the other is a
 *   density matrix), the statevector must not be distributed _unless_ the density matrix is distributed.
 *   The CPU vs GPU deployments however are permitted to disagree. These requirements are again
 *   consistent with the automatic deployments of the createQureg() and createDensityQureg() functions.
 * 
 * @equivalences
 * 
 * - When both @p qureg and @p other are statevectors, this function is equivalent to calling
 *   calcInnerProduct() and squaring the absolute value of the result.
 *   ```
     qcomp prod = calcInnerProduct(qureg, other);
     qreal fid = pow(abs(prod), 2);
 *   ```
 * - When one of @p qureg or @p other is a statevector in the computational basis state @f$\ket{i}@f$
 *   (e.g. as can be produced via initClassicalState()), this function is slower but equivalent to 
 *   finding directly the probability of the basis state.
 *   ```
     // initClassicalState(other, index);

     qreal fid = calcProbOfBasisState(qureg, index);
 *   ```
 *
 * @myexample
 * ```
   // rho = |psi><psi|
   Qureg psi = createQureg(5);
   Qureg rho = createDensityQureg(5);
   initRandomPureState(psi);
   initPureState(rho, psi);

   qreal fid0 = calcFidelity(rho, psi); // = 1

   mixDepolarising(rho, 0, 0.5);
   qreal fid1 = calcFidelity(rho, psi); // < 1
 * ```
 *
 * @param[in] qureg a state
 * @param[in] other another state containing an equal number of qubits.
 * @returns The fidelity between @p qureg and @p other.
 * @throws @validationerror
 * - if @p qureg or @p other is uninitialised.
 * - if @p qureg and @p other contain a different number of qubits.
 * - if @p qureg and @p other are incompatible deployed.
 * - if both @p qureg and @p other are density matrices (as is not yet supported).
 * - if @p qureg or @p other is unnormalised such that the calculated fidelity is non-real.
 * @notyetvalidated
 * @see
 * - calcInnerProduct()
 * - calcDistance()
 * @author Tyson Jones
 */
qreal calcFidelity(Qureg qureg, Qureg other);


/** Calculates one of three distance measures between @p qureg and @p other, depending
 * upon whether one or both are density matrices. These are the Hilbert-Schmidt distance,
 * Bures distance and purified distance.
 *
 * @formulae
 * 
 * - When both @p qureg and @p other are statevectors (respectively @f$\ket{\psi}@f$ and 
 *   @f$\ket{\phi}@f$), this function returns the **Bures distance** defined as
 *   @f[
         d_B\left(\ket{\psi},\ket{\phi}\right) = \sqrt{2 - 2 \left| \braket{\phi}{\psi} \right|}
 *   @f]
 *   where @f$\left| \braket{\phi}{\psi} \right|@f$ is the square-root of the fidelity
 *   between @f$\ket{\psi}@f$ and @f$\ket{\phi}@f$ as would be computed by calcFidelity().
 *
 * - When both @p qureg and @p other are density matrices (respectively @f$\mathbf{\rho}@f$
 *   and @f$\mathbf{\sigma}@f$), this function returns the **Hilbert-Schmidt distance** defined as
 *   @f[
         d_{HS}\left(\mathbf{\rho}, \mathbf{\sigma}\right) 
            = 
            \sqrt{ \tr{
               \left| \mathbf{\rho} - \mathbf{\sigma} \right|^2
            } }
            =
            \sqrt{
               \sum\limits_{ij} \left| \rho_{ij} - \sigma_{ij} \right|^2
            }.
 *   @f]
 *
 * - When one of @p qureg or @p other is a statevector @f$\svpsi@f$, and the other is a density
 *   matrix @f$\dmrho@f$, this function returns the **purified distance** defined as
 *   @f[
         d_p\left(\svpsi,\dmrho\right) = \sqrt{ 1 - \brapsi \dmrho \svpsi }
 *   @f]
 *   where @f$\brapsi \dmrho \svpsi@f$ is the fidelity as returned by calcFidelity().
 * 
 * @constraints
 * 
 * - The output of this function is always real, which is always mathematically satisfied by the
 *   Hilbert-Schmidt distance, but may be violated by the Bures and purified distances when the
 *   input Qureg are not normalised, or otherwise due to numerical imprecision. Postcondition
 *   validation of the Bures distance will check that
 *   @f[
         \left| \braket{\phi}{\psi} \right| \le 1 + \valeps
 *   @f]
 *   while the purified distance validation will check that
 *   @f[
         \left| \, \im{ \brapsi \dmrho \svpsi } \, \right| \le \valeps, \\
         \re{ \brapsi \dmrho \svpsi } \le 1 + \valeps,
 *   @f]
 *   where @f$\valeps@f$ is the validation epsilon, adjustable via setValidationEpsilon().
 * 
 * - Even when the above postcondition validation is disabled, the Bures and purified distance
 *   calculations will respectively replace @f$\left| \braket{\phi}{\psi} \right|@f$ and 
 *   @f$\re{ \brapsi \dmrho \svpsi }@f$ which exceed @f$1@f$ with value @f$1@f$, and the imaginary
 *   component of @f$\brapsi \dmrho \svpsi@f$ is discarded.
 * 
 * - When @p qureg and @p other are _both_ statevectors, or _both_ density matrices, then _both_ or
 *   _neither_ must be GPU-accelerated. That is, their CPU vs GPU deployments must agree. They are
 *   permitted to differ in distribution however. Such considerations are only relevant when
 *   creating the registers using createCustomQureg(), since the automatic deployments of createQureg()
 *   and createDensityQureg() will always agree.
 * 
 * - When @p qureg and @p other dimensionally _differ_ (i.e. one is a statevector while the other is a
 *   density matrix), the statevector must not be distributed _unless_ the density matrix is distributed.
 *   The CPU vs GPU deployments however are permitted to disagree. These requirements are again
 *   consistent with the automatic deployments of the createQureg() and createDensityQureg() functions.
 * 
 * @equivalences
 * 
 * - When both @p qureg and @p other are statevectors, this function wraps calcInnerProduct().
 *   ```
     qcomp prod = calcInnerProduct(qureg, other); // <qureg|other>
     qreal mag = abs(prod);
     mag = (mag > 1)? 1 : mag;
     qreal dist = std::sqrt(2 - 2 * mag);
 *   ```
 *
 * - When @p qureg is a density matrix and @p other is a statevector, this function wraps calcInnerProduct()
 *   as a complex-valued proxy for calcFidelity().
 *   ```
     qcomp prod = calcInnerProduct(other, qureg); // <other|qureg|other>
     qreal re = real(prod);
     re = (re > 1)? 1 : re;
     qreal dist = sqrt(1 - re);
 *   ```
 *
 * @myexample
 * ```
   Qureg rho1 = createDensityQureg(5);
   Qureg rho2 = createDensityQureg(5);

   initRandomMixedState(rho1, 10);
   setQuregToClone(rho2, rho1);
   qreal distA = calcDistance(rho1, rho2); // = 0

   initRandomMixedState(rho2, 10);
   qreal distB = calcDistance(rho1, rho2); // > 0
 * ```
 *
 * @param[in] qureg a state
 * @param[in] other another state containing an equal number of qubits
 * @returns The distance between @p qureg and @p other, according to the above measures.
 * @throws @validationerror
 * - if @p qureg or @p other is uninitialised.
 * - if @p qureg and @p other contain a different number of qubits.
 * - if @p qureg and @p other are incompatible deployed.
 * - if @p qureg or @p other is unnormalised such that the Bures or purified distances would be non-real.
 * @notyetvalidated
 * @see
 * - calcInnerProduct()
 * - calcFidelity()
 * @author Tyson Jones
 */
qreal calcDistance(Qureg qureg, Qureg other);


/** @} */



/** 
 * @defgroup calc_partialtrace Partial trace
 * @brief Functions for calculating reduced density matrices, creating a new output Qureg.
 * @{
 */


/** Creates and populates a new Qureg which is a reduced density matrix resulting from tracing out 
 * the specified qubits of @p qureg. This should be later freed by the user like all Qureg.
 * 
 * Note that the deployments of the output Qureg (i.e. whether multithreaded, GPU-accelerated and
 * distributed) will match those of @p qureg. It is ergo intended that this function is used to
 * trace out few qubits, and may show worsening performance when tracing many qubits.
 * 
 * The ordering of @p traceOutQubits has no effect, and the ordering of the remaining qubits in
 * the output Qureg match their original relative ordering in @p qureg.
 * 
 * @formulae
 * 
 * Let @f$\dmrho_{\text{in}} = @f$ @p qureg and let @f$\vec{t} = @f$ @p traceOutQubits which is a list of
 * length @f$n = @f$ @p numTraceQubits.
 * 
 * This function returns a new Qureg @f$\dmrho_{\text{out}}@f$ which satisfies
 * @f[
        \dmrho_{\text{out}} = \text{Tr}_{\vec{t}} \left( \dmrho_{\text{in}} \right)
        =
        \sum\limits_i^{2^n} 
        (\hat{\id} \otimes \bra{i}_{\vec{t}} ) \,
         \dmrho_{\text{in}} \,
        (\hat{\id} \otimes \ket{i}_{\vec{t}} )
 * @f]
 * where @f$\ket{i}_{\vec{t}}@f$ notates the @f$i@f$-th basis state (in any orthonormal basis) of the
 * targeted qubits, and @f$(\hat{\id} \otimes \ket{i}_{\vec{t}})@f$ notates interleaved identity operators
 * upon the non-targeted qubits.
 * 
 * Given an @f$N@f$-qubit Qureg @f$\dmrho_{\text{in}}@f$, the output @f$\dmrho_{\text{out}}@f$ contains
 * @f$N-n@f$ qubits.
 * 
 * @constraints
 * 
 * - The given @p qureg must be a density matrix. It is however straightforward to prepare a density matrix
 *   from a statevector.
 *   ```
     // let qureg be the intended initial statevector

     Qureg temp = createDensityQureg(qureg.numQubits);
     initPureState(temp, qureg);

     Qureg reduced = calcPartialTrace(temp, traceOutQubits, numTraceQubits);
     destroyQureg(temp);
 *   ```
 * 
 * - When @p qureg is distributed, the returned Qureg will also be distributed, which imposes a minimum on
 *   the number of qubits contained within; @f$\log_2(W)@f$ where @f$W@f$ is the number of distributed nodes
 *   (or "world size"). This imposes a maximum upon @p traceOutQubits of
 *   ```
 *   numTraceQubits <= qureg.numQubits - qureg.logNumNodes
 *   ```
 * 
 * @equivalences
 * 
 * - The function calcReducedDensityMatrix() is entirely equivalent, but conveniently permits specifying
 *   a list of which qubits to _retain_ during partial tracing. 
 * 
 * - The functions setQuregToPartialTrace() and setQuregToReducedDensityMatrix() are also equivalent but
 *   permit overwriting an existing Qureg.
 *  
 * @myexample
 * 
 * ```
   Qureg state = createDensityQureg(5);
   initRandomMixedState(state, 10);
   reportQureg(state);

   int qubits[] = {0,2,4};
   Qureg reduced = calcPartialTrace(state, qubits, 3);
   reportQureg(reduced);

   // state's qubits {1,3} have become reduced's qubits {0,1}
 * ```
 * 
 * @param[in] qureg          a density matrix which is not modified.
 * @param[in] traceOutQubits a list of qubits to trace out and ergo from the output Qureg.
 * @param[in] numTraceQubits the length of @p traceOutQubits.
 * @returns A new, smaller Qureg initialised to the reduced density matrix of @p qureg.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p numTraceQubits is less than one.
 * - if @p numTraceQubits is equal or greater than the number of qubits in @p qureg.
 * - if @p qureg is distributed and @p numTraceQubits exceeds `qureg.numQubits - qureg.logNumNodes`.
 * - if the system contains insufficient RAM (or VRAM) to store the new Qureg in any deployment.
 * - if any memory allocation of the output Qureg unexpectedly fails.
 * @throws seg-fault
 * - if @p traceOutQubits is not a list of length @p numTraceQubits.
 * @notyetvalidated
 * @see
 * - calcReducedDensityMatrix()
 * - setQuregToPartialTrace()
 * - setQuregToReducedDensityMatrix()
 * @author Tyson Jones
 */
Qureg calcPartialTrace(Qureg qureg, int* traceOutQubits, int numTraceQubits);


/** Creates and populates a new Qureg which is a reduced density matrix of @p qureg,
 * retaining only the specified qubits and tracing out all others.
 * 
 * Note that the deployments of the output Qureg (i.e. whether multithreaded, GPU-accelerated and
 * distributed) will match those of @p qureg. It is ergo intended that this function is used to
 * preserve most qubits of @p qureg, and may show worsening performance when retaining only few.
 * 
 * > [!CAUTION]
 * > The ordering of @p retainQubits has no effect on the output state. The ordering of the
 * > retained qubits will match their original, relative ordering in @p qureg.
 *
 * @formulae
 * 
 * This function is entirely equivalent to calcPartialTrace() except that here the _retained_ qubits
 * are specified, whereas calcPartialTrace() accepts those to be traced out.
 * 
 * Let @f$\dmrho_{\text{in}} = @f$ @p qureg, @f$\vec{r} = @f$ @p retainQubits, and let @f$\vec{q}@f$
 * be a list containing _all_ qubits of @p qureg. This function partially traces out all qubits in
 * list @f$\vec{t} = \vec{q} \setminus \vec{r}@f$, and returns a new Qureg @f$\dmrho_{\text{out}}@f$ 
 * which satisfies
 * @f[
        \dmrho_{\text{out}} = \text{Tr}_{\vec{t}} \left( \dmrho_{\text{in}} \right)
        =
        \sum\limits_i^{2^n} 
        (\hat{\id} \otimes \bra{i}_{\vec{t}} ) \,
         \dmrho_{\text{in}} \,
        (\hat{\id} \otimes \ket{i}_{\vec{t}} )
 * @f]
 * where @f$\ket{i}_{\vec{t}}@f$ notates the @f$i@f$-th basis state (in any orthonormal basis) of the
 * qubits in @f$\vec{t}@f$, and @f$(\hat{\id} \otimes \ket{i}_{\vec{t}})@f$ notates interleaved identity
 * operators upon the qubits in @f$\vec{r}@f$.
 * 
 * @constraints
 * 
 * - The given @p qureg must be a density matrix. It is however straightforward to prepare a density matrix
 *   from a statevector.
 *   ```
     // let qureg be the intended initial statevector

     Qureg temp = createDensityQureg(qureg.numQubits);
     initPureState(temp, qureg);

     Qureg reduced = calcReducedDensityMatrix(temp, retainQubits, numRetainQubits);
     destroyQureg(temp);
 *   ```
 * 
 * - When @p qureg is distributed, the returned Qureg will also be distributed, which imposes a minimum on
 *   the number of qubits contained within; @f$\log_2(W)@f$ where @f$W@f$ is the number of distributed nodes
 *   (or "world size"). This imposes bounds upon @p numRetainQubits of
 *   ```
 *   qureg.logNumNodes <= numRetainQubits <= qureg.numQubits - 1
 *   ```
 *
 * @equivalences
 * 
 * - The function calcPartialTrace() is entirely equivalent, but permits directly specifying the qubits to
 *   be traced out.
 * 
 * - The functions setQuregToPartialTrace() and setQuregToReducedDensityMatrix() are also equivalent but
 *   permit overwriting an existing Qureg.
 *  
 * @myexample
 * 
 * ```
   Qureg state = createDensityQureg(5);
   initRandomMixedState(state, 10);
   reportQureg(state);

   int qubits[] = {1,3};
   Qureg reduced = calcReducedDensityMatrix(state, qubits, 2);
   reportQureg(reduced);

   // state's qubits {1,3} have become reduced's qubits {0,1}
 * ```
 * 
 * @param[in] qureg            a density matrix.
 * @param[in] retainQubits    a list of qubits to retain in the reduced density matrix (at shifted, contiguous indices).
 * @param[in] numRetainQubits the length of @p retainQubits.
 * @returns A new Qureg containing @p numRetainQubits qubits, initialised to the reduced density matrix of @p qureg.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p numRetainQubits is less than one.
 * - if @p numRetainQubits is equal or greater than the number of qubits in @p qureg.
 * - if @p qureg is distributed and @p numRetainQubits is less than `qureg.logNumNodes`.
 * - if the system contains insufficient RAM (or VRAM) to store the new Qureg in any deployment.
 * - if any memory allocation of the output Qureg unexpectedly fails.
 * @throws seg-fault
 * - if @p retainQubits is not a list of length @p numRetainQubits.
 * @notyetvalidated
 * @see
 * - calcPartialTrace()
 * - setQuregToPartialTrace()
 * - setQuregToReducedDensityMatrix()
 * @author Tyson Jones
 */
Qureg calcReducedDensityMatrix(Qureg qureg, int* retainQubits, int numRetainQubits);


/** @} */


// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because they pass or
 * return qcomp primitives by-value (rather than by pointer).
 * This is prohibited because the C and C++ ABI does not agree
 * on a complex type, though C's _Complex has the same memory
 * layout as C++'s std::complex<>. To work around this, the 
 * below functions have a C-compatible wrapper defined in
 * wrappers.h which passes/receives the primitives by pointer;
 * a qcomp ptr can be safely passed from the C++ source binary
 * the user's C binary. We manually add these functions to their
 * respective Doxygen doc groups defined above
 */


/** @ingroup calc_comparisons
 * 
 * Calculates the inner product of state @p qureg with @p other. 
 *
 * @formulae
 * 
 * - When both @p qureg and @p other are statevectors (respectively @f$\ket{\psi}@f$ and 
 *   @f$\ket{\phi}@f$), this function returns
 *   @f[
         \braket{\psi}{\phi} = \sum\limits_i \psi_i^* \phi_i
 *   @f]
 *   where @f$\psi_i@f$ and @f$\phi_i@f$ are the @f$i@f$-th amplitudes of @f$\ket{\psi}@f$ 
 *   (@p qureg) and  @f$\ket{\phi}@f$ (@p other) respectively, and @f$\alpha^*@f$ notates
 *   the complex conjugate of scalar @f$\alpha@f$.
 * 
 * - When both @p qureg and @p other are density matrices (respectively @f$\mathbf{\rho}@f$
 *   and @f$\mathbf{\sigma}@f$), this function returns
 *   @f[
         \tr{ \rho^\dagger \sigma } = \sum\limits_{ij} {\rho_{ij}}^* \, \sigma_{ij}.
 *   @f]
 * 
 * - When @p qureg is a density matrix @f$\dmrho@f$ and @p other is a statevector @f$\ket{\phi}@f$,
 *   this function returns
 *   @f[
         \bra{\phi} \dmrho^\dagger \ket{\phi}.
 *   @f]
 *
 * - When @p qureg is a statevector @f$\svpsi@f$ and @p other is a density matrix @f$\mathbf{\sigma}@f$,
 *   this function returns
 *   @f[
         \brapsi \mathbf{\sigma} \svpsi.
 *   @f]
 *
 * @constraints
 * 
 * - When @p qureg and @p other are _both_ statevectors, or _both_ density matrices, then _both_ or
 *   _neither_ must be GPU-accelerated. That is, their CPU vs GPU deployments must agree. They are
 *   permitted to differ in distribution however. Such considerations are only relevant when
 *   creating the registers using createCustomQureg(), since the automatic deployments of createQureg()
 *   and createDensityQureg() will always agree.
 * 
 * - When @p qureg and @p other dimensionally _differ_ (i.e. one is a statevector while the other is a
 *   density matrix), the statevector must not be distributed _unless_ the density matrix is distributed.
 *   The CPU vs GPU deployments however are permitted to disagree. These requirements are again
 *   consistent with the automatic deployments of the createQureg() and createDensityQureg() functions.
 *
 * @myexample
 * ```
   Qureg rho1 = createDensityQureg(5);
   Qureg rho2 = createDensityQureg(5);

   // rho1 = rho2 = |psi><psi|
   initRandomPureState(rho1);
   setQuregToClone(rho2, rho1);
   qcomp prodA = calcInnerProduct(rho1, rho2); // = 1

   // rho1 = rho2 = sum_i prob_i |psi_i><psi_i|
   initRandomMixedState(rho1, 10);
   setQuregToClone(rho2, rho1);
   qcomp prodB = calcInnerProduct(rho1, rho2); // < 1, real

   // rho1 != rho2
   initRandomMixedState(rho2, 10);
   qcomp prodC = calcInnerProduct(rho1, rho2); // abs < 1, complex
 * ```
 *
 * @param[in] qureg a state
 * @param[in] other another state with an equal number of qubits
 * @returns The inner product of @p qureg with @p other.
 * @throws @validationerror
 * - if @p qureg or @p other is uninitialised.
 * - if @p qureg and @p other contain a different number of qubits.
 * - if @p qureg and @p other are incompatibly deployed.
 * @notyetvalidated
 * @see
 * - calcDistance()
 * - calcFidelity()
 * @author Tyson Jones
 */
qcomp calcInnerProduct(Qureg qureg, Qureg other);


/** @ingroup calc_expec
 * 
 * Calculates the expectation value of the given permittedly non-Hermitian operator @p sum 
 * - a weighted sum of Pauli strings with complex weights - under the given state @p qureg, 
 * which is not modified.
 * 
 * @formulae
 * 
 * This function is mathematically equivalent to calcExpecPauliStrSum(), _except_ that here a
 * complex scalar is returned. This permits obtaining the full scalar when @p sum contains non-real
 * weights, and/or when @p qureg is unnormalised.
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    PauliStrSum sum = createInlinePauliStrSum(R"(
        0.123 + 3.5i  ZIZIZI
        1.234 - 1E-5i XYZXZ
        -1E-2         IIIII
    )");

    // prints "expec: 0.113+3.5i"
    qcomp expec = calcExpecNonHermitianPauliStrSum(qureg, sum);
    reportScalar("expec", expec);  
 * ```
 *
 * @param[in] qureg the permittedly unnormalised reference state.
 * @param[in] sum   the permittedly non-Hermitian operator.
 * @returns The permittedly complex expectation value.
 * @throws @validationerror
 * - if @p qureg or @p sum are uninitialised.
 * - if any PauliStr in @p sum targets a higher-index qubit than exists in @p qureg.
* @notyetvalidated
 * @see
 * - calcExpecPauliStrSum()
 * @author Tyson Jones
 */
qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum); 


/** @ingroup calc_expec
 * 
 * Calculates the expectation value of the given permittedly non-Hermitian operator @p matr,
 * under the given state @p qureg, without modifying it. 
 * 
 * @formulae
 * 
 * This function is mathematically equivalent to calcExpecFullStateDiagMatr(), _except_ that here a
 * complex scalar is returned. This permits obtaining the full scalar when @p sum contains non-real
 * elements, and/or when @p qureg is unnormalised.
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    initPlusState(qureg);

    FullStateDiagMatr matr = createFullStateDiagMatr(qureg.numQubits);

    // profanely inefficient per-element initialisation
    for (int n=0; n<matr.numElems; n++) {
        qcomp elem = getQcomp(n, n+1);
        setFullStateDiagMatr(matr, n, &elem, 1);
    }

    // prints "expec: 15.5+16.5i"
    qcomp expec = calcExpecNonHermitianFullStateDiagMatr(qureg, matr);
    reportScalar("expec", expec);
 * ```
 *
 * @param[in] qureg the permittedly unnormalised reference state.
 * @param[in] matr  the permittedly non-Hermitian operator.
 * @returns The permittedly complex expectation value.
 * @throws @validationerror
 * - if @p qureg or @p matr are uninitialised.
 * - if @p matr does not match the dimension of @p qureg
 * - if @p matr is distributed but @p qureg is not
* @notyetvalidated
 * @see
 * - calcExpecFullStateDiagMatr()
 * - calcExpecFullStateDiagMatrPower()
 * - calcExpecNonHermitianFullStateDiagMatrPower()
 * @author Tyson Jones
 */
qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);


/** @ingroup calc_expec
 * 
 * Calculates the expectation value of the given permittedly non-Hermitian operator @p matrix,
 * raised to the arbitrary complex @p exponent, under the given state @p qureg, which is not modified.
 * 
 * @formulae
 * 
 * This function is mathematically equivalent to calcExpecFullStateDiagMatrPower(), _except_ that 
 * here a complex scalar is returned, in addition to @p exponent being permittedly complex.
 * This permits obtaining the full scalar when @p qureg is unnormalised or @p matrix (after being
 * raised to @p exponent) is non-Hermitian.
 *
 * @myexample
 * ```
    Qureg qureg = createQureg(5);
    initPlusState(qureg);

    FullStateDiagMatr matrix = createFullStateDiagMatr(qureg.numQubits);

    // profanely inefficient per-element initialisation
    for (int n=0; n<matrix.numElems; n++) {
        qcomp elem = getQcomp(n, n+1);
        setFullStateDiagMatr(matrix, n, &elem, 1);
    }

    qcomp exponent = 3+4_i;

    // prints "expec: -257.26-613.8i"
    qcomp expec = calcExpecNonHermitianFullStateDiagMatrPower(qureg, matrix, exponent);
    reportScalar("expec", expec);
 * ```
 *
 * @param[in] qureg     the permittedly unnormalised reference state.
 * @param[in] matrix    the permittedly non-Hermitian operator.
 * @param[in] exponent  the permittedly complex exponent.
 * @returns The permittedly complex expectation value.
 * @throws @validationerror
 * - if @p qureg or @p matrix are uninitialised.
 * - if @p matrix does not match the dimension of @p qureg
 * - if @p matrix is distributed but @p qureg is not
* @notyetvalidated
 * @see
 * - calcExpecFullStateDiagMatr()
 * @author Tyson Jones
 */
qcomp calcExpecNonHermitianFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

#include <vector>


/// @ingroup calc_prob
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see calcProbOfMultiQubitOutcome()
qreal calcProbOfMultiQubitOutcome(Qureg qureg, std::vector<int> qubits, std::vector<int> outcomes);


/// @ingroup calc_prob
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
/// @cppvectoroverload
/// @see calcProbsOfAllMultiQubitOutcomes()
std::vector<qreal> calcProbsOfAllMultiQubitOutcomes(Qureg qureg, std::vector<int> qubits);


/// @ingroup calc_partialtrace
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see calcPartialTrace()
Qureg calcPartialTrace(Qureg qureg, std::vector<int> traceOutQubits);


/// @ingroup calc_partialtrace
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see calcReducedDensityMatrix()
Qureg calcReducedDensityMatrix(Qureg qureg, std::vector<int> retainQubits);


#endif // __cplusplus


#endif // CALCULATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
