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
 * Let @f$ \pstr = @f$ @p str.
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
 * - The returned value is always real, even when @p qureg is an unnormalised density matrix, in
 *   which case the imaginary component of the above expression is neglected.
 *   The full complex value can be obtained using calcExpecNonHermitianPauliStrSum().
 * 
 * @equivalence
 * - When @p str is general, this function is equivalent to calling calcExpecPauliStrSum() with a 
 *   PauliStrSum composed of only a single PauliStr term and a unity coefficient.
 * - When @p str @f$ = \id^\otimes @f$, the output is equivalent to that of calcTotalProb().
 * 
 * @myexample
 * ```
    Qureg qureg = createQureg(4);
    PauliStr str = getPauliStr("XYZ");

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
 * @throws invalidQuESTInputError()
 * - if @p qureg is uninitialised.
 * - if @p str contains a (non-identity) Pauli upon a higher-index qubit than exists in @p qureg.
 * - if the output (with unreturned imaginary component) is not approximately real.
 * @notvalidated
 * @author Tyson Jones
 */
qreal calcExpecPauliStr(Qureg qureg, PauliStr str);


/** Calculates the expectation value of the given Hermitian observable @p sum - a weighted sum of 
 * Pauli strings - under the given state @p qureg, without modifying it. 
 * 
 * @formulae
 * Let @f$ \hat{H} = @f$ @p sum.
 * - When @p qureg is a statevector @f$\svpsi@f$, this function returns
 *   @f[ 
    \brapsi \hat{H} \svpsi \in \mathbb{R}.
 *   @f]
 * - When @p qureg is a density matrix @f$\dmrho@f$, this function returns the real component of
 *   @f[ 
     \tr{ \hat{H} \dmrho }
 *   @f]
 *   which is the exact expectation value when @f$\dmrho@f$ is physical (specifically Hermitian).
 * 
 * @constraints
 * - Hermiticity of @p sum requires that every coefficient within is real. 
 *   Validation will check @p sum is _approximately_ Hermitian, i.e. that
 *   @f[ 
     |\im{c}| \le \valeps
 *   @f]
 *   for all @f$c \in @f$ `sum.coeffs`. Adjust @f$\valeps@f$ using setValidationEpsilon().
 * - The returned value is always real, and the imaginary component is neglected even when 
 *   Hermiticity validation is relaxed and/or @p qureg is an unnormalised density matrix. 
 *   The full complex value can be obtained using calcExpecNonHermitianPauliStrSum().
 * 
 * @equivalence
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
 * @param[in] qureg the reference state.
 * @param[in] sum   the observable operator.
 * @returns The real component of the expectation value.
 * @throws invalidQuESTInputError()
 * - if @p qureg or @p sum are uninitialised.
 * - if any PauliStr in @p sum targets a higher-index qubit than exists in @p qureg.
 * - if @p sum is not approximately Hermitian.
 * - if the output (with unreturned imaginary component) is not approximately real.
* @notvalidated
 * @see
 * - calcExpecNonHermitianPauliStrSum()
 * - calcExpecFullStateDiagMatr()
 * @author Tyson Jones
 */
qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum);


/// @notdoced
/// @notvalidated
qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);


/// @notdoced
/// @notvalidated
qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matr, qreal exponent);


/** @} */



/** 
 * @defgroup calc_prob Probabilities
 * @brief Functions for non-destructively calculating the probabilities of measurement outcomes.
 * @{
 */


/// @notdoced
/// @notvalidated
qreal calcProbOfBasisState(Qureg qureg, qindex index);


/// @notdoced
/// @notvalidated
qreal calcProbOfQubitOutcome(Qureg qureg, int qubit, int outcome);


/// @notdoced
/// @notvalidated
qreal calcProbOfMultiQubitOutcome(Qureg qureg, int* qubits, int* outcomes, int numQubits);


/// @notdoced
/// @notvalidated
void calcProbsOfAllMultiQubitOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits);


/** @} */



/** 
 * @defgroup calc_properties Properties
 * @brief Functions for calculating single-state properties like normalisation and purity.
 * @{
 */


/// @notdoced
/// @notvalidated
qreal calcTotalProb(Qureg qureg);


/// @notdoced
/// @notvalidated
qreal calcPurity(Qureg qureg);


/** @} */



/** 
 * @defgroup calc_comparisons Comparisons
 * @brief Functions for comparing multiple quantum states.
 * @{
 */


/// @notdoced
/// @notvalidated
qreal calcFidelity(Qureg qureg, Qureg other);

/// @notdoced
/// @notvalidated
qreal calcDistance(Qureg qureg1, Qureg qureg2);


/** @} */



/** 
 * @defgroup calc_partialtrace Partial trace
 * @brief Functions for calculating reduced density matrices, creating a new output Qureg.
 * @{
 */


/// @notdoced
/// @notvalidated
Qureg calcPartialTrace(Qureg qureg, int* traceOutQubits, int numTraceQubits);


/// @notdoced
/// @notvalidated
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


/// @ingroup calc_comparisons
/// @notdoced
/// @notvalidated
qcomp calcInnerProduct(Qureg qureg1, Qureg qureg2);


/// @ingroup calc_expec
/// @notdoced
/// @notvalidated
qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum); 


/// @ingroup calc_expec
/// @notdoced
/// @notvalidated
qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);


/// @ingroup calc_expec
/// @notdoced
/// @notvalidated
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
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
qreal calcProbOfMultiQubitOutcome(Qureg qureg, std::vector<int> qubits, std::vector<int> outcomes);


/// @ingroup calc_prob
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
std::vector<qreal> calcProbsOfAllMultiQubitOutcomes(Qureg qureg, std::vector<int> qubits);


/// @ingroup calc_partialtrace
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
Qureg calcPartialTrace(Qureg qureg, std::vector<int> traceOutQubits);


/// @ingroup calc_partialtrace
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
Qureg calcReducedDensityMatrix(Qureg qureg, std::vector<int> retainQubits);


#endif // __cplusplus


#endif // CALCULATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
