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


// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


/// @notvalidated
qreal calcExpecPauliStr(Qureg qureg, PauliStr str);

/// @notvalidated
qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum);

/// @notvalidated
qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);

/// @notvalidated
qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matr, qreal exponent);


/// @notvalidated
qreal calcTotalProb(Qureg qureg);

/// @notvalidated
qreal calcProbOfBasisState(Qureg qureg, qindex index);

/// @notvalidated
qreal calcProbOfQubitOutcome(Qureg qureg, int qubit, int outcome);

/// @notvalidated
qreal calcProbOfMultiQubitOutcome(Qureg qureg, int* qubits, int* outcomes, int numQubits);

/// @notvalidated
void  calcProbsOfAllMultiQubitOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits);


/// @notvalidated
qreal calcPurity(Qureg qureg);

/// @notvalidated
qreal calcFidelity(Qureg qureg, Qureg other);

/// @notvalidated
qreal calcDistance(Qureg qureg1, Qureg qureg2);


/// @notvalidated
Qureg calcPartialTrace(Qureg qureg, int* traceOutQubits, int numTraceQubits);

/// @notvalidated
Qureg calcReducedDensityMatrix(Qureg qureg, int* retainQubits, int numRetainQubits);


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
 * the user's C binary. 
 */

/// @notvalidated
qcomp calcInnerProduct(Qureg qureg1, Qureg qureg2);

/// @notvalidated
qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum); 

/// @notvalidated
qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);

/// @notvalidated
qcomp calcExpecNonHermitianFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);


#endif // CALCULATIONS_H

/** @} (end doxygen defgroup) */
