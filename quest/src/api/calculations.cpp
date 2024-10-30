/** @file
 * API definitions for calculating properties of quantum states,
 * such as probabilities and expectation values.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/calculations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/bitwise.hpp"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions



/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because of limited
 * interoperability of the qcomp type. See calculations.h 
 * for more info. We here define a C++-only signature (with
 * name-mangling), and a C-friendly wrapper which passes by
 * pointer; the C-friendly interface in wrappers.h which itself
 * wrap this.
 */


qcomp calcInnerProduct(Qureg qureg1, Qureg qureg2) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_calcInnerProduct(Qureg bra, Qureg ket, qcomp* out) {

    *out = calcInnerProduct(bra, ket);
}


qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_calcExpecNonHermitianPauliStrSum(qcomp* out, Qureg qureg, PauliStrSum sum) {

    *out = calcExpecNonHermitianPauliStrSum(qureg, sum);
}


qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_calcExpecNonHermitianFullStateDiagMatr(qcomp* out, Qureg qureg, FullStateDiagMatr matr) {

    *out = calcExpecNonHermitianFullStateDiagMatr(qureg, matr);
}



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {


/*
 * EXPECTED VALUES
 */

qreal calcExpecPauliStr(Qureg qureg, PauliStr str) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum) {

    // TODO
    valdidate_pauliStrSumIsHermitian(sum, __func__);

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr) {

    // TODO
    validate_matrixIsHermitian(matr, __func__);

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


/*
 * PROBABILITIES
 */


qreal calcProbOfBasisState(Qureg qureg, qindex index) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


qreal calcProbOfQubitOutcome(Qureg qureg, int qubit, int outcome) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


qreal calcProbOfMultiQubitOutcome(Qureg qureg, int* qubits, int* outcomes, int numQubits) {

    // TODO
    error_functionNotImplemented(__func__);
}


void calcProbsOfAllMultiQubitOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {
}


qreal calcTotalProb(Qureg qureg) {

    // TODO
    error_functionNotImplemented(__func__);    
    return -1;
}



/*
 * MEASURES
 */


qreal calcPurity(Qureg qureg) {

    // TODO: how can I have an equivalent that works when density-matrix qureg is unnormalised?
    // Can it call calcInnerProduct(qureg,qureg)??

    // TODO
    error_functionNotImplemented(__func__);    
    return -1;
}


qreal calcFidelity(Qureg qureg, Qureg other) {

    // TODO: how can I have an equivalent that works when density-matrix qureg is unnormalised?
    // Just direct users to use calcInnerProduct(qureg, other)??

    // TODO
    error_functionNotImplemented(__func__);    
    return -1;
}


qreal calcHilbertSchmidtDistance(Qureg qureg1, Qureg qureg2) {

    // TODO: we could eval this for pure sttes too - waht is it?


    // TODO
    error_functionNotImplemented(__func__);    
    return -1;
}



/*
 * PARTIAL TRACE
 */

Qureg calcPartialTrace(Qureg qureg, int* qubits, int numQubits) {

    // BIG TODO
    error_functionNotImplemented(__func__);    
    return createQureg(0);
}


} // end de-name mangler
