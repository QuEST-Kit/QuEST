/** @file
 * API definitions for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/channels.h"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions

// enable invocation by both C and C++ binaries
extern "C" {



void mixDephasing(Qureg qureg, int qubit, qreal prob) {

    // TODO
    error_functionNotImplemented(__func__);
}


void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    // TODO
    error_functionNotImplemented(__func__);
}


void mixDepolarising(Qureg qureg, int qubit, qreal prob) {

    // TODO
    error_functionNotImplemented(__func__);
}


void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    // TODO
    error_functionNotImplemented(__func__);
}


void mixDamping(Qureg qureg, int qubit, qreal prob) {

    // TODO
    error_functionNotImplemented(__func__);
}


void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {

    // TODO
    error_functionNotImplemented(__func__);
}


void mixQureg(Qureg qureg, Qureg other, qreal prob) {

    // TODO:
    // other can be another density matrix, so:
    // qureg ->  (1-prob)qureg + (prob)other
    // 
    // OR a statevector, so that
    // qureg -> (1-prob)qureg + prob|other><other|
    //
    // if the latter is easy enough to distribute (which I suspect it is)

    // TODO
    error_functionNotImplemented(__func__);
}


void mixKrausMap(Qureg qureg, int* qubits, int numQubits, KrausMap map) {

    // TODO
    error_functionNotImplemented(__func__);
}


} // end de-mangler