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

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


void mixDephasing(Qureg qureg, int qubit, qreal prob);

void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob);

void mixDepolarising(Qureg qureg, int qubit, qreal prob);

void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob);

void mixDamping(Qureg qureg, int qubit, qreal prob);

void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ);

void mixQureg(Qureg qureg, Qureg other, qreal prob);

void mixKrausMap(Qureg qureg, int* qubits, int numQubits, KrausMap map);


// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // DECOHERENCE_H

/** @} (end doxygen defgroup) */
