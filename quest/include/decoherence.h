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


/// @notdoced
/// @notvalidated
void mixDephasing(Qureg qureg, int qubit, qreal prob);

/// @notdoced
/// @notvalidated
void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob);

/// @notdoced
/// @notvalidated
void mixDepolarising(Qureg qureg, int qubit, qreal prob);

/// @notdoced
/// @notvalidated
void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob);

/// @notdoced
/// @notvalidated
void mixDamping(Qureg qureg, int qubit, qreal prob);

/// @notdoced
/// @notvalidated
void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ);

/// @notdoced
/// @notvalidated
void mixQureg(Qureg qureg, Qureg other, qreal prob);

/// @notdoced
/// @notvalidated
void mixKrausMap(Qureg qureg, int* qubits, int numQubits, KrausMap map);

/// @notdoced
/// @notvalidated
void mixSuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop);


// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // DECOHERENCE_H

/** @} (end doxygen defgroup) */
