/** @file
 * API signatures for initialising Quregs into 
 * particular states. Note when a Qureg is GPU-
 * accelerated, these functions only update the
 * state in GPU memory; the CPU amps are unchanged.
 * 
 * @author Tyson Jones
 * 
 * @defgroup initialisations Initialisations
 * @ingroup api
 * @brief Functions for preparing Quregs in particular states.
 * @{
 */

#ifndef INITIALISATIONS_H
#define INITIALISATIONS_H

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/// @notvalidated
void initBlankState(Qureg qureg);

/// @notvalidated
void initZeroState(Qureg qureg);

/// @notvalidated
void initPlusState(Qureg qureg);

/// @notvalidated
/// @nottested
void initPureState(Qureg qureg, Qureg pure);

/// @notvalidated
void initClassicalState(Qureg qureg, qindex stateInd);

/// @notvalidated
void initDebugState(Qureg qureg);

/// @notvalidated
void initArbitraryPureState(Qureg qureg, qcomp* amps);

/// @notvalidated
void initRandomPureState(Qureg qureg);

/// @notvalidated
void initRandomMixedState(Qureg qureg, qindex numPureStates);



/// @notvalidated
void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);

/// @notvalidated
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);

/// @notvalidated
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);

/// @nottested
void setQuregToClone(Qureg targetQureg, Qureg copyQureg);

/// @nottested
void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2);

/// @notvalidated
qreal setQuregToRenormalized(Qureg qureg);

/// @notvalidated
void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum);

/// @nottested
void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);

/// @nottested
void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // INITIALISATIONS_H

/** @} (end doxygen defgroup) */
