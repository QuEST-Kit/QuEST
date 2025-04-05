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



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup init_states States
 * @brief Functions for initialising Qureg into physical states.
 * @{
 */


/// @notdoced
/// @notvalidated
void initBlankState(Qureg qureg);


/// @notdoced
/// @notvalidated
void initZeroState(Qureg qureg);


/// @notdoced
/// @notvalidated
void initPlusState(Qureg qureg);


/// @notdoced
/// @notvalidated
/// @nottested
void initPureState(Qureg qureg, Qureg pure);


/// @notdoced
/// @notvalidated
void initClassicalState(Qureg qureg, qindex stateInd);


/// @notdoced
/// @notvalidated
void initDebugState(Qureg qureg);


/// @notdoced
/// @notvalidated
void initArbitraryPureState(Qureg qureg, qcomp* amps);


/// @notdoced
/// @notvalidated
void initRandomPureState(Qureg qureg);


/// @notdoced
/// @notvalidated
void initRandomMixedState(Qureg qureg, qindex numPureStates);


/** @} */



/** 
 * @defgroup init_amps Amplitudes
 * @brief Functions for overwriting Qureg amplitudes.
 * @{
 */


/// @notdoced
/// @notvalidated
void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);


/// @notdoced
/// @notvalidated
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);


/// @notdoced
/// @notvalidated
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);


/// @notdoced
/// @nottested
void setQuregToClone(Qureg targetQureg, Qureg copyQureg);


/// @notdoced
/// @nottested
void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2);


/// @notdoced
/// @notvalidated
qreal setQuregToRenormalized(Qureg qureg);


/// @notdoced
/// @notvalidated
void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum);


/// @notdoced
/// @nottested
void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);


/// @notdoced
/// @nottested
void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits);


/** @} */



// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We 
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

#include <vector>


/// @ingroup init_amps
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void setQuregAmps(Qureg qureg, qindex startInd, std::vector<qcomp> amps);


/// @ingroup init_amps
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, std::vector<std::vector<qcomp>> amps);


/// @ingroup init_amps
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, std::vector<qcomp> amps);


/// @ingroup init_amps
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void setQuregToPartialTrace(Qureg out, Qureg in, std::vector<int> traceOutQubits);


/// @ingroup init_amps
/// @nottested
/// @notdoced
/// @notvalidated
/// @cpponly
void setQuregToReducedDensityMatrix(Qureg out, Qureg in, std::vector<int> retainQubits);


#endif // __cplusplus



#endif // INITIALISATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
