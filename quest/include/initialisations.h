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


/// @notyetdoced
/// @notyetvalidated
void initBlankState(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
void initZeroState(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
void initPlusState(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
/// @notyettested
void initPureState(Qureg qureg, Qureg pure);


/// @notyetdoced
/// @notyetvalidated
void initClassicalState(Qureg qureg, qindex stateInd);


/// @notyetdoced
/// @notyetvalidated
void initDebugState(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
void initArbitraryPureState(Qureg qureg, qcomp* amps);


/// @notyetdoced
/// @notyetvalidated
void initRandomPureState(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
void initRandomMixedState(Qureg qureg, qindex numPureStates);


/** @} */



/** 
 * @defgroup init_amps Amplitudes
 * @brief Functions for overwriting Qureg amplitudes.
 * @{
 */


/// @notyetdoced
/// @notyetvalidated
void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);


/// @notyetdoced
/// @notyetvalidated
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);


/// @notyetdoced
/// @notyetvalidated
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);


/// @notyetdoced
/// @notyettested
void setQuregToClone(Qureg targetQureg, Qureg copyQureg);


/// @notyetdoced
/// @notyettested
void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2);


/// @notyetdoced
/// @notyetvalidated
qreal setQuregToRenormalized(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum);


/// @notyetdoced
/// @notyettested
void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);


/// @notyetdoced
/// @notyettested
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
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
void setQuregAmps(Qureg qureg, qindex startInd, std::vector<qcomp> amps);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, std::vector<std::vector<qcomp>> amps);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, std::vector<qcomp> amps);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
void setQuregToPartialTrace(Qureg out, Qureg in, std::vector<int> traceOutQubits);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
void setQuregToReducedDensityMatrix(Qureg out, Qureg in, std::vector<int> retainQubits);


#endif // __cplusplus



#endif // INITIALISATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
