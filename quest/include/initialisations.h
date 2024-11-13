/** @file
 * API signatures for initialisaing Qureg states.
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



void initBlankState(Qureg qureg);

void initZeroState(Qureg qureg);

void initPlusState(Qureg qureg);

void initPureState(Qureg qureg, Qureg pure);

void initClassicalState(Qureg qureg, qindex stateInd);

void initDebugState(Qureg qureg);

void initArbitraryState(Qureg qureg, qcomp* amps);

void initRandomPureState(Qureg qureg);



void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);

void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);

void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum);

void setQuregToClone(Qureg targetQureg, Qureg copyQureg);

void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2);

qreal setQuregToRenormalized(Qureg qureg);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // INITIALISATIONS_H