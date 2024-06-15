/** @file
 * API signatures for initialisaing Qureg states.
 */

#ifndef INITIALISATIONS_H
#define INITIALISATIONS_H

#include "quest/include/qureg.h"

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



void initBlankState(Qureg qureg);

void initZeroState(Qureg qureg);

void initPlusState(Qureg qureg);

void initClassicalState(Qureg qureg, qindex stateInd);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // INITIALISATIONS_H