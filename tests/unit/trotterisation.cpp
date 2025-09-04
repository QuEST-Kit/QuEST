/** @file
 * Unit tests of the trotterisation module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unittrotter Trotterisation
 * @ingroup unittests
 */

#include "quest.h"



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[trotterisation]"



/**
 * @todo
 * UNTESTED FUNCTIONS
 */

void applyTrotterizedNonUnitaryPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps);

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);

void applyTrotterizedControlledPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps);

void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps);

void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps);

void applyTrotterizedUnitaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal time, int order, int reps);

void applyTrotterizedImaginaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal tau, int order, int reps);

void applyTrotterizedNoisyTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal* damps, PauliStr* jumps, int numJumps, qreal time, int order, int reps);
