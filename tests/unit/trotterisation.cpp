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

void applyNonUnitaryTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps);

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps);

void applyControlledTrotterizedPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps);

void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps);

void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps);
