/** @file
 * API signatures for debugging QuEST behaviour, or controlling
 * input validation.
 */

#ifndef DEBUG_H
#define DEBUG_H

#include "quest/include/types.h"

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


/*
 * SEEDING
 */

void setSeeds(unsigned* seeds, int numSeeds);
void setSeedsToDefault();

void getSeeds(unsigned* seeds);
int getNumSeeds();


/*
 * VALIDATION
 */

void invalidQuESTInputError(const char* msg, const char* func);

void setValidationOn();
void setValidationOff();

void setValidationEpsilonToDefault();
void setValidationEpsilon(qreal eps);
qreal getValidationEpsilon();


/*
 * REPORTING
 */

void setMaxNumReportedItems(qindex numRows, qindex numCols);
void setMaxNumReportedSigFigs(int numSigFigs);


/*
 * CACHING
 */

qindex getGpuCacheSize();
void clearGpuCache();


/*
 * ENVIRONMENT
 */

void getEnvironmentString(char str[200]);


// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // DEBUG_H