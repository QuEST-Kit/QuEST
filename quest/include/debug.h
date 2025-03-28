/** @file
 * API signatures for debugging QuEST behaviour, 
 * controlling input validation, changing reporter
 * parameters or seeding random generation.
 * 
 * @author Tyson Jones
 *
 * @defgroup debug Debug
 * @ingroup api
 * @brief Utilities for controlling QuEST behaviour such as seeding, input validation and printing.
 * @{
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

/// @notdoced
void setSeeds(unsigned* seeds, int numSeeds);

/// @notdoced
void setSeedsToDefault();

/// @notdoced
void getSeeds(unsigned* seeds);

/// @notdoced
int getNumSeeds();


/*
 * VALIDATION
 */

/// @notdoced
void invalidQuESTInputError(const char* msg, const char* func);

/// @notdoced
void setValidationOn();

/// @notdoced
void setValidationOff();

/// @notdoced
void setValidationEpsilonToDefault();

/// @notdoced
void setValidationEpsilon(qreal eps);

/// @notdoced
qreal getValidationEpsilon();


/*
 * REPORTING
 */

/// @notdoced
/// @nottested
void setMaxNumReportedItems(qindex numRows, qindex numCols);

/// @notdoced
void setMaxNumReportedSigFigs(int numSigFigs);

/// @notdoced
void setNumReportedNewlines(int numNewlines);


/*
 * CACHING
 */

/// @notdoced
qindex getGpuCacheSize();

/// @notdoced
void clearGpuCache();


/*
 * ENVIRONMENT
 */

/// @notdoced
/// @nottested
void getEnvironmentString(char str[200]);


// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // DEBUG_H

/** @} (end doxygen defgroup) */
