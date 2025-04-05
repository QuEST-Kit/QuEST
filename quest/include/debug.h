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



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup debug_seed Seeding
 * @brief Functions for seeding QuEST's random generators.
 * @details Re-seeding with identical seeds will determine all of QuEST's subsequent 
 *          random outputs (such as measurement and random state preparation), and can
 *          be done at any stage of execution. When seeding is not explicitly performed,
 *          QuEST will attempt to use a cryptographically secure pseudorandom number generator
 *          (CSPRNG) if locally available, else fall back to a standard PRNG, via using
 *          the standard C++ `random_device` class.
 * @{
 */


/// @notdoced
void setSeeds(unsigned* seeds, int numSeeds);


/// @notdoced
void setSeedsToDefault();


/// @notdoced
void getSeeds(unsigned* seeds);


/// @notdoced
int getNumSeeds();


/** @} */



/** 
 * @defgroup debug_validation Validation
 * @brief Functions to control QuEST's user-input validation.
 * @details These can be used to adjust the precision with which properties like unitarity 
 *          are checked/enforced, or otherwise disable all input validation (e.g. is the
 *          given qubit index valid?). Note passing erroneous input while validation is 
 *          disabled can result in runtime errors like segmentation faults. 
 * @{
 */


/// @notdoced
void setInputErrorHandler(void (*callback)(const char* func, const char* msg));


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


/** @} */



/** 
 * @defgroup debug_reporting Reporting
 * @brief Functions to control how QuEST's reporters display and truncate information.
 * @{
 */


/// @notdoced
/// @nottested
void setMaxNumReportedItems(qindex numRows, qindex numCols);


/// @notdoced
void setMaxNumReportedSigFigs(int numSigFigs);


/// @notdoced
void setNumReportedNewlines(int numNewlines);


/** @} */



/** 
 * @defgroup debug_cache Caching
 * @brief Functions to control temporary memory used by the QuEST process.
 * @{
 */


/// @notdoced
qindex getGpuCacheSize();


/// @notdoced
void clearGpuCache();


/** @} */



/** 
 * @defgroup debug_info Info
 * @brief Functions for getting debugging information.
 * @{
 */


/// @notdoced
/// @nottested
void getEnvironmentString(char str[200]);


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


/// @ingroup debug_seed
/// @nottested
/// @notdoced
/// @cpponly
void setSeeds(std::vector<unsigned> seeds);


/// @ingroup debug_seed
/// @nottested
/// @notdoced
/// @cpponly
std::vector<unsigned> getSeeds();


#endif // __cplusplus



#endif // DEBUG_H

/** @} */ // (end file-wide doxygen defgroup)
