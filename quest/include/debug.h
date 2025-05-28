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


/// @notyetdoced
void setSeeds(unsigned* seeds, int numSeeds);


/// @notyetdoced
void setSeedsToDefault();


/// @notyetdoced
void getSeeds(unsigned* seeds);


/// @notyetdoced
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


/** @notyetdoced
 *
 * @see
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/setting_errorhandler.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/setting_errorhandler.cpp) examples
 */
void setInputErrorHandler(void (*callback)(const char* func, const char* msg));


/// @notyetdoced
void setValidationOn();


/// @notyetdoced
void setValidationOff();


/// @notyetdoced
void setValidationEpsilonToDefault();


/// @notyetdoced
void setValidationEpsilon(qreal eps);


/// @notyetdoced
qreal getValidationEpsilon();


/** @} */



/** 
 * @defgroup debug_reporting Reporting
 * @brief Functions to control how QuEST's reporters display and truncate information.
 * @{
 */


/// @notyetdoced
/// @notyettested
void setMaxNumReportedItems(qindex numRows, qindex numCols);


/** @notyetdoced
 * > This function does not affect the significant figures in printed memory sizes
 * > (e.g. `5.32 KiB`) which is always shown with three significant figures 
 * > (or four when in bytes, e.g. `1023 bytes`).
 */
void setMaxNumReportedSigFigs(int numSigFigs);


/// @notyetdoced
void setNumReportedNewlines(int numNewlines);


/** 
 * @notyetdoced
 * @notyettested
 * @myexample
 * ```
   PauliStr str = getInlinePauliStr("XYZ", {0,10,20});
   reportPauliStr(str);

   setReportedPauliChars(".xyz");
   reportPauliStr(str);
 * ```
 */
void setReportedPauliChars(const char* paulis);


/** 
 * @notyetdoced
 * @notyettested
 * @myexample
 * ```
   PauliStr str = getInlinePauliStr("XYZ", {0,10,20});

   setReportedPauliStrStyle(0);
   reportPauliStr(str);

   setReportedPauliStrStyle(1);
   reportPauliStr(str);
 * ```
 */
void setReportedPauliStrStyle(int style);


/** @} */



/** 
 * @defgroup debug_cache Caching
 * @brief Functions to control temporary memory used by the QuEST process.
 * @{
 */


/// @notyetdoced
qindex getGpuCacheSize();


/// @notyetdoced
void clearGpuCache();


/** @} */



/** 
 * @defgroup debug_info Info
 * @brief Functions for getting debugging information.
 * @{
 */


/// @notyetdoced
/// @notyettested
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
/// @notyettested
/// @notyetdoced
/// @cppvectoroverload
/// @see setSeeds()
void setSeeds(std::vector<unsigned> seeds);


/// @ingroup debug_seed
/// @notyettested
/// @notyetdoced
/// @cpponly
/// @see getSeeds()
std::vector<unsigned> getSeeds();


#endif // __cplusplus



#endif // DEBUG_H

/** @} */ // (end file-wide doxygen defgroup)
