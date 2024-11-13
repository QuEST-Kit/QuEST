/** @file
 * API definitions for debugging QuEST behaviour, or controlling
 * input validation.
 */

#include "quest/include/types.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/randomiser.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <vector>
#include <limits>

// enable invocation by both C and C++ binaries
extern "C" {



/*
 * SEEDING
 */


void setSeeds(unsigned* seeds, int numSeeds) {
    validate_envIsInit(__func__);
    validate_randomSeeds(seeds, numSeeds, __func__);

    rand_setSeeds(std::vector<unsigned>(seeds, seeds+numSeeds));
}

void setSeedsToDefault() {
    validate_envIsInit(__func__);

    rand_setSeedsToDefault();
}


int getNumSeeds() {
    validate_envIsInit(__func__);

    return rand_getNumSeeds();
}

void getSeeds(unsigned* seeds) {
    validate_envIsInit(__func__);

    auto vec = rand_getSeeds();
    auto num = rand_getNumSeeds();

    for (int i=0; i<num; i++)
        seeds[i] = vec[i];
}



/*
 * VALIDATION
 */


void setValidationOn() {
    validate_envIsInit(__func__);
    
    validateconfig_enable();
}

void setValidationOff() {
    validate_envIsInit(__func__);

    validateconfig_disable();
}


void setValidationEpsilon(qreal eps) {
    validate_envIsInit(__func__);
    validate_newEpsilonValue(eps, __func__);

    validateconfig_setEpsilon(eps);
}

void setValidationEpsilonToDefault() {
    validate_envIsInit(__func__);

    validateconfig_setEpsilonToDefault();
}

qreal getValidationEpsilon() {
    validate_envIsInit(__func__);

    return validateconfig_getEpsilon();
}



/*
 * REPORTER CONFIGURATION
 */


void setMaxNumReportedItems(qindex numRows, qindex numCols) {
    validate_envIsInit(__func__);
    validate_newMaxNumReportedScalars(numRows, numCols, __func__);

    // replace 0 values (indicating no truncation) with max-val,
    // since there can never be max(qindex)-many amps
    qindex max = std::numeric_limits<qindex>::max();
    numRows = (numRows == 0)? max : numRows;
    numCols = (numCols == 0)? max : numCols;

    printer_setMaxNumPrintedScalars(numRows, numCols);
}


void setMaxNumReportedSigFigs(int numSigFigs) {
    validate_envIsInit(__func__);
    validate_newMaxNumReportedSigFigs(numSigFigs, __func__);

    printer_setMaxNumPrintedSigFig(numSigFigs);
} 



/*
 * GPU CACHE
 */


qindex getGpuCacheSize() {
    validate_envIsInit(__func__);

    if (getQuESTEnv().isGpuAccelerated)
        return gpu_getCacheMemoryInBytes();

    // safely returns 0 if not GPU accelerated
    return 0;
}


void clearGpuCache() {
    validate_envIsInit(__func__);

    // safely do nothing if not GPU accelerated
    if (getQuESTEnv().isGpuAccelerated)
        gpu_clearCache();
}



} // end de-name mangler
