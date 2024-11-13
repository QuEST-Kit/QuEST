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

// enable invocation by both C and C++ binaries
extern "C" {



/*
 * SEEDING
 */


void setSeeds(unsigned* seeds, int numSeeds) {

    rand_setSeeds(std::vector<unsigned>(seeds, seeds+numSeeds));
}

void setSeedsToDefault() {

    rand_setSeedsToDefault();
}


int getNumSeeds() {

    return rand_getNumSeeds();
}

void getSeeds(unsigned* seeds) {

    auto vec = rand_getSeeds();
    auto num = rand_getNumSeeds();

    for (int i=0; i<num; i++)
        seeds[i] = vec[i];
}



/*
 * VALIDATION
 */


void setValidationOn() {
    validateconfig_enable();
}

void setValidationOff() {
    validateconfig_disable();
}


void setValidationEpsilon(qreal eps) {
    validate_newEpsilonValue(eps, __func__);

    validateconfig_setEpsilon(eps);
}

void setValidationEpsilonToDefault() {

    validateconfig_setEpsilonToDefault();
}

qreal getValidationEpsilon() {
    return validateconfig_getEpsilon();
}



/*
 * REPORTERS
 */


void setMaxNumReportedItems(qindex numRows, qindex numCols) {
    validate_newMaxNumReportedScalars(numRows, numCols, __func__);

    // replace 0 values (indicating no truncation) with max-val,
    // since there can never be max(qindex)-many amps
    qindex max = std::numeric_limits<qindex>::max();
    numRows = (numRows == 0)? max : numRows;
    numCols = (numCols == 0)? max : numCols;

    printer_setMaxNumPrintedScalars(numRows, numCols);
}



/*
 * GPU CACHE
 */


qindex getGpuCacheSize() {

    if (getQuESTEnv().isGpuAccelerated)
        return gpu_getCacheMemoryInBytes();

    // safely returns 0 if not GPU accelerated
    return 0;
}


void clearGpuCache() {

    // safely do nothing if not GPU accelerated
    if (getQuESTEnv().isGpuAccelerated)
        gpu_clearCache();
}



} // end de-name mangler
