/** @file
 * API definitions for debugging QuEST behaviour, 
 * controlling input validation, changing reporter
 * parameters or seeding random generation.
 * 
 * @author Tyson Jones
 */

#include "quest/include/types.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/randomiser.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <vector>
#include <limits>

using std::vector;



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {


/*
 * SEEDING
 */


void setQuESTSeeds(unsigned* seeds, int numSeeds) {
    validate_envIsInit(__func__);
    validate_randomSeeds(seeds, numSeeds, __func__);

    // consults only root-node seeds
    rand_setSeeds(vector<unsigned>(seeds, seeds+numSeeds));
}

void setQuESTSeedsToDefault() {
    validate_envIsInit(__func__);

    rand_setSeedsToDefault();
}


int getQuESTNumSeeds() {
    validate_envIsInit(__func__);

    return rand_getNumSeeds();
}

void getQuESTSeeds(unsigned* seeds) {
    validate_envIsInit(__func__);

    auto vec = rand_getSeeds();
    auto num = rand_getNumSeeds();

    for (int i=0; i<num; i++)
        seeds[i] = vec[i];
}



/*
 * VALIDATION
 */

void setQuESTInputErrorHandler(void (*callback)(const char*, const char*)) {
    validate_envIsInit(__func__);

    validateconfig_setErrorHandler(callback);
}

void setQuESTValidationOn() {
    validate_envIsInit(__func__);
    
    validateconfig_enable();
}

void setQuESTValidationOff() {
    validate_envIsInit(__func__);

    // disables all validation and computation
    // of matrix properties like isUnitary. Also
    // means pre-computed matrix properties are
    // ignored. It does not however erase pre-
    // computed properties; subsequently restoring 
    // validation will not necessitate re-eval.

    validateconfig_disable();
}


void setQuESTValidationEpsilon(qreal eps) {
    validate_envIsInit(__func__);
    validate_newEpsilonValue(eps, __func__);

    validateconfig_setEpsilon(eps);
    util_setEpsilonSensitiveHeapFlagsToUnknown();
}

void setQuESTValidationEpsilonToDefault() {
    validate_envIsInit(__func__);

    validateconfig_setEpsilonToDefault();
    util_setEpsilonSensitiveHeapFlagsToUnknown();
}

qreal getQuESTValidationEpsilon() {
    validate_envIsInit(__func__);

    return validateconfig_getEpsilon();
}



/*
 * REPORTER CONFIGURATION
 */


void setQuESTMaxNumReportedItems(qindex numRows, qindex numCols) {
    validate_envIsInit(__func__);
    validate_newMaxNumReportedScalars(numRows, numCols, __func__);

    // replace 0 values (indicating no truncation) with max-val,
    // since there can never be max(qindex)-many amps
    qindex max = std::numeric_limits<qindex>::max();
    numRows = (numRows == 0)? max : numRows;
    numCols = (numCols == 0)? max : numCols;

    printer_setMaxNumPrintedScalars(numRows, numCols);
}


void setQuESTMaxNumReportedSigFigs(int numSigFigs) {
    validate_envIsInit(__func__);
    validate_newMaxNumReportedSigFigs(numSigFigs, __func__);

    printer_setMaxNumPrintedSigFig(numSigFigs);
}


void setQuESTNumReportedNewlines(int numNewlines) {
    validate_envIsInit(__func__);
    validate_newNumReportedNewlines(numNewlines, __func__);

    printer_setNumTrailingNewlines(numNewlines);
}


void setQuESTReportedPauliChars(const char* paulis) {
    validate_envIsInit(__func__);
    validate_numPauliChars(paulis, __func__);

    printer_setPauliChars(paulis);
}


void setQuESTReportedPauliStrStyle(int flag) {
    validate_envIsInit(__func__);
    validate_reportedPauliStrStyleFlag(flag, __func__);

    printer_setPauliStrFormat(flag);
}



/*
 * GPU CACHE
 */


qindex getQuESTGpuCacheSize() {
    validate_envIsInit(__func__);

    if (getQuESTEnv().isGpuAccelerated)
        return gpu_getCacheMemoryInBytes();

    // safely returns 0 if not GPU accelerated
    return 0;
}


void clearQuESTGpuCache() {
    validate_envIsInit(__func__);

    // safely do nothing if not GPU accelerated
    if (getQuESTEnv().isGpuAccelerated)
        gpu_clearCache();
}


} // end de-name mangler



/*
 * C++ OVERLOADS
 */


void setQuESTSeeds(vector<unsigned> seeds) {
    setQuESTSeeds(seeds.data(), seeds.size());
}

vector<unsigned> getQuESTSeeds() {
    validate_envIsInit(__func__);

    // allocate temp vector, and pedantically validate successful
    vector<unsigned> out;
    int numSeeds = rand_getNumSeeds();
    auto callback = [&]() { validate_tempListAllocSucceeded(false, numSeeds, sizeof(unsigned), __func__); };
    util_tryAllocVector(out, numSeeds, callback);

    getQuESTSeeds(out.data());
    return out;
}
