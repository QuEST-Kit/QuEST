/** @file
 * API definitions for debugging QuEST behaviour, or controlling
 * input validation.
 */

#include "quest/src/core/validation.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/gpu/gpu_config.hpp"

// enable invocation by both C and C++ binaries
extern "C" {



/*
 * VALIDATION
 */


void setValidationOn() {
    validate_enable();
}

void setValidationOff() {
    validate_disable();
}



/*
 * REPORTERS
 */


void setNumReportedItems(qindex num) {
    validate_numReportedItems(num, __func__);

    printer_setMaxNumPrintedItems(num);
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
