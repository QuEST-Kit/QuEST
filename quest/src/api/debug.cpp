/** @file
 * API definitions for debugging QuEST behaviour, or controlling
 * input validation.
 */

#include "quest/src/core/validation.hpp"
#include "quest/src/core/printer.hpp"


// enable invocation by both C and C++ binaries
extern "C" {


void setValidationOn() {
    validate_enable();
}

void setValidationOff() {
    validate_disable();
}

void setNumReportedItems(qindex num) {
    validate_numReportedItems(num, __func__);

    printer_setMaxNumPrintedItems(num);
}


} // end de-name mangler
