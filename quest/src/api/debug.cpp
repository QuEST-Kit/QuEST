/** @file
 * API definitions for debugging QuEST behaviour, or controlling
 * input validation.
 */

#include "quest/src/core/validation.hpp"
#include "quest/src/core/formatter.hpp"


// enable invocation by both C and C++ binaries
extern "C" {


void setValidationOn() {
    validate_enable();
}

void setValidationOff() {
    validate_disable();
}

void setNumReportedMatrixElems(qindex num) {
    validate_numReportedMatrixElems(num, __func__);

    form_setMaxNumPrintedMatrixElems(num);
}


} // end de-name mangler
