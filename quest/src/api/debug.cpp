/** @file
 * API definitions for debugging QuEST behaviour, or controlling
 * input validation.
 */

#include "quest/src/core/validation.hpp"


// enable invocation by both C and C++ binaries
extern "C" {


void setValidationOn() {
    validate_enable();
}

void setValidationOff() {
    validate_disable();
}


} // end de-name mangler
