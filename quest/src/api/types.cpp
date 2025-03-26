/** @file
 * API convenience functions related to types.
 * 
 * @author Tyson Jones
 */

#include "quest/include/types.h"

#include "quest/src/core/printer.hpp"
#include "quest/src/core/validation.hpp"

#include <string>

using std::string;


void reportScalar(string label, string num) {
    validate_envIsInit(__func__);

    print(label + (label.empty()? "" : ": ") + num);
}

void reportScalar(string      label, qcomp num) { reportScalar(label,         printer_toStr(num)); }
void reportScalar(string      label, qreal num) { reportScalar(label,         printer_toStr(num)); }
void reportScalar(const char* label, qreal num) { reportScalar(string(label), printer_toStr(num)); }

extern "C" {
    void reportScalar       (const char* label, qcomp num) { reportScalar(string(label), printer_toStr(num)); }
    void _reportScalar_real (const char* label, qreal num) { reportScalar(string(label), printer_toStr(num)); }
}
