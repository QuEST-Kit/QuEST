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



/*
 * STRINGS
 */

void reportStr(std::string str) {
    validate_envIsInit(__func__);

    print(str);
    print_newlines();
}

extern "C" void reportStr(const char* str) {
    reportStr(string(str));
}



/*
 * SCALARS
 */

void reportScalar(string label, string numstr) {
    validate_envIsInit(__func__);

    // harmlessly re-validates
    reportStr(label + (label.empty()? "" : ": ") + numstr);
}

void reportScalar(string      label, qcomp num) { reportScalar(label,         printer_toStr(num)); }
void reportScalar(string      label, qreal num) { reportScalar(label,         printer_toStr(num)); }
void reportScalar(const char* label, qreal num) { reportScalar(string(label), printer_toStr(num)); }

extern "C" {
    void reportScalar       (const char* label, qcomp num) { reportScalar(string(label), printer_toStr(num)); }
    void _reportScalar_real (const char* label, qreal num) { reportScalar(string(label), printer_toStr(num)); }
}
