/** @file
 * API convenience functions related to types.
 * 
 * @author Tyson Jones
 */

#include "quest/include/types.h"

#include "quest/src/core/printer.hpp"
#include "quest/src/core/validation.hpp"



extern "C" void reportQcomp(qcomp num) {
    validate_envIsInit(__func__);

    print(num);
}
