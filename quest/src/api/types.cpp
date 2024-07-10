/** @file
 * API definitions of convenience functions for using QuEST's numerical types
 */

#include "quest/include/types.h"

#include "quest/src/core/formatter.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"

#include <iostream>


extern "C" void reportQcomp(qcomp num) {
    validate_envIsInit(__func__);

    // only root node reports (but no synch necesary)
    if (comm_getRank() != 0)
        return;

    std::cout << form_str(num) << std::endl;
}
