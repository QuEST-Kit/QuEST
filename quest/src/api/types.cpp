/** @file
 * API definitions of convenience functions for using QuEST's numerical types
 */

#include "quest/include/types.h"

#include "quest/src/core/formatter.hpp"

#include <iostream>


extern "C" void reportQcomp(qcomp num) {

    std::cout << form_str(num) << std::endl;
}
