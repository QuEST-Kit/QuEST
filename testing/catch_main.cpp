/** @file
 * This file is left empty so that catch doesn't need 
 * slow (~16s) recompilation each time unit tests are edited
 *
 * @author Tyson Jones
 */

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

/** Redefinition of QuEST_validation's invalidQuESTInputError function, called when a 
 * user passes an incorrect parameter (e.g. an negative qubit index). This is 
 * redefined here to, in lieu of printing and exiting, throw a C++ exception
 * which can be caught (and hence unit tested for) by Catch2
 */
 extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
     throw errMsg;
 } 