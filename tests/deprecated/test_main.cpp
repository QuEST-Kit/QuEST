/** @file
 * This file is left mostly empty so that catch doesn't need 
 * slow (~16s) recompilation each time unit tests are edited
 *
 * @author Tyson Jones
 */


/** Use our modified Catch in custom-main mode (main defined below).
 * catch.hpp was modified to, in distributed mode, output only once.
 */
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"


#define INCLUDE_DEPRECATED_FUNCTIONS 1
#define DISABLE_DEPRECATION_WARNINGS 1

#include "quest.h"
#include "test_utilities.hpp"

#include <stdexcept>


/** Redefinition of QuEST_validation's invalidQuESTInputError function, called when a 
 * user passes an incorrect parameter (e.g. a negative qubit index). This is 
 * redefined here to, in lieu of printing and exiting, throw a C++ exception
 * which can be caught (and hence unit tested for) by Catch2
 */
 extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
  
    throw std::runtime_error(errMsg);
 }


/** Explicit declaration of main to create (destroy) the QuESTEnv before (after)
 * invoking the Catch unit tests 
 */
int main(int argc, char* argv[]) {

  initQuESTEnv();
  setRandomTestStateSeeds();

  int result = Catch::Session().run( argc, argv );

  finalizeQuESTEnv();
  return result;
}
