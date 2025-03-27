/** @file
 * Entry-point for the ported tests of QuEST's deprecated
 * v3 API.
 * 
 * This file was originally written for catch2 v2, though has
 * since been refactored for compatibility with catch2 v3. The
 * comments however have not been updated and may mislead.
 *
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */


/** Use our modified Catch in custom-main mode (main defined below).
 * catch.hpp was modified to, in distributed mode, output only once.
 */
#include <catch2/catch_session.hpp>


#define INCLUDE_DEPRECATED_FUNCTIONS 1
#define DISABLE_DEPRECATION_WARNINGS 1

#include "quest/include/quest.h"
#include "test_utilities.hpp"

#include <stdexcept>


/* Redefinition of QuEST_validation's invalidQuESTInputError function, called when a 
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
