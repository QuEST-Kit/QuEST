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



/*
 * recast QuEST errors into exceptions which Catch2 can intercept
 */

/// @private 
extern "C" void validationErrorHandler(const char* errFunc, const char* errMsg) {

  throw std::runtime_error(std::string(errFunc) + ": " + std::string(errMsg));
}


/** Explicit declaration of main to create (destroy) the QuESTEnv before (after)
 * invoking the Catch unit tests 
 */
int main(int argc, char* argv[]) {

  initQuESTEnv();
  setInputErrorHandler(validationErrorHandler);
  setRandomTestStateSeeds();

  int result = Catch::Session().run( argc, argv );

  finalizeQuESTEnv();
  return result;
}
