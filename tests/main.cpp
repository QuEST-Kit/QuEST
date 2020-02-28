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

#include "QuEST.h"

/** The global QuESTEnv instance, to be created and destroyed once in this 
 * main(), so that the MPI environment is correctly created once when running 
 * distributed unit tests 
 */
QuESTEnv QUEST_ENV;

/** Redefinition of QuEST_validation's invalidQuESTInputError function, called when a 
 * user passes an incorrect parameter (e.g. an negative qubit index). This is 
 * redefined here to, in lieu of printing and exiting, throw a C++ exception
 * which can be caught (and hence unit tested for) by Catch2
 */
 extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
     throw errMsg;
 }
 
/** Explicit declaration of main to create (destroy) the QuESTEnv before (after)
 * invoking the Catch unit tests 
 */
int main( int argc, char* argv[] ) {
  QUEST_ENV = createQuESTEnv();
  int result = Catch::Session().run( argc, argv );
  destroyQuESTEnv(QUEST_ENV);
  return result;
}