/** @file
 * Defensively designed functions for checking internal preconditions, 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 */

#ifndef ERRORS_HPP
#define ERRORS_HPP

#include <string>



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(std::string msg, std::string illFormedVar);

void error_validationMessageVarNotSubstituted(std::string msg, std::string var);

void error_validationMessageContainedUnsubstitutedVars(std::string msg);



/*
 * COMMUNICATION ERRORS
 */

void error_commAlreadyInit();



/*
 * GPU ERRORS
 */

void error_gpuQueriedButGpuNotCompiled();



#endif // ERRORS_HPP