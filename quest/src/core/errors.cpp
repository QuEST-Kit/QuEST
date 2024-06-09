/** @file
 * Defensively designed functions for checking internal preconditions, 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 */

#include "quest/src/comm/communication.hpp"

#include <iostream>
#include <string>



/*
 * INTERNAL ERROR RESPONSE
 */

void raiseInternalError(std::string errorMsg) {

    if (comm_getRank() == 0)
        std::cout 
            << "\n\n"
            << "A fatal internal QuEST error occurred. "
            << errorMsg << " "
            << "Please report this to the developers. QuEST will now exit..."
            << "\n"
            << std::endl;

    exit(EXIT_FAILURE);
}



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(std::string msg, std::string illFormedVar) {

    raiseInternalError("User input validation failed and an error string was attemptedly prepared, but an ill-formed variable was attempted substitution. This variable was \"" + illFormedVar + "\" and was attempted substitution into message:\n" + msg + "\n");
}

void error_validationMessageVarNotSubstituted(std::string msg, std::string var) {

    raiseInternalError("User input validation failed and an error string was attemptedly prepared. However, the internal variable \"" + var + "\" was unable to be found and substituted into the message, which was:\n" + msg + "\n");
}

void error_validationMessageContainedUnsubstitutedVars(std::string msg) {

    raiseInternalError("User input validation failed and an error string was prepared. However, the message contained unexpected (and potentially ill-formed) unsubstituted variables. The message was:\n" + msg + "\n");
}