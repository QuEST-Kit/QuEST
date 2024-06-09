/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#include "quest/src/core/errors.hpp"

#include <iostream>
#include <cstdlib>
#include <map>



/*
 * INVALID INPUT ERROR MESSAGES
 * which can contain variables with syntax ${VAR1} ${VAR2}, substituted at error-throw with
 * runtime parameters via assertThat(..., {{"${VAR1}",1}, {"${VAR2}",2}}, ...)
 */

namespace report {

    std::string EXAMPLE = 
        "Example variables were ${X} and ${MY_FLAG}";

}



/*
 * INVALID INPUT RESPONSE BEHAVIOUR
 */

// default C/C++ compatible error response is to simply exit in fail state
void default_invalidQuESTInputError(const char* msg, const char* func) {

    if (comm_getRank() == 0)
        std::cout 
            << std::endl
            << "QuEST encountered a validation error during function '" << func << "':\n"
            << msg << "\nExiting..." 
            << std::endl;

    exit(EXIT_FAILURE);
}

// enable default error response to be user-overriden as a weak symbol (even in C, and on Windows)
extern "C" {

    #ifndef _WIN32
        #pragma weak invalidQuESTInputError
        void invalidQuESTInputError(const char* msg, const char* func) {
            default_invalidQuESTInputError(msg, func);
        }
    #elif defined(_WIN64)
        #pragma comment(linker, "/alternatename:invalidQuESTInputError=default_invalidQuESTInputError")   
    #else
        #pragma comment(linker, "/alternatename:_invalidQuESTInputError=_default_invalidQuESTInputError")
    #endif

} // end C++ de-mangler



/*
 * UTILITIES
 */

std::string getStringWithSubstitutedVars(std::string oldStr, std::initializer_list<std::map<std::string, int>::value_type> varsAndVals) {

    std::string newStr = oldStr;

    for (auto varAndVal : varsAndVals) {

        std::string var = std::get<0>(varAndVal);
        std::string val = std::to_string(std::get<1>(varAndVal));

        if (var[0] != '$' || var[1] != '{' || var.back() != '}' )
            error_validationMessageVarWasIllFormed(newStr, var);

        size_t ind = newStr.find(var);
        if (ind == std::string::npos)
            error_validationMessageVarNotSubstituted(newStr, var);

        newStr = newStr.replace(ind, var.length(), val);
    }

    // assert there is no $ left in the strings
    if (newStr.find("$") != std::string::npos)
        error_validationMessageContainedUnsubstitutedVars(newStr);

    return newStr;
}

void assertThat(bool valid, std::string msg, const char* func) {

    if (!valid)
        invalidQuESTInputError(msg.c_str(), func);
}

void assertThat(bool valid, std::string msg, std::initializer_list<std::map<std::string, int>::value_type> vars, const char* func) {

    std::string newMsg = getStringWithSubstitutedVars(msg, vars);
    assertThat(valid, newMsg, func);
}
