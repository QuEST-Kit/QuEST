/** @file
 * Examples of setting the invalid error input handler in C++14,
 * overriding the default behaviour of immediately exiting.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>


void myErrorHandlerA(const char* errFunc, const char* errMsg) {

    std::string func(errFunc);
    std::string msg(errMsg);

    std::cout 
        << "an error?? in the '" << errFunc << "' function?? with message '" << msg << "'! "
        << "how queer!! ive never seen such a thing - i must throw an exception post-haste!!!"
        << std::endl
        << std::endl;

    // exception forces control-flow out of QuEST env, safely back to user.
    // without this, control-flow would return to the QuEST backend and likely
    // cause internal integrity checks to fail, or segmentation faults
    throw std::runtime_error(std::string(errFunc) + ": " + std::string(errMsg));
}


void myErrorHandlerB(const char* errFunc, const char* errMsg) {

    std::string func(errFunc);
    std::string msg(errMsg);

    std::cout 
        << "Function '" << func << "' threw '" << msg << "'." 
        << std::endl
        << "We will now exit completely. Good day!" 
        << std::endl
        << std::endl;

    exit(0);
}


int main() {
    initQuESTEnv();

    setInputErrorHandler(myErrorHandlerA);

    try {
        Qureg qureg = createQureg(-123);
    } catch (std::runtime_error& e) {
        std::cout 
            << "Error caught! Function aborted, but execution continues." 
            << std::endl
            << std::endl;
    }

    setInputErrorHandler(myErrorHandlerB);
    initQuESTEnv(); // illegal to recall
    
    std::cout << "this will never be reached, because myErrorHandlerB exits!" << std::endl;

    finalizeQuESTEnv();
}
