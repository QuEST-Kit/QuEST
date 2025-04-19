/** @file
 * Examples of setting the invalid error input handler in C11,
 * overriding the default behaviour of immediately exiting.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <stdio.h>
#include <stdlib.h>


void myErrorHandler(const char* errFunc, const char* errMsg) {
    printf("Ruh-roh, Raggy! Function '%s' has reported '%s'.\n", errFunc, errMsg);
    printf("We will now be very good children and exit immediately!\n");
    exit(0);
}


int main() {
    initQuESTEnv();
    setInputErrorHandler(myErrorHandler);

    Qureg qureg = createQureg(-123);

    finalizeQuESTEnv();
}
