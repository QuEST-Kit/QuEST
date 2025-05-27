/** @file
 * Shows the output of reportQuESTEnv() on
 * Github Action runners
 * 
 * @author Tyson Jones
*/

#include "quest.h"

int main() {
    initQuESTEnv();
    reportQuESTEnv();
    finalizeQuESTEnv();
}
