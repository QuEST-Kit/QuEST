/** @file
 * A minimum C/C++-agnostic example of running
 * QuEST, reporting the execution environment
 * and preparing 20-qubit random pure state.
 * 
 * @author Tyson Jones
*/

#include "quest.h"

int main(void) {

    initQuESTEnv();
    reportQuESTEnv();

    Qureg qureg = createForcedQureg(20);
    reportQuregParams(qureg);

    initRandomPureState(qureg);
    reportQureg(qureg);

    qreal prob = calcTotalProb(qureg);
    reportScalar("Total probability", prob);
    
    destroyQureg(qureg);
    finalizeQuESTEnv();

    return 0;
}
