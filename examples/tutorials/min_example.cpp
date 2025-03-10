/** @file
 * A minimum C/C++-agnostic example of running
 * QuEST, reporting the execution environment
 * and preparing 20-qubit plus state.
 * 
 * @author Tyson Jones
*/

#include "quest.h"

int main(void) {
  initQuESTEnv();
  reportQuESTEnv();

  Qureg qureg = createQureg(20);
  reportQuregParams(qureg);

  initPlusState(qureg);
  reportQureg(qureg);

  destroyQureg(qureg);
  finalizeQuESTEnv();

  return 0;
}