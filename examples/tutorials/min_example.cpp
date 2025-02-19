#include "quest/include/quest.h"

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