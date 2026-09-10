#include <quest.h>

#if defined(_OPENMP)
#error "QuEST must not export private OpenMP compilation flags"
#endif

int main(void) {
    initCustomQuESTEnv(0, 0, 0);
    Qureg qureg = createQureg(1);
    initZeroState(qureg);
    qreal probability = calcTotalProb(qureg);
    destroyQureg(qureg);
    finalizeQuESTEnv();
    return probability == (qreal) 1 ? 0 : 1;
}
