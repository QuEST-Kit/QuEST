#include <quest.h>

#if defined(_OPENMP)
#error "QuEST must not export private OpenMP compilation flags"
#endif

#if defined(_MSVC_LANG)
#if _MSVC_LANG != 201402L
#error "QuEST must not raise the consumer's C++14 requirement"
#endif
#elif __cplusplus != 201402L
#error "QuEST must not raise the consumer's C++14 requirement"
#endif

int main() {
    initCustomQuESTEnv(0, 0, 0);
    auto qureg = createQureg(1);
    initZeroState(qureg);
    qreal probability = calcTotalProb(qureg);
    destroyQureg(qureg);
    finalizeQuESTEnv();
    return probability == qreal{1} ? 0 : 1;
}
