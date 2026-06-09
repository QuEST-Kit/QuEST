// End-to-end check against the REAL QuEST library (the actual modified function),
// adapted from the recipe in issue #598. Applies a large all-ones-ish CompMatr to a
// fixed random state and reports amp[0] at high precision. Linked once against the
// default (Kahan) build and once against the QUEST_DENSE_ACCUM_NAIVE build to show
// the modified function's accuracy difference end-to-end.
#include "quest.h"
#include <cstdio>
#include <vector>

int main() {
    initQuESTEnv();
    int numTargets = 12;
    int targets[] = {0,1,2,3,4,5,6,7,8,9,10,11};

    CompMatr matrix = createCompMatr(numTargets);
    // adversarial row 0: large +/- pair that cancels, plus many O(1) terms.
    for (qindex i=0; i<matrix.numRows; i++)
        for (qindex j=0; j<matrix.numRows; j++)
            matrix.cpuElems[i][j] = 1;
    matrix.cpuElems[0][0] = 1e18;
    matrix.cpuElems[0][matrix.numRows-1] = -1e18;
    syncCompMatr(matrix);

    Qureg qureg = createQureg(numTargets);
    // fixed deterministic state: all amps = 1 (unnormalised; validation off)
    std::vector<qcomp> ones(qureg.numAmps, qcomp(1,0));
    setQuregAmps(qureg, 0, ones);

    setQuESTValidationEpsilon(0);
    applyCompMatr(qureg, targets, numTargets, matrix);
    qcomp amp = getQuregAmp(qureg, 0);
    // expected exact amp[0] = (1e18) + (numRows-2)*1 + (-1e18) = numRows-2
    printf("amp0.real = %.20g   (exact = %lld)\n", (double)real(amp), (long long)matrix.numRows-2);

    destroyCompMatr(matrix);
    destroyQureg(qureg);
    finalizeQuESTEnv();
    return 0;
}
