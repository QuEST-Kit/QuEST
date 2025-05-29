#include "quest.h"

int main() {
    initQuESTEnv();

    Qureg qureg = createQureg(3);
    PauliStr str = getInlinePauliStr("XYZ", {0,1,2});
    qcomp angle = .4 + .8i;

    initPlusState(qureg);
    applyNonUnitaryPauliGadget(qureg, str, angle);

    qreal norm = calcTotalProb(qureg);
    reportScalar("norm", norm);

    finalizeQuESTEnv();
    return 0;
}