#include "quest.h"

int main() {
    initQuESTEnv();

    Qureg qureg = createQureg(3);
    PauliStr str = getInlinePauliStr("XYZ", {0,1,2});
    qcomp angle = getQcomp(.4, .8);

    initPlusState(qureg);
    applyNonUnitaryPauliGadget(qureg, str, angle);

    qreal norm = calcTotalProb(qureg);
    reportScalar("norm", norm);

    finalizeQuESTEnv();
    return 0;
}