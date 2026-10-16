#include "quest/include/quest.h"
#include <stdio.h>
#include <vector>
#include <cmath>

int main() {

    initQuESTEnv();

    Qureg qureg = createQureg(4);  // 16 amplitudes, indices 0-15
    initZeroState(qureg);          // |0000>, amp[0] = 1

    // flip qubit 0 → |0001>, amp[1] = 1
    applyPauliX(qureg, 0);

    // SWAP qubit 0 ↔ qubit 3: |0001> → |1000>, amp[8] = 1
    std::vector<int> a = {0};
    std::vector<int> b = {3};
    applyMultiSwap(qureg, a, b);

    // print non-zero amps
    printf("Statevector after SWAP(0,3) on |0001>:\n");
    for (int i = 0; i < 16; i++) {
        qcomp amp = getQuregAmp(qureg, i);
        qreal re = std::real(amp);
        qreal im = std::imag(amp);
        if (re*re + im*im > 0.01)
            printf("  amp[%d] = (%.3f, %.3f)\n", i, re, im);
    }

    printf("EXPECT: amp[8] = (1.000, 0.000)\n");
    qcomp result = getQuregAmp(qureg, 8);
    printf("%s\n", (std::real(result) > 0.99) ? "PASS" : "FAIL");

    destroyQureg(qureg);
    finalizeQuESTEnv();
    return 0;
}