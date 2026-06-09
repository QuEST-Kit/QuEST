/* benchmark_749.cpp
 *
 * Micro-benchmark for QuEST issue #749 ("Profile and optimise-away small GPU
 * allocations"). It times the two families of multi-qubit GPU operations that
 * previously copied a qubit-index list to the device on every call:
 *
 *   - applyMultiQubitProjector()      -> thrust_*_multiQubitProjector_sub
 *   - calcProbOfMultiQubitOutcome()   -> thrust_*_calcProbOfMultiQubitOutcome_sub
 *
 * For a sweep of (small) Qureg sizes it reports the mean wall-clock time per
 * call (in microseconds) as CSV on stdout. Build it once against the baseline
 * (unmodified origin/devel) and once against the optimised branch; the
 * accompanying analyze.py compares the two CSVs and plots the speedup.
 *
 * The GPU is forced on (useGpuAccel=1) and distribution/multithreading are off,
 * so we isolate the single-GPU code path. syncQuESTEnv() is called before each
 * timed region so we measure completed GPU work, not just async dispatch.
 *
 * Usage:  ./bench_749 [minQubits] [maxQubits] [numTargs] [reps]
 *   defaults: minQubits=4 maxQubits=20 numTargs=3 reps=2000
 */

#include "quest.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using clk = std::chrono::high_resolution_clock;

static double timeProjector(Qureg qureg, const std::vector<int>& qubits,
                            const std::vector<int>& outcomes, int reps) {

    // collapse once so subsequent identical projections are well-defined (prob=1)
    initPlusState(qureg);
    applyMultiQubitProjector(qureg, qubits, outcomes);
    syncQuESTEnv();

    auto t0 = clk::now();
    for (int r = 0; r < reps; r++)
        applyMultiQubitProjector(qureg, qubits, outcomes);
    syncQuESTEnv();
    auto t1 = clk::now();

    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
    return ns / reps / 1e3; // microseconds per call
}

static double timeProbOfOutcome(Qureg qureg, const std::vector<int>& qubits,
                                const std::vector<int>& outcomes, int reps) {

    initPlusState(qureg);
    volatile qreal sink = 0;
    sink += calcProbOfMultiQubitOutcome(qureg, qubits, outcomes);
    syncQuESTEnv();

    auto t0 = clk::now();
    for (int r = 0; r < reps; r++)
        sink += calcProbOfMultiQubitOutcome(qureg, qubits, outcomes);
    syncQuESTEnv();
    auto t1 = clk::now();

    double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
    (void) sink;
    return ns / reps / 1e3; // microseconds per call
}

int main(int argc, char** argv) {

    int minQubits = (argc > 1) ? atoi(argv[1]) : 4;
    int maxQubits = (argc > 2) ? atoi(argv[2]) : 20;
    int numTargs  = (argc > 3) ? atoi(argv[3]) : 3;
    int reps      = (argc > 4) ? atoi(argv[4]) : 2000;

    initQuESTEnv();

    // CSV header (printed once; '#' lines are ignored by analyze.py)
    printf("numQubits,numTargs,proj_us,prob_us\n");

    for (int n = minQubits; n <= maxQubits; n++) {

        // statevector forced onto the GPU only
        Qureg qureg = createCustomQureg(n, 0, /*distrib*/0, /*gpu*/1, /*mt*/0);

        std::vector<int> qubits(numTargs), outcomes(numTargs, 0);
        for (int i = 0; i < numTargs; i++)
            qubits[i] = i;

        double proj = timeProjector(qureg, qubits, outcomes, reps);
        double prob = timeProbOfOutcome(qureg, qubits, outcomes, reps);

        printf("%d,%d,%.4f,%.4f\n", n, numTargs, proj, prob);
        fflush(stdout);

        destroyQureg(qureg);
    }

    finalizeQuESTEnv();
    return 0;
}
