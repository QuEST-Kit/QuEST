#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "quest/include/quest.h"

int main() {
    initQuESTEnv();

    std::vector<int> n_qubits = {24};
    std::vector<int> k_swaps = {3, 5, 8};

    std::cout << "===================================================================\n";
    std::cout << "               QuEST v4 SWAP Fusion Benchmark                      \n";
    std::cout << "===================================================================\n";
    std::cout << std::left << std::setw(12) << "Qubits (n)"
              << std::setw(12) << "Pairs (k)"
              << std::setw(18) << "Sequential (s)"
              << std::setw(18) << "Fused (s)"
              << "Speedup\n";
    std::cout << "-------------------------------------------------------------------\n";

    for (int n : n_qubits) {
        for (int k : k_swaps) {

            std::vector<int> targs1(k), targs2(k);
            for (int i = 0; i < k; ++i) {
                targs1[i] = 2 * i;
                targs2[i] = 2 * i + 1;
            }

            // --- Sequential: create, run, destroy ---
            Qureg qureg_seq = createQureg(n);
            initZeroState(qureg_seq);
            auto start_seq = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < k; ++i)
                applySwap(qureg_seq, targs1[i], targs2[i]);
            auto end_seq = std::chrono::high_resolution_clock::now();
            destroyQureg(qureg_seq);  // FREE before allocating next

            // --- Fused: create, run, destroy ---
            Qureg qureg_fused = createQureg(n);
            initZeroState(qureg_fused);
            auto start_fused = std::chrono::high_resolution_clock::now();
            applyMultiSwap(qureg_fused, targs1, targs2);
            auto end_fused = std::chrono::high_resolution_clock::now();
            destroyQureg(qureg_fused);  // FREE immediately

            std::chrono::duration<double> diff_seq   = end_seq   - start_seq;
            std::chrono::duration<double> diff_fused = end_fused - start_fused;
            double speedup = diff_seq.count() / diff_fused.count();

            std::cout << std::left << std::setw(12) << n
                      << std::setw(12) << k
                      << std::fixed << std::setprecision(6)
                      << std::setw(18) << diff_seq.count()
                      << std::setw(18) << diff_fused.count()
                      << std::setprecision(2) << speedup << "x\n";
        }
    }

    finalizeQuESTEnv();
    return 0;
}