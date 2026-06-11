/** @file
 * Quick benchmark for BMI2-assisted bit-index helpers.
 *
 * @author tzh476
 */

#include "quest/src/core/bitwise.hpp"

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

static volatile qindex sinkValue = 0;

template <size_t N>
qindex makeMask(const std::array<int, N>& indices, qindex pattern) {
    qindex mask = 0;
    for (size_t i=0; i<N; i++)
        if ((pattern >> i) & 1)
            mask |= QINDEX_ONE << indices[i];
    return mask;
}

template <size_t N>
double benchGet(const std::string& name, const std::array<int, N>& indices, const std::vector<qindex>& inputs, qindex ampMask) {
    constexpr qindex numIterations = 5000000;
    constexpr int numReps = 5;

    size_t inputMask = inputs.size() - 1;
    double best = std::numeric_limits<double>::max();

    for (int r=0; r<numReps; r++) {
        qindex acc = static_cast<qindex>(0x13579BDF);
        auto start = std::chrono::steady_clock::now();

        for (qindex i=0; i<numIterations; i++) {
            qindex n = (inputs[static_cast<size_t>(i) & inputMask] + acc) & ampMask;
            acc ^= getValueOfBits(n, indices.data(), static_cast<int>(N)) + (i & 7);
        }

        auto end = std::chrono::steady_clock::now();
        sinkValue ^= acc;

        double nsPerCall = std::chrono::duration<double, std::nano>(end - start).count() / static_cast<double>(numIterations);
        best = std::min(best, nsPerCall);
    }

    std::cout << std::left << std::setw(30) << name << " " << std::fixed << std::setprecision(3) << best << " ns/call\n";
    return best;
}

template <size_t N>
double benchInsert(const std::string& name, const std::array<int, N>& indices, const std::vector<qindex>& inputs, qindex valueMask, qindex insertedMask) {
    constexpr qindex numIterations = 5000000;
    constexpr int numReps = 5;

    size_t inputMask = inputs.size() - 1;
    double best = std::numeric_limits<double>::max();

    for (int r=0; r<numReps; r++) {
        qindex acc = static_cast<qindex>(0x2468ACE0);
        auto start = std::chrono::steady_clock::now();

        for (qindex i=0; i<numIterations; i++) {
            qindex n = (inputs[static_cast<size_t>(i) & inputMask] + acc) & valueMask;
            acc ^= insertBitsWithMaskedValues(n, indices.data(), static_cast<int>(N), insertedMask) + (i & 15);
        }

        auto end = std::chrono::steady_clock::now();
        sinkValue ^= acc;

        double nsPerCall = std::chrono::duration<double, std::nano>(end - start).count() / static_cast<double>(numIterations);
        best = std::min(best, nsPerCall);
    }

    std::cout << std::left << std::setw(30) << name << " " << std::fixed << std::setprecision(3) << best << " ns/call\n";
    return best;
}

int main() {
#if defined(QUEST_USE_BMI2_INTRINSICS)
    std::cout << "BMI2 intrinsics: enabled\n";
#else
    std::cout << "BMI2 intrinsics: disabled\n";
#endif

    std::vector<qindex> inputs(1 << 15);
    qindex state = static_cast<qindex>(0x123456789ABCDEFULL);
    for (qindex& input : inputs) {
        state = state * static_cast<qindex>(0x5851F42D4C957F2DULL) + static_cast<qindex>(0x14057B7EF767814FULL);
        input = state;
    }

    qindex nineQubitMask = (QINDEX_ONE << 9) - QINDEX_ONE;
    const std::array<int, 2> inds2 = {2, 7};
    const std::array<int, 5> inds5 = {0, 2, 4, 6, 8};
    const std::array<int, 6> inds6 = {0, 1, 3, 5, 7, 8};

    benchGet("getValueOfBits 2 bits", inds2, inputs, nineQubitMask);
    benchGet("getValueOfBits 5 bits", inds5, inputs, nineQubitMask);
    benchGet("getValueOfBits 6 bits", inds6, inputs, nineQubitMask);

    benchInsert("insertBitsWithMask 2 bits", inds2, inputs, (QINDEX_ONE << 7) - QINDEX_ONE, makeMask(inds2, 0b01));
    benchInsert("insertBitsWithMask 5 bits", inds5, inputs, (QINDEX_ONE << 4) - QINDEX_ONE, makeMask(inds5, 0b10101));
    benchInsert("insertBitsWithMask 6 bits", inds6, inputs, (QINDEX_ONE << 3) - QINDEX_ONE, makeMask(inds6, 0b101011));

    std::cout << "sink: " << sinkValue << "\n";
    return 0;
}
