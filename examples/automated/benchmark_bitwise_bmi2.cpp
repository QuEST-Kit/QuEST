/** @file
 * A quick, self-contained micro-benchmark of the BMI2 PEXT/PDEP fast paths added for issue #717,
 * comparing them against the original scalar bit gather/scatter loops. It prints per-call timings
 * so QuEST's CI can compare the speedup across its tested platforms and compilers.
 *
 * The two scalar routines below mirror getValueOfBits() and insertBitsWithMaskedValues() from
 * quest/src/core/bitwise.hpp; the BMI2 routines are the single-instruction _pext_u64 / _pdep_u64
 * paths. This file deliberately depends on nothing but the C++ standard library (and <immintrin.h>
 * when targeting x86 BMI2), so it compiles and runs on every platform — emitting the scalar
 * timings alone where BMI2 is unavailable, never raising SIGILL.
 *
 * Build note: this target is compiled with -mbmi2 (see examples/automated/CMakeLists.txt) so the
 * intrinsic path is enabled; the QuEST library itself enables -mbmi2 the same way in the top-level
 * CMakeLists.txt. Whether the fast path was compiled in is printed at runtime.
 *
 * @author (issue #717 contribution)
 */

#include <cstdint>
#include <cstdio>
#include <chrono>

#if defined(__BMI2__) && (defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86))
  #include <immintrin.h>
  #define BENCH_USE_BMI2
#endif

using std::uint64_t;

// --- scalar references (mirroring quest/src/core/bitwise.hpp) -------------------------------------

// getValueOfBits: gather the bits at the given (strictly increasing) positions into the low bits.
static inline uint64_t scalarGather(uint64_t number, const int* inds, int n) {
    uint64_t value = 0;
    for (int i=0; i<n; i++)
        value |= ((number >> inds[i]) & 1ULL) << i;
    return value;
}

// insertBitsWithMaskedValues: spread number's low bits into the positions NOT named by inds (i.e.
// insert a 0 at each increasing index), then OR in the precomputed value mask.
static inline uint64_t scalarScatter(uint64_t number, const int* inds, int n, uint64_t valueMask) {
    uint64_t r = number;
    for (int i=0; i<n; i++) {
        uint64_t lo = r & ((1ULL << inds[i]) - 1);
        uint64_t hi = r & ~((1ULL << inds[i]) - 1);
        r = (hi << 1) | lo;
    }
    return valueMask | r;
}

static inline uint64_t makePosMask(const int* inds, int n) {
    uint64_t m = 0;
    for (int i=0; i<n; i++)
        m |= 1ULL << inds[i];
    return m;
}

// --- timing harness ------------------------------------------------------------------------------

static double nsPerCall(uint64_t iters, double seconds) {
    return 1e9 * seconds / (double) iters;
}

template <typename F>
static double timeMin(uint64_t iters, int reps, F&& fn) {
    double best = 1e300;
    for (int r=0; r<reps; r++) {
        auto t0 = std::chrono::steady_clock::now();
        fn(iters);
        auto t1 = std::chrono::steady_clock::now();
        double s = std::chrono::duration<double>(t1 - t0).count();
        if (s < best) best = s;
    }
    return best;
}

int main() {

    printf("QuEST issue #717 - BMI2 PEXT/PDEP bitwise micro-benchmark\n");
#ifdef BENCH_USE_BMI2
    printf("BMI2 fast path: ACTIVE (compiled with -mbmi2)\n\n");
#else
    printf("BMI2 fast path: INACTIVE (x86 BMI2 not targeted; scalar timings only)\n\n");
#endif

    const uint64_t iters = 8000000;   // keeps total runtime well under a second
    const int reps = 3;
    const int counts[] = {3, 6};      // representative qubit-arity per gate

    printf("%-8s %-4s %14s %14s %10s\n", "op", "k", "scalar ns/call", "bmi2 ns/call", "speedup");

    for (int ci=0; ci<2; ci++) {
        int k = counts[ci];

        // a fixed, strictly-increasing index set and a value mask consistent with it
        int inds[8];
        for (int i=0; i<k; i++) inds[i] = 3*i + 1;
        uint64_t posMask = makePosMask(inds, k);
        uint64_t valueMask = posMask & 0xA5A5A5A5A5A5A5A5ULL;

        volatile uint64_t sink = 0;

        // ---- gather (getValueOfBits) ----
        double sg = timeMin(iters, reps, [&](uint64_t N){
            uint64_t acc = 0;
            for (uint64_t n=0; n<N; n++) acc ^= scalarGather(n, inds, k);
            sink ^= acc;
        });
#ifdef BENCH_USE_BMI2
        double bg = timeMin(iters, reps, [&](uint64_t N){
            uint64_t acc = 0;
            for (uint64_t n=0; n<N; n++) acc ^= (uint64_t) _pext_u64(n, posMask);
            sink ^= acc;
        });
        printf("%-8s %-4d %14.3f %14.3f %9.2fx\n", "gather", k,
               nsPerCall(iters, sg), nsPerCall(iters, bg), sg/bg);
#else
        printf("%-8s %-4d %14.3f %14s %10s\n", "gather", k, nsPerCall(iters, sg), "-", "-");
#endif

        // ---- scatter (insertBitsWithMaskedValues) ----
        double ss = timeMin(iters, reps, [&](uint64_t N){
            uint64_t acc = 0;
            for (uint64_t n=0; n<N; n++) acc ^= scalarScatter(n, inds, k, valueMask);
            sink ^= acc;
        });
#ifdef BENCH_USE_BMI2
        double bs = timeMin(iters, reps, [&](uint64_t N){
            uint64_t acc = 0;
            for (uint64_t n=0; n<N; n++) acc ^= (valueMask | (uint64_t) _pdep_u64(n, ~posMask));
            sink ^= acc;
        });
        printf("%-8s %-4d %14.3f %14.3f %9.2fx\n", "scatter", k,
               nsPerCall(iters, ss), nsPerCall(iters, bs), ss/bs);
#else
        printf("%-8s %-4d %14.3f %14s %10s\n", "scatter", k, nsPerCall(iters, ss), "-", "-");
#endif

#ifdef BENCH_USE_BMI2
        // sanity: the intrinsic and scalar paths must agree (bit-for-bit) for these sorted indices
        bool ok = true;
        for (uint64_t n=0; n<4096 && ok; n++) {
            if ((uint64_t)_pext_u64(n, posMask) != scalarGather(n, inds, k)) ok = false;
            if ((valueMask | (uint64_t)_pdep_u64(n, ~posMask)) != scalarScatter(n, inds, k, valueMask)) ok = false;
        }
        printf("           (k=%d results verified bit-identical to scalar: %s)\n", k, ok ? "yes" : "NO");
#endif
        (void) sink;
    }

    return 0;
}
