/*
 * Cost/benefit benchmark for compensated (Kahan) vs naive accumulation in the
 * dense-matrix inner-product kernel of cpu_statevec_anyCtrlAnyTargDenseMatr_sub().
 *
 * This harness replicates the accumulation performed by that QuEST function:
 * each output amplitude is  sum_j  matr[k][j] * cache[j], where j ranges over the
 * 2^numTargets matrix columns. We reproduce that serial inner loop once with naive
 * summation and once with complex Kahan summation. NOTE: the harness uses
 * std::complex<Real> rather than QuEST's base_qcomp type; the two are arithmetically
 * identical here (both perform plain componentwise IEEE real/imaginary arithmetic, so
 * complex Kahan degenerates to two independent real Kahan sums), but this is not
 * literally QuEST's struct. The real-library accuracy claims (regression test, e2e
 * driver) exercise the actual base_qcomp type. Precisions:
 *
 *     fp1  ->  std::complex<float>       (24-bit mantissa)
 *     fp2  ->  std::complex<double>      (53-bit mantissa)
 *     fp4  ->  std::complex<long double> (QuEST's "quad" precision)
 *
 * IMPORTANT HONESTY NOTE on fp4:
 *   On this Apple-Silicon macOS host, BOTH AppleClang and Homebrew GCC report
 *   sizeof(long double)==8 with a 53-bit mantissa -- i.e. long double IS double.
 *   Therefore QuEST's fp4 build is NOT a genuine quadruple-precision build here,
 *   and fp4 accuracy is expected to equal fp2. We report fp4 anyway (it is what a
 *   user compiling QUEST_FLOAT_PRECISION=4 on this machine actually gets) but flag
 *   it as not-truly-quad.
 *
 *   For the ACCURACY GROUND-TRUTH we therefore do NOT use long double. We use GCC's
 *   __float128 (113-bit mantissa, genuine quad) via quadmath as the reference, and
 *   measure each precision's error against that genuinely-higher-precision result.
 *
 * Build (GCC required for __float128 / quadmath):
 *   g++-15 -O3 -std=c++17 bench_kahan.cpp -o bench_kahan -lquadmath
 */

#include <complex>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <quadmath.h>

using clk = std::chrono::steady_clock;

// ---- genuine-quad complex (reference) built on __float128 -------------------
struct qcomp128 {
    __float128 re, im;
};
static inline qcomp128 operator*(qcomp128 a, qcomp128 b) {
    return { a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re };
}
static inline qcomp128 operator+(qcomp128 a, qcomp128 b) {
    return { a.re+b.re, a.im+b.im };
}

// naive accumulation, templated on the working complex type T (== std::complex<...>)
template <typename T, typename Elem, typename Cache>
T accumulate_naive(const Elem& row, const Cache& cache, long dim) {
    T sum(0, 0);
    for (long j = 0; j < dim; j++)
        sum += row[j] * cache[j];
    return sum;
}

// compensated (Kahan) accumulation -- mirrors the QuEST function body exactly
template <typename T, typename Elem, typename Cache>
T accumulate_kahan(const Elem& row, const Cache& cache, long dim) {
    T sum(0, 0);
    T compensation(0, 0);
    for (long j = 0; j < dim; j++) {
        T product     = row[j] * cache[j];
        T corrected   = product - compensation;
        T next        = sum + corrected;
        compensation  = (next - sum) - corrected;
        sum           = next;
    }
    return sum;
}

// reference accumulation in genuine quad (113-bit) -- naive is fine since the
// reference type so vastly out-precisions the candidates that its own rounding
// is negligible relative to fp1/fp2 error.
qcomp128 accumulate_ref(const std::vector<qcomp128>& row, const std::vector<qcomp128>& cache, long dim) {
    qcomp128 sum{0, 0};
    for (long j = 0; j < dim; j++)
        sum = sum + (row[j] * cache[j]);
    return sum;
}

template <typename Real>
struct Result {
    double abs_err_naive;   // |naive   - ref|
    double abs_err_kahan;   // |kahan   - ref|
    double ms_naive;        // ms per full apply (all output amps)
    double ms_kahan;
};

// Run one (precision, numTargets) cell. We build a random unitary-ish dense matrix
// and a random input vector, in genuine quad, then DOWN-CAST to the working precision
// so that naive/kahan see identically-rounded inputs and differ only in the sum.
template <typename Real>
Result<Real> run_cell(int numTargets, int reps, unsigned seed) {
    using T = std::complex<Real>;
    const long dim = 1L << numTargets;          // matrix is dim x dim, vector length dim

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    // reference-precision matrix (row-major) and input cache
    std::vector<qcomp128> matr128((size_t)dim * dim);
    std::vector<qcomp128> cache128(dim);

    // To exercise cancellation realistically, scale entries across many magnitudes:
    // mix O(1) terms with a few very large +/- pairs that must cancel.
    for (long i = 0; i < dim; i++) {
        double mag = std::pow(10.0, uni(rng) * 6.0);   // 1e-6 .. 1e6
        cache128[i] = { (__float128)(uni(rng) * mag), (__float128)(uni(rng) * mag) };
    }
    for (long r = 0; r < dim; r++)
        for (long c = 0; c < dim; c++) {
            double mag = std::pow(10.0, uni(rng) * 6.0);
            matr128[(size_t)r*dim + c] = { (__float128)(uni(rng) * mag), (__float128)(uni(rng) * mag) };
        }

    // down-cast to working precision
    std::vector<T> matr((size_t)dim * dim);
    std::vector<T> cache(dim);
    for (size_t i = 0; i < matr128.size(); i++)
        matr[i] = T((Real)(double)matr128[i].re, (Real)(double)matr128[i].im);
    for (long i = 0; i < dim; i++)
        cache[i] = T((Real)(double)cache128[i].re, (Real)(double)cache128[i].im);

    // genuine-quad reference for every output amp
    std::vector<qcomp128> ref(dim);
    {
        std::vector<qcomp128> rowbuf(dim);
        for (long k = 0; k < dim; k++) {
            for (long j = 0; j < dim; j++)
                rowbuf[j] = matr128[(size_t)k*dim + j];
            ref[k] = accumulate_ref(rowbuf, cache128, dim);
        }
    }

    // helper to view row k of the working-precision matrix
    auto rowptr = [&](long k) { return &matr[(size_t)k*dim]; };

    // accuracy: worst-case (max over output amps) absolute error vs quad reference
    double err_naive = 0, err_kahan = 0;
    for (long k = 0; k < dim; k++) {
        T sn = accumulate_naive<T>(rowptr(k), cache, dim);
        T sk = accumulate_kahan<T>(rowptr(k), cache, dim);
        double dn = std::hypot((double)((__float128)(Real)sn.real() - ref[k].re),
                               (double)((__float128)(Real)sn.imag() - ref[k].im));
        double dk = std::hypot((double)((__float128)(Real)sk.real() - ref[k].re),
                               (double)((__float128)(Real)sk.imag() - ref[k].im));
        if (dn > err_naive) err_naive = dn;
        if (dk > err_kahan) err_kahan = dk;
    }

    // runtime: full apply = compute all dim output amps; repeat for stable timing
    volatile double sink = 0;
    auto time_variant = [&](bool kahan) {
        // warm up
        for (long k = 0; k < dim; k++) {
            T s = kahan ? accumulate_kahan<T>(rowptr(k), cache, dim)
                        : accumulate_naive<T>(rowptr(k), cache, dim);
            sink += (double)s.real();
        }
        auto t0 = clk::now();
        for (int rep = 0; rep < reps; rep++)
            for (long k = 0; k < dim; k++) {
                T s = kahan ? accumulate_kahan<T>(rowptr(k), cache, dim)
                            : accumulate_naive<T>(rowptr(k), cache, dim);
                sink += (double)s.real();
            }
        auto t1 = clk::now();
        double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        return total_ms / reps;   // ms per full apply
    };

    Result<Real> res;
    res.abs_err_naive = err_naive;
    res.abs_err_kahan = err_kahan;
    res.ms_naive = time_variant(false);
    res.ms_kahan = time_variant(true);
    (void)sink;
    return res;
}

template <typename Real>
void run_precision(const char* label, const int* targs, int ntargs) {
    for (int t = 0; t < ntargs; t++) {
        int nt = targs[t];
        // more reps for tiny sizes to get stable timing; fewer for big sizes
        long dim = 1L << nt;
        int reps = (int)std::max(1L, (long)(2'000'000 / (dim * dim)));
        Result<Real> r = run_cell<Real>(nt, reps, 0xC0FFEEu + nt);
        // CSV: precision,numTargets,dim,reps,err_naive,err_kahan,ms_naive,ms_kahan
        printf("%s,%d,%ld,%d,%.6e,%.6e,%.6e,%.6e\n",
               label, nt, dim, reps,
               r.abs_err_naive, r.abs_err_kahan, r.ms_naive, r.ms_kahan);
        fflush(stdout);
    }
}

int main() {
    // report what long double actually is on this host, for the record
    fprintf(stderr, "host: sizeof(float)=%zu sizeof(double)=%zu sizeof(long double)=%zu "
                    "ld_mant=%d  sizeof(__float128)=%zu f128_mant=%d\n",
            sizeof(float), sizeof(double), sizeof(long double),
            (int)__LDBL_MANT_DIG__, sizeof(__float128), (int)FLT128_MANT_DIG);

    const int targs[] = {2, 4, 6, 8, 10, 12, 13, 14};
    const int ntargs = sizeof(targs)/sizeof(targs[0]);

    printf("precision,numTargets,dim,reps,err_naive,err_kahan,ms_naive,ms_kahan\n");
    run_precision<float>       ("fp1", targs, ntargs);  // single
    run_precision<double>      ("fp2", targs, ntargs);  // double
    run_precision<long double> ("fp4", targs, ntargs);  // QuEST "quad" = double on this host
    return 0;
}
