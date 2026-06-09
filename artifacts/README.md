# Kahan vs naive accumulation in `cpu_statevec_anyCtrlAnyTargDenseMatr_sub()`

Benchmark + analysis for issue #598 / PR #784. Measures the runtime and accuracy
effect of compensated (Kahan) summation vs the original uncompensated accumulation
in the CPU dense-matrix inner-product loop, across the three QuEST `qcomp`
precisions, on a single CPU.

## Files

| file | purpose |
| --- | --- |
| `bench_kahan.cpp` | Standalone harness replicating the exact accumulation kernel (naive, Kahan, and a genuine-quad `__float128` reference) for fp1/fp2/fp4. |
| `plot_kahan.py` | Generates the three figures from `results.csv`. |
| `e2e_compmatr.cpp` | End-to-end check against the real built QuEST library (Kahan build vs `QUEST_DENSE_ACCUM_NAIVE` build) using the recipe from issue #598. |
| `results.csv` | Raw benchmark data. |
| `kahan_accuracy.png`, `kahan_runtime.png`, `kahan_costbenefit.png` | Figures. |

## How the function was made selectable

`cpu_statevec_anyCtrlAnyTargDenseMatr_sub()` now wraps the accumulation in a
compile-time switch. The default is compensated (Kahan) summation (preserving the
PR behaviour and its regression test); defining `QUEST_DENSE_ACCUM_NAIVE` reverts to
the original `sum += elem * cache[j]` for apples-to-apples benchmarking. Example:

```bash
cmake -S . -B build -D QUEST_FLOAT_PRECISION=2 -D QUEST_ENABLE_OMP=OFF -D QUEST_ENABLE_MPI=OFF
# naive variant for benchmarking:
cmake -S . -B build_naive ... -D CMAKE_CXX_FLAGS="-DQUEST_DENSE_ACCUM_NAIVE"
```

## Reproduce

```bash
# benchmark (GCC needed for __float128 genuine-quad reference)
g++-15 -O3 -std=c++17 bench_kahan.cpp -o bench_kahan -lquadmath
./bench_kahan > results.csv 2> host_info.txt
python3 plot_kahan.py
```

## The three-precision reality on this host (IMPORTANT)

On this Apple-Silicon macOS host, **both AppleClang and Homebrew GCC report
`sizeof(long double) == 8` with a 53-bit mantissa** — i.e. `long double` is just
`double`. QuEST's fp4 build (which uses `long double`) is therefore **NOT a genuine
quadruple-precision build here**, and its accuracy is identical to fp2. This is
visible in `results.csv`: every fp4 error equals the corresponding fp2 error exactly.

Consequently:
- The fp1/fp2/fp4 **runtime** numbers are all valid (they are what those builds
  actually do on this machine).
- For the **accuracy ground truth** we do NOT use `long double`. We use GCC's
  `__float128` (113-bit mantissa, genuine quad, via quadmath) as the reference and
  measure each precision's error against it.

## Key results (worst-case abs error vs 113-bit `__float128` reference)

Adversarial ill-conditioned matrix (entries spanning 1e-6 .. 1e6). `ms` is per
single-CPU `applyCompMatr`; lower is better.

| precision | targets | err naive | err Kahan | err reduction | ms naive | ms Kahan | slowdown |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: |
| fp1 | 8  | 1.57e5 | 4.90e4 | 3.2x | 0.038 | 0.141 | 3.7x |
| fp1 | 14 | 9.54e6 | 3.21e5 | 30x  | 184   | 611   | 3.3x |
| fp2 | 8  | 3.49e-4 | 5.28e-5 | 6.6x | 0.032 | 0.133 | 4.2x |
| fp2 | 14 | 3.29e-2 | 6.54e-4 | 50x  | 175   | 608   | 3.5x |
| fp4 | 14 | 3.29e-2 | 6.54e-4 | 50x  | 173   | 608   | 3.5x |  *(== fp2; not true quad)*

## Conclusion

- **Accuracy benefit grows with the number of targets**, exactly as the issue
  anticipated: negligible at 2–4 targets (no cancellation yet), then 1–2 orders of
  magnitude error reduction by 12–14 targets, for both single and double precision.
- **Cost is roughly constant at ~3.3–4.2x slowdown** at every non-trivial size. The
  Kahan loop is compute-bound here (4 complex ops per term vs 1), and the inner loop
  is not memory-bandwidth bound at these matrix sizes, so the overhead does not
  amortise away. It is *not* "time-free".
- **End-to-end** the modified function turns a fully cancelled `amp0 = 0` (naive)
  into the correct `~4096` (Kahan) for a 12-target adversarial `CompMatr`.

So Kahan is worthwhile for large or ill-conditioned dense matrices where the
accuracy matters, but a blanket 3–4x slowdown is too costly for small/well-conditioned
matrices. A size threshold (e.g. enable Kahan only above ~6–8 targets) is the natural
compromise; the `QUEST_DENSE_ACCUM_NAIVE` toggle is the simplest hook for that.
