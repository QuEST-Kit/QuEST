# Benchmark driver for issue #749 — small GPU allocations

> **Note for maintainers:** this directory is a throw-away benchmarking aid shared
> at the reviewer's request (PR #783). It is **not** wired into CMake/CI and can be
> deleted before merge — the squash will keep it out of the permanent history.

`benchmark_749.cpp` times the two multi-qubit GPU operations targeted by #749, which
previously copied a qubit-index list to the device (`cudaMalloc` + `cudaMemcpyAsync`
+ `cudaFree`) on *every* call:

- `applyMultiQubitProjector`     → `*_multiQubitProjector_sub`
- `calcProbOfMultiQubitOutcome`  → `*_calcProbOfMultiQubitOutcome_sub`

It forces the single-GPU path (`useGpuAccel=1`, distribution/multithreading off),
calls `syncQuESTEnv()` around each timed region (so we measure *completed* GPU work),
and prints per-call wall time (µs) as CSV.

## Build & run

It compiles against QuEST through the built-in `USER_SOURCE_NAMES` mechanism — no
extra include/link flags needed:

```bash
# from a QuEST checkout (this branch):
cmake -S . -B build_bench \
    -D QUEST_ENABLE_CUDA=ON \
    -D CMAKE_CUDA_ARCHITECTURES=120 \
    -D CMAKE_BUILD_TYPE=Release \
    -D USER_SOURCE_NAMES=benchmarks/benchmark_749.cpp \
    -D USER_OUTPUT_EXE_NAME=bench_749
cmake --build build_bench --target bench_749 -j

# args: [minQubits=4] [maxQubits=20] [numTargs=3] [reps=2000]
./build_bench/bench_749 4 20 3 2000
```

To get the before/after numbers, build the same file once against a clean
`origin/devel` checkout (baseline) and once against this branch (optimised).

For the CUDA-API counts, wrap a shorter run under Nsight Systems and tally the
runtime API calls, e.g.:

```bash
nsys profile --trace=cuda -o trace_749 ./build_bench/bench_749 8 12 3 200
nsys stats --report cuda_api_sum trace_749.nsys-rep   # cudaMalloc/Free/MemcpyAsync/LaunchKernel
```

## Results on the author's machine (RTX PRO 6000 Blackwell, CUDA 13.0, sm_120)

CUDA runtime API call counts (N = 8…12, 200 reps, both ops):

| CUDA runtime call | baseline | optimised | change |
|---|--:|--:|--:|
| `cudaMalloc`      | 3020 | 1010 | −66% |
| `cudaFree`        | 3021 | 1011 | −66% |
| `cudaMemcpyAsync` | 3015 | 1005 | −67% |
| `cudaLaunchKernel`| 2025 | 2025 | unchanged |

Per-call wall time (µs), `numTargs = 3`, 2000 reps:

| N | projector base | projector opt | speedup | prob base | prob opt | speedup |
|--:|--:|--:|--:|--:|--:|--:|
| 4  | 12.41 | 6.48 | 1.92× | 20.35 | 14.59 | 1.39× |
| 8  | 11.74 | 6.40 | 1.83× | 19.44 | 14.18 | 1.37× |
| 12 | 11.85 | 6.52 | 1.82× | 19.79 | 14.44 | 1.37× |
| 16 | 12.22 | 6.84 | 1.79× | 28.63 | 23.21 | 1.23× |
| 20 | 58.41 | 13.74 | 4.25× | 67.99 | 61.99 | 1.10× |

The residual ~1000 `cudaMalloc` in the optimised build is `thrust::reduce`'s own
internal temporary inside `calcProb` (inherent to Thrust, out of scope).
