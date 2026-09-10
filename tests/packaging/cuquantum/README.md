# cuQuantum finder checks

Run the synthetic filesystem discovery fixtures without CUDA hardware or an SDK:

```sh
cmake -S tests/packaging/cuquantum -B build/cuquantum-fixtures
ctest --test-dir build/cuquantum-fixtures --output-on-failure
```

Each case uses a fresh configuration and a separate prefix containing spaces.
The SDK headers and library artifacts are fixture files; CUDA imported targets
are supplied by an isolated test module. These checks test discovery and target
interfaces, not binary ABI or device execution.

Run the real SDK compile/link check using a C++ compiler only:

```sh
cmake -S tests/packaging/cuquantum/real-sdk -B build/cuquantum-real \
  -DCUQUANTUM_ROOT=/path/to/cuquantum -DCUDAToolkit_ROOT=/path/to/cuda
cmake --build build/cuquantum-real
ctest --test-dir build/cuquantum-real --output-on-failure
```

The executable calls `custatevecGetVersion()` and checks the component major
version against the header. It needs the SDK shared libraries at runtime but no
GPU operations. Full QuEST CUDA installation and relocation checks live in the
parent packaging test suite.

The finder exports `CUQUANTUM_cuStateVec_VERSION`,
`CUQUANTUM_cuTensorNet_VERSION`, and `CUQUANTUM_cuDensityMat_VERSION` from each
component's own header. It deliberately rejects package-level version requests:
these headers do not supply the overall cuQuantum SDK release number.
