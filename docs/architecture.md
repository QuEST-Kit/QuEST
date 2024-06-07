
All user-visible API signatures are contained in `include/`, divided into semantic submodules (like `calculations.h` and `qureg.h`), but all exposed by `quest.h`. They are all strictly `C` _and_ `C++` compatible, hence their `.h` file extension.

The source code within `src/` is divided between four subdirectories, listed below in order of increasing control flow depth. All code is parsed strictly by `C++`, hence all files have `.cpp` and `.hpp` extensions.
- `api/` contains definitions of the entire API, divided semantically similarly to `include/`, and which call only `core/` functions.
- `core/` contains internal non-simulation functions like user-input validaters, non-accelerated maths, and dispatchers to hardware-accelerated backends.
- `comm/` contains functions needed for exchanging data between distributed nodes, as invoked by the `core/` layer before hardware-accelerated backends.
- `cpu/` constitutes the multithreaded CPU backend, containing OpenMP-accelerated subroutines and AVX intrinsics.
- `gpu/` constitutes the GPU backend, containing CUDA-accelerated subroutines, GPU hardware queriers, interfaces to CUDA libraries (like Thrust and cuQuantum), and wrappers for AMD/HIP compatibility. Note all source files therein end in `.cpp` in lieu of `.cu`, because they are permittedly parsed by non-CUDA compilers when disabling GPU-acceleration at compile-time.

The control flow from the user interface (`quest.h`) to the hardware-accelerated simulation subroutines is mostly captured by:

- `include/quest.h`
- `api/*.cpp`
- `core/validation.cpp`
  - `core/memory.cpp`
  - `gpu/config.cpp`
- `core/localiser.cpp`
  - `comm/communication.cpp`
- `core/accelerator.cpp`
  - `cpu/omp_subroutines.cpp`
  - `gpu/cuda_subroutines.cpp`

Every API function first validates the user given inputs via `validation.cpp`, which in-turn may consult hardware facilities like available RAM (via `core/memory.cpp`) and VRAM (via `gpu/config.cpp`). Simulation API functions will then invoke `localiser.cpp` which checks whether distributed data exchange is necessary in order for all subsequently needed data to become locally available. If so, it invokes functions of `communication.cpp`, then proceeds. Stage `accelerator.cpp` chooses whether to process the (potentially just received) data using the CPU (accelerating with `OpenMP`) or the GPU (accelerating with `CUDA`).