# 🏗️  Architecture

<!--
  Explanation of QuEST's software architecture
  (this comment must be under the title for valid doxygen rendering)
  
  @author Tyson Jones
-->

All user-visible API signatures are contained in `include/`, divided into semantic submodules (like `calculations.h` and `qureg.h`), but all exposed by `quest.h`. They are all strictly `C` _and_ `C++` compatible, hence their `.h` file extension.

The source code within `src/` is divided between five subdirectories, listed below in order of increasing control flow depth. All code is parsed strictly by `C++`, hence all files have `.cpp` and `.hpp` extensions.
- `api/` 
  > contains definitions of the API, directly callable by users (in `C` _or_ `C++`, so it contains de-mangling guards). Functions are divided between files therein similarly to `include/`, and they call only `core/` functions.
- `core/`
  > contains internal non-simulation functions like user-input validaters, non-accelerated maths, and dispatchers to hardware-accelerated backends.
- `comm/` 
  > contains functions needed for exchanging data between distributed nodes, as invoked by the `core/` layer before hardware-accelerated backends.
- `cpu/` 
  > constitutes the multithreaded CPU backend, containing OpenMP-accelerated subroutines and (potentially) AVX intrinsics.
- `gpu/` 
  > constitutes the GPU backend, containing CUDA-accelerated subroutines, GPU hardware queriers, interfaces to CUDA libraries (like Thrust and cuQuantum), and wrappers for AMD/HIP compatibility. Note that some files therein have suffix `.cpp` in lieu of `.cu`, because they are permittedly parsed by non-CUDA compilers when disabling GPU-acceleration at compile-time.

The control flow from the user interface (`quest.h`) to the hardware-accelerated simulation subroutines is mostly captured by:

- `include/quest.h`
- `api/*`
- `core/validation`
  - `core/memory`
  - `gpu/gpu_config`
  - `comm/comm_config`
- `core/utilities`
- `core/localiser`
  - `comm/comm_routines`
- `core/accelerator`
  - `cpu/cpu_subroutines`
  - `gpu/gpu_subroutines`
    - `gpu/gpu_kernels`
    - `gpu/gpu_thrust`
    - `gpu/gpu_cuquantum`


Every API function first validates the user given inputs via `core/validation.cpp`, which in-turn may consult hardware facilities like available RAM (via `core/memory.cpp`) and VRAM (via `gpu/gpu_config.cpp`), or the distributed configuration (via `comm/comm_config.cpp`). Simulation API functions will then invoke `core/localiser.cpp` which checks whether distributed data exchange is necessary in order for all subsequently needed data to become locally available. If so, it invokes functions of `comm/comm_routines.cpp`, then proceeds. Stage `core/accelerator.cpp` chooses whether to process the (potentially just received) data using the CPU (accelerating with `OpenMP`) or the GPU (accelerating with `CUDA`). GPU-acceleration can involve dispatching to custom kernels (in `gpu/cpu_kernels.cpp`), or Thrust routines (in `gpu/gpu_thrust.hpp`), or alternatively to a cuQuantum routine (in `gpu/gpu_cuquantum.hpp`) if optionally compiled. At any call depth, functions within `core/errors.cpp` will be called to confirm preconditions.