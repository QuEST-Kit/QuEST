# 🔧  Compilers

<!--
  A summary of necessary compilers to use QuEST's
  various backend parallelisation deployments
  (this comment must be under the title for valid doxygen rendering)
  
  @author Tyson Jones
-->

QuEST separates compilation of the _frontend_, _backend_ and the _tests_, which have progressively stricter compiler requirements.
This page details the specialised compilers necessary to enable specific features hardware accelerators, and lists such compilers which are
known to be compatible with QuEST.


<!-- 
    we are using explicit <a>, rather than markdown links,
    for Doxygen compatibility. It cannot handle [](#sec)
    links, and its <a> anchors are not scoped to files, so
    we here prefix each name with the filename. Grr!
-->

> **TOC**:
> - <a href="#compilers_frontend">Frontend</a>
> - <a href="#compilers_backend">Backend</a>
>    * <a href="#compilers_comm">Comm</a>
>    * <a href="#compilers_cpu">cpu</a>
>    * <a href="#compilers_gpu">gpu</a>
>    * <a href="#compilers_comm-gpu">comm + gpu</a>
>    * <a href="#compilers_gpu-cuquantum">gpu + cuquantum</a>
> - <a href="#compilers_tests">Tests</a>

> **See also**:
> - [`compile.md`](compile.md) for a guide to compiling QuEST.
> - [`cmake.md`](cmake.md) for the full list of passable compiler variables.
> - [`qtechtheory.org`](https://quest.qtechtheory.org/download/) for help downloading the compilers listed on this page.

> [!TIP]
> QuEST's [Github Actions](https://github.com/QuEST-Kit/QuEST/actions/workflows/compile.yml) regularly test QuEST compilation using a broad combination of compilers; presently `108` combinations! This can provide a clue as to which modern compiler versions are supported, and be a concrete example of _how_ to compile in a sanitised, virtual setting. Check out the [`compile.yml`](/.github/workflows/compile.yml) workflow.



---------------


<!-- permit doxygen to reference section -->
<a id="compilers_frontend"></a>

## Frontend

[![Languages](https://img.shields.io/badge/C-11-ff69b4.svg)](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3631.pdf)
[![Languages](https://img.shields.io/badge/C++-14-ff69b4.svg)](https://isocpp.org/wiki/faq/cpp14)

User code can be written in either `C11` or `C++14`, and has so far been tested with compilers
- [clang](https://clang.llvm.org/) 15
- [gnu](https://gcc.gnu.org/) 10-14
- [msvc](https://learn.microsoft.com/en-us/cpp/build/reference/compiling-a-c-cpp-program?view=msvc-170) 19

> The standards are imposed by the QuEST header's use of `C11` [generics](https://en.cppreference.com/w/c/language/generic) and `C++14` complex arithmetic overloads. Each can be relaxed to enable compatibility with `C99` and `C++11` by simple modifications to the headers - ask us for help! 


---------------


<!-- permit doxygen to reference section -->
<a id="compilers_backend"></a>

## Backend

[![Languages](https://img.shields.io/badge/C++-17-ff69b4.svg)](https://en.cppreference.com/w/cpp/17)

The backend is divided into subdirectories [`api/`](/quest/src/api), [`core/`](/quest/src/core), [`comm/`](/quest/src/comm),  [`cpu/`](/quest/src/cpu) and [`gpu/`](/quest/src/gpu). All can be compiled with a generic `C++17` compiler, but enabling distribution, multithreading and GPU-acceleration requires using specialised compilers for the latter three. Each can be toggled and compiled independently. Note however that tightly-coupled multi-GPU simulations (`comm` + `gpu`) can be accelerated using bespoke compilers, and use of [cuQuantum](https://developer.nvidia.com/cuquantum-sdk) requires modern compilers (`gpu + cuquantum`), detailed below.


<!-- permit doxygen to reference section -->
<a id="compilers_comm"></a>

### comm

Enabling distribution requires compiling `comm/` with an [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface)-compatible compiler, which has so far been tested with
- [openmpi](https://www-lb.open-mpi.org/software/ompi/v4.0/) 4
- [openmpi](https://www.open-mpi.org/software/ompi/v5.0/) 5
- [mpich](https://www.mpich.org/) 4
- [msmpi](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi) 10
- [impi](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html) 3

when wrapping all previously mentioned compilers.


<!-- permit doxygen to reference section -->
<a id="compilers_cpu"></a>

### cpu

Enabling multithreading requires compiling `cpu/` with an [OpenMP](https://www.openmp.org/)-compatible compiler. Versions
- [openmp](https://www.openmp.org/specifications/) 2.0
- [openmp](https://www.openmp.org/specifications/) 4.5

have been explicitly tested, as used by the aforementioned compilers.


<!-- permit doxygen to reference section -->
<a id="compilers_gpu"></a>

### gpu

Enabling acceleration on NVIDIA or AMD GPUs requires compiling `gpu/` with a [CUDA](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/) or [ROCm](https://rocm.docs.amd.com/en/docs-6.0.2/) compiler respectively. These must be compatible with [Thrust](https://developer.nvidia.com/thrust) and [rocThrust](https://github.com/ROCm/rocThrust) respectively. QuEST v4 has been so far tested with
- [cuda](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html) 11
- [cuda](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html) 12


<!-- permit doxygen to reference section -->
<a id="compilers_comm-gpu"></a>

### comm + gpu

Simultaneously emabling both distribution _and_ GPU-acceleration is possible with use of the separate compilers above. However, simulation can be accelerated by using a [CUDA-aware MPI](https://developer.nvidia.com/blog/introduction-cuda-aware-mpi/) compiler, enabling QuEST to use [GPUDirect](https://developer.nvidia.com/gpudirect) and avoid superfluous exchanges of CPU and GPU memories. So far, QuEST has been tested with:
- [UCX](https://openucx.org/) 1.13
- [UCX](https://openucx.org/) 1.15


<!-- permit doxygen to reference section -->
<a id="compilers_gpu-cuquantum"></a>

### gpu + cuquantum

Enabling [cuQuantum](https://developer.nvidia.com/cuquantum-sdk) on NVIDIA GPUs with [compute-capability](https://developer.nvidia.com/cuda-gpus) >= `7.0` requires use of a modern CUDA compiler, specifically
- [cuda](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html) 11
- [cuda](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html) 12


---------------


<!-- permit doxygen to reference section -->
<a id="compilers_tests"></a>

## Tests

[![Languages](https://img.shields.io/badge/C++-20-ff69b4.svg)](https://en.cppreference.com/w/cpp/20)

QuEST's [`tests/`](/tests/) make use of several `C++20` features and may not be compatible with older compilers. So far, they have been tested with
- [clang](https://clang.llvm.org/) 15
- [gnu](https://gcc.gnu.org/) 13
- [gnu](https://gcc.gnu.org/) 14
- [msvc](https://learn.microsoft.com/en-us/cpp/build/reference/compiling-a-c-cpp-program?view=msvc-170) 19
