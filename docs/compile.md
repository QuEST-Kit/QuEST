# 🛠️  Compile

<!--
  Instructions for compiling QuEST with CMake
  (this comment must be under the title for valid doxygen rendering)

  @author Tyson Jones

  Some notes about this guide:
  - we will always use a build directory called 'build'
  - we will use spaces around cmake argnames and values for clarity, e.g.
    cmake -B build -D ENABLE_CUDA=ON
  - we will demonstrate the simplest and visually clear (and likely sub-optimal) 
    use-cases before progressively more visually complicated examples
-->

QuEST can be compiled with [CMake](https://cmake.org/) to make a standalone executable, or an exported library, or a library installed on the system. 
Compiling is configured with variables supplied by the [`-D` flag](https://cmake.org/cmake/help/latest/command/add_definitions.html) to the [CMake CLI](https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#command-line-cmake-tool). This page details _how_ to compile QuEST for varying purposes and hardwares.


<!-- 
    we are using explicit <a>, rather than markdown links,
    for Doxygen compatibility. It cannot handle [](#sec)
    links, and its <a> anchors are not scoped to files, so
    we here prefix each name with the filename. Grr!
-->

> **TOC**:
> - <a href="#compile_basic">Basic</a>
> - <a href="#compile_optimising">Optimising</a>
> - <a href="#compile_linking">Linking</a>
> - <a href="#compile_configuring">Configuring</a>
>    * <a href="#compile_precision">Precision</a>
>    * <a href="#compile_compilers">Compilers</a>
>    * <a href="#compile_flags">Flags</a>
> - <a href="#compile_examples">Examples</a>
> - <a href="#compile_tests">Tests</a>
>    * <a href="#compile_v4">v4</a>
>    * <a href="#compile_v3">v3</a>
> - <a href="#compile_multithreading">Multithreading</a>
> - <a href="#compile_gpu-acceleration">GPU-acceleration</a>
>    * <a href="#compile_nvidia">NVIDIA</a>
>    * <a href="#compile_amd">AMD</a>
> - <a href="#compile_cuquantum">cuQuantum</a>
> - <a href="#compile_distribution">Distribution</a>
> - <a href="#compile_multi-gpu">Multi-GPU</a>

> **See also**:
> - [`cmake.md`](cmake.md) for the full list of passable compiler variables.
> - [`compilers.md`](compilers.md) for a list of compatible and necessary compilers.
> - [`qtechtheory.org`](https://quest.qtechtheory.org/download/) for help downloading the necessary compilers.
> - [`launch.md`](launch.md) for a guide to executing the compiled application.

> [!TIP]
> QuEST's [Github Actions](https://github.com/QuEST-Kit/QuEST/actions/workflows/compile.yml) regularly test QuEST compilation using a broad combination of deployment settings; presently `108` combinations! The [`compile.yml`](/.github/workflows/compile.yml) workflow can serve as a concrete example of how to compile QuEST in a sanitised, virtual setting.


------------------


<!-- permit doxygen to reference section -->
<a id="compile_basic"></a>

## Basic

Compilation is a two-step process which can generate lots of temporary files and so should be performed in a `build/` folder to avoid clutter. From the `QuEST/` root, run (in terminal):
```bash
# configure
cmake -B build

# build
cmake --build build
```
or more safely from within the `QuEST/build/` folder (as we from here assume):
```bash
# configure
cmake ..

# build
cmake --build .
```

> [!TIP]
> Speed up building by passing [`-j`](https://cmake.org/cmake/help/latest/manual/cmake.1.html#cmdoption-cmake-build-j) or [`--parallel`](https://cmake.org/cmake/help/latest/manual/cmake.1.html#cmdoption-cmake-build-j)
> ```bash
> cmake --build . --parallel
> ```

With no additional arguments, these commands compile [`min_example.cpp`](/examples/tutorials/min_example.cpp) into an executable `min_example` in the `build` folder which can be run via
```bash
./min_example
```

You should expect to see
```
QuEST execution environment:
  [precision]
    qreal.................double (8 bytes)
    qcomp.................std::__1::complex<double> (16 bytes)
    qindex................long long int (8 bytes)
    ...
```


How _boring_! We must pass additional arguments in order to link QuEST to our own code; build other examples; run the unit tests; enable compiler optimisations; enable hardware acceleration; and integrate additional libraries and backends.


------------------



<!-- permit doxygen to reference section -->
<a id="compile_optimising"></a>

## Optimising

QuEST's source code is careful to enable a myriad of optimisations such as [inlining](https://en.wikipedia.org/wiki/Inline_expansion), [loop unrolling](https://en.wikipedia.org/wiki/Loop_unrolling), [auto-vectorisation](https://en.wikipedia.org/wiki/Automatic_vectorization) and [cache optimisations](https://en.wikipedia.org/wiki/Cache_replacement_policies). To utilise them fully, we must instruct our compilers to enable them; like we might do with the [`-O3`](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) flag when invoking a compiler like `gcc` directly.

On most platforms (with the exception of Windows), this is automatic with the commands above, but can otherwise be forced by specifying [`CMAKE_BUILD_TYPE`](https://cmake.org/cmake/help/latest/guide/tutorial/Packaging%20Debug%20and%20Release.html) at _configure time_:

```bash
# configure
cmake .. -D CMAKE_BUILD_TYPE=Release
```

When compiling on **Windows** however (using [Visual Studio](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#visual-studio-generators)), or otherwise using a "[_multi-config generator_](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#other-generators)", we must always supply the build type at _build time_ via [`config`](https://cmake.org/cmake/help/latest/manual/cmake.1.html#cmdoption-cmake-build-config):
```bash
# build
cmake --build . --config Release
```

Otherwise, such generators may default to the `Debug` configuration which can produce executables over `10x` slower than `Release`!

It is always safe to specify either flag even when not used, so one can gaurantee optimisations are enabled by cautiously performing:
```bash
# configure
cmake .. -D CMAKE_BUILD_TYPE=Release

# build
cmake --build . --config Release
```

Read more about CMake generator configurations [here](https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#build-configurations).

> [!TIP]
> Re-configuring a previously compiled CMake project will preserve any manually set variables so they do not need to be re-specified.
> Ergo subsequently running
> ```bash
> # configure
> cmake ..
> ```
> will still incorporate `-D CMAKE_BUILD_TYPE=Release` as previously set.

> [!WARNING]
> The above tip does _not_ apply to _re-building_, for which the `--config Release` _must_ be re-specified (on Windows)



------------------


<!-- permit doxygen to reference section -->
<a id="compile_linking"></a>

## Linking

QuEST can be pre-compiled and later linked to other binaries, _or_ compiled directly alongside the user's source code. 
We focus on the latter use-case, common among scientists when writing simulation scripts. Users seeking to integrate QuEST into larger stacks are likely already familiar with linking libraries through CMake and should check out [`cmake.md`](cmake.md) directly.

To compile a `C` or `C++` file such as
```C
/* myfile.c */

#include "quest.h"

int main() {
    initQuESTEnv();
    finalizeQuESTEnv(); // changed my mind
    return 0;
}
```
simply specify variables `USER_SOURCE` and `OUTPUT_EXE` at _configure time_:
```bash
# configure
cmake .. -D USER_SOURCE=myfile.c -D OUTPUT_EXE=myexec
```
where 
- `myfile.c` is your `C` source file (or `myfile.cpp` if using `C++`).
- `myexec` is the output executable name, which will be saved in `build`.

To compile multiple dependent files, such as
```cpp
/* myfile.cpp */

#include "quest.h"

extern void myfunc();

int main() {
    initQuESTEnv();
    myfunc();
    finalizeQuESTEnv();
    return 0;
}
```
```cpp
/* otherfile.cpp */

#include <stdio.h>

void myfunc() {
    printf("hello quworld!\n");
}
```
simply separate them by `;` in `USER_SOURCE`, wrapped in quotations:
```bash
# configure
cmake .. -D USER_SOURCE="myfile.cpp;otherfile.cpp" -D OUTPUT_EXE=myexec
```


Building then proceeds as normal, e.g.
```bash
# build
cmake --build . --parallel --config Release
```
and the executable can thereafter be run (from within `build`) via
```bash
./myexec
```

You can pass compiler and linker flags needed by your source files through the [`CMAKE_C_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html), [`CMAKE_CXX_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html) and [`CMAKE_EXE_LINKER_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_EXE_LINKER_FLAGS.html) CMake flags as detailed in the <a href="#compile_flags">below section</a>. Note however that if your configuration becomes complicated or your source code requires different `C`/`C++` standards than the QuEST source, you should consider separately compiling QuEST then linking it
to your project as a library!

------------------


<!-- permit doxygen to reference section -->
<a id="compile_configuring"></a>

## Configuring


<!-- permit doxygen to reference section -->
<a id="compile_precision"></a>

### Precision

QuEST's numerical precision can be configured at compile-time, informing what _type_, and ergo how many _bytes_, are used to represent each `qreal` (a floating-point real number) and `qcomp` (a complex amplitude). This affects the memory used by each `Qureg`, but also the user-facing `qreal` and `qcomp` types, as detailed below. Reducing the precision accelerates QuEST at the cost of worsened numerical accuracy. 

Precision is set at configure-time using the `FLOAT_PRECISION` [cmake variable](cmake.md), taking on the values `1`, `2` (default) or `4`.
For example
```bash
# configure
cmake .. -D FLOAT_PRECISION=1
```

The values inform types:

| Value | Precision | `qreal` | size      | `C` (`gcc`) `qcomp` | `C` (`msvc`) `qcomp` | `C++` `qcomp` | size |
|-------|-----------|---------|-----------|---------------------|----------------------|---------------|------|
| `1`   | Single    | `float` | `4` bytes | `float _Complex`    | `_Fcomplex`  | `std::complex<float>` | `8` bytes |
| `2`   | Double    | `double` | `8` bytes | `double _Complex` | `_Dcomplex` | `std::complex<double>` | `16` bytes |
| `4`   | Quadruple*      | `long double` | `<= 16` bytes | `long double _Complex` | `_Lcomplex` | `std::complex<long double>` | `<= 32` bytes |

> [!WARNING]
> While the size of `float` and `double` are fixed by [IEEE 754 format](https://www.intel.com/content/www/us/en/docs/programmable/683242/current/ieee-754-format.html), the size of the `long double` is platform and compiler dependent, and not necessarily a genuine [quadruple-precision float](https://en.wikipedia.org/wiki/Quadruple-precision_floating-point_format). For example, the size of `long double` in Clang can be set to `64`, `80` or `128` bits at [compile-time](https://clang.llvm.org/docs/ClangCommandLineReference.html#long-double-options). Never hardcode the size; always use [`sizeof`](https://en.wikipedia.org/wiki/Sizeof)!



> [!NOTE]
> When enabling <a href="#compile_gpu-acceleration">GPU-acceleration</a>, the precision _must_ be set to `1` or `2` since GPUs do not support quad precision.

<!-- permit doxygen to reference section -->
<a id="compile_compilers"></a>

### Compilers

If multiple compilers are installed, you can choose which to use to compile your `C` and `C++` sources (the latter including the QuEST source) with respective configure-time commands:
```bash
# configure
cmake .. -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++
```
replacing `gcc` and `g++` with e.g. [`clang`](https://clang.llvm.org/), [`cl`](https://learn.microsoft.com/en-us/cpp/build/reference/compiler-options?view=msvc-170), [`icc`](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-8/compiler-commands.html), [`ibm-clang`](https://www.ibm.com/docs/en/open-xl-c-cpp-zos/1.1.0?topic=compiler-command-line-syntax), or aliases for specific versions like `gcc-8.5`.

These compilers will also be used as the _host compilers_ (around which bespoke compilers _wrap_) when enabling GPU-acceleration or distribution.

> [!IMPORTANT]
> It is _not_ correct to specify GPU and MPI compilers, like `nvcc` or `mpicc`, via the above flags. See the respective <a href="#compile_gpu-acceleration">GPU</a> and <a href="#compile_distribution">MPI</a> sections.




<!-- permit doxygen to reference section -->
<a id="compile_flags"></a>

### Flags



Additional flags needed by your files can be passed to the `C` and `C++` compilers, and the linker (respectively), at _configuration time_ via
- [`CMAKE_C_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)
- [`CMAKE_CXX_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)
- [`CMAKE_EXE_LINKER_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_EXE_LINKER_FLAGS.html)

For example,
```bash
# configure
cmake .. -D CMAKE_C_FLAGS="-D MYMACRO=5" -D CMAKE_EXE_LINKER_FLAGS="-lm"
```

Such flags are listed in [`cmake.md`](cmake.md).
However, if your configuration is any more complicated or your source code requires different `C`/`C++` 
standards than the QuEST source, you should consider separately compiling QuEST then linking it
to your source code as a library!


QuEST itself accepts a variety of its preprocessors (mostly related to testing) to be overriden by compiler flags, passed through custom CMake variables, as detailed in [`cmake.md`](cmake.md).



------------------


<!-- permit doxygen to reference section -->
<a id="compile_examples"></a>

## Examples

To compile all of QuEST's [`examples/`](/examples/), use
```bash
# configure
cmake .. -D BUILD_EXAMPLES=ON

# build
cmake --build .
```
The executables will be saved in the (current) `build` directory, in a sub-directory structure mimicking the [`examples/`](/examples/) folder. They can be run by e.g.
```bash
./examples/matrices/cpp_initialisation
```
as elaborated upon in [`launch.md`](launch.md#tests).
<!-- @todo the above link fails in Doxygen; it's too stupid to recognise the section ref -->



------------------


<!-- permit doxygen to reference section -->
<a id="compile_tests"></a>

## Tests


<!-- permit doxygen to reference section -->
<a id="compile_v4"></a>

### v4

To compile QuEST's latest unit and integration tests, use

```bash
# configure
cmake .. -D ENABLE_TESTING=ON

# build
cmake --build .
```
This will compile an executable `tests` in subdirectory `build/tests/`, which can be run as explained in [`launch.md`](launch.md#tests).
<!-- @todo the above link fails in Doxygen; it's too stupid to recognise the section ref -->



<!-- permit doxygen to reference section -->
<a id="compile_v3"></a>

### v3

QuEST's deprecated v3 API has its own unit tests which can be additionally compiled (_except_ on Windows) via
```bash
# configure
cmake .. -D ENABLE_TESTING=ON -D ENABLE_DEPRECATED_API=ON

# build
cmake --build .
```
and run as explained in [`launch.md`](launch.md#v3).
<!-- @todo the above link fails in Doxygen; it's too stupid to recognise the section ref -->



------------------


<!-- permit doxygen to reference section -->
<a id="compile_multithreading"></a>

## Multithreading

Multithreading allows multiple cores of a CPU, or even multiple connected CPUs, to cooperatively perform and ergo accelerate QuEST's expensive functions. Practically all modern computers have the capacity for, and benefit from, multithreading. Note it requires that the CPUs have shared memory (such as through [NUMA](https://learn.microsoft.com/en-us/windows/win32/procthread/numa-support)) and so ergo live in the same machine. CPUs on _different_ machines, connected via a network, can be parallelised over using <a href="#compile_distribution">distribution</a>.

QuEST uses [OpenMP](https://www.openmp.org/) to perform multithreading, so accelerating QuEST over multiple CPUs or cores requires a compiler integrated with OpenMP. This is true of almost all major compilers - see a list of tested compilers in [`compilers.md`](compilers.md#cpu).
<!-- @todo the above link fails in Doxygen; it's too stupid to recognise the section ref -->


> [!IMPORTANT]  
> Using [`Clang`](https://clang.llvm.org/) on MacOS requires use of the `libomp` library, obtainable via [Homebrew](https://brew.sh/):
> ```bash
> brew install libomp
> ```
> and which must be exposed to Clang prior to compilation via
> ```
> export OpenMP_ROOT=$(brew --prefix)/opt/libomp
> ```

To compile with multithreading, simply enable it during configuration:
```bash
# configure
cmake .. -D ENABLE_MULTITHREADING=ON

# build
cmake --build .
```
This is in fact the default behaviour!

The number of threads over which to parallelise QuEST's execution is chosen through setting environment variables, like [`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html), immediately before execution. See [`launch.md`](launch.md#multithreading) for a general guide on multithreaded deployment.
<!-- @todo the above link fails in Doxygen; it's too stupid to recognise the section ref -->



------------------


<!-- permit doxygen to reference section -->
<a id="compile_gpu-acceleration"></a>

## GPU-acceleration

QuEST's core functions perform simple mathematical transformations on very large arrays, and are ergo well suited to parallelisation using general purpose GPUs. This involves creating persistent memory in the GPU VRAM which mirrors the ordinary CPU memory in RAM, and dispatching the transformations to the GPU, updating the GPU memory. The greater number of cores and massive internal memory bandwidth of the GPU can make this extraordinarily faster than using multithreading. 

QuEST supports parallelisation using both NVIDIA GPUs (using CUDA) and AMD GPUs (using HIP). Using either requires obtaining a specialised compiler and passing some GPU-specific compiler flags.



<!-- permit doxygen to reference section -->
<a id="compile_nvidia"></a>

### NVIDIA

> TODO!
> - CUDA-compatible GPGPU
> - nvcc compiler
> - nvidia-smi
> - minimum compute-capability

Check the CUDA compiler is installed correctly via
```bash
nvcc --version
```

To compile your QuEST application with CUDA-acceleration, specify both
```bash
# configure
cmake .. -D ENABLE_CUDA=ON -D CMAKE_CUDA_ARCHITECTURES=$CC
```
where `$CC` is your GPU's [compute capability](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#compute-capabilities) (excluding the full-stop) which you can look up [here](https://developer.nvidia.com/cuda-gpus). 
For example, compiling for the [NVIDIA A100](https://www.nvidia.com/en-us/data-center/a100/) looks like:
```bash
# configure
cmake .. -D ENABLE_CUDA=ON -D CMAKE_CUDA_ARCHITECTURES=80
```


<!-- the below link fails in Doxygen - it's too stupid to recognise the section ref -->
> [!CAUTION]
> Setting the wrong compute capability will cause silently erroneous results. Always run the [unit tests](launch.md#tests) after compiling for the first time to confirm it was set correctly.


Building then proceeds as normal, e.g.
```bash
# build
cmake --build . --parallel
```

<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
The compiled executable can be run like any other, though the GPU behaviour can be prior configured with environment variables. See [`launch.md`](launch.md#gpu-acceleration) for a general guide on GPU-accelerated deployment.




<!-- permit doxygen to reference section -->
<a id="compile_amd"></a>

### AMD

> TODO!
> - ROCm
> - ENABLE_HIP
> - CMAKE_HIP_ARCHITECTURES


To compile your QuEST application with HIP-acceleration, specify both
```bash
# configure
cmake .. -D ENABLE_HIP=ON -D CMAKE_HIP_ARCHITECTURES=$TN
```
where `$TN` is your AMD GPU's [LLVM target name](https://rocm.docs.amd.com/en/latest/reference/gpu-arch-specs.html#glossary). You can look this up [here](https://rocm.docs.amd.com/en/latest/reference/gpu-arch-specs.html), or find the names of all of your local GPUs by running the [ROCM agent enumerator](https://rocm.docs.amd.com/projects/rocminfo/en/latest/how-to/use-rocm-agent-enumerator.html) command, i.e.
```bash
rocm_agent_enumerator -name
```
For example, compiling for the [AMD Instinct MI210 accelerator](https://www.amd.com/en/products/accelerators/instinct/mi200/mi210.html) looks like:
```bash
# configure
cmake .. -D ENABLE_HIP=ON -D CMAKE_HIP_ARCHITECTURES=gfx90a
```


<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
> [!CAUTION]
> Setting the wrong LLVM target name can cause silently erroneous results. Always run the [unit tests](launch.md#tests) after compiling for the first time to confirm it was set correctly.


Building then proceeds as normal, e.g.
```bash
# build
cmake --build . --parallel
```

<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
The compiled executable can be run like any other, though the GPU behaviour can be prior configured with environment variables. See [`launch.md`](launch.md#gpu-acceleration) for a general guide on GPU-accelerated deployment.



------------------

<!-- permit doxygen to reference section -->
<a id="compile_cuquantum"></a>

## cuQuantum

When compiling for NVIDIA GPUs, you can choose to optionally enable [_cuQuantum_](https://docs.nvidia.com/cuda/cuquantum/latest/index.html). This will replace some of QuEST's custom GPU functions with [_cuStateVec_](https://docs.nvidia.com/cuda/cuquantum/latest/custatevec/index.html) routines which are likely to use tailored optimisations for your particular GPU and ergo run faster.


> [!IMPORTANT]
> cuStateVec is presently only supported on Linux, with CUDA `11` or `12`, and modern GPUs with a compute capability equal or above `7.0` (`Volta`, `Turing`, `Ampere`, `Ada`, `Hopper`, `Blackwell`). Check the updated requirements [here](https://docs.nvidia.com/cuda/cuquantum/latest/custatevec/index.html).

Using the cuQuantum backend requires separately downloading cuStateVec, as detailed [here](https://docs.nvidia.com/cuda/cuquantum/latest/getting-started/index.html). Note it is _not_ necessary to download _cuTensorNet_ and _cuDensityMatrix_ which are not used. We recommend downloading...
- on a personal computer, via the [NVIDIA developer site](https://developer.nvidia.com/cuQuantum-downloads).
- on a remote machine, via `wget` as detailed [here](https://docs.nvidia.com/cuda/cuquantum/latest/getting-started/index.html#using-archive).


After download and installation, and before compiling, you must set the `CUQUANTUM_ROOT` environment variable to the cuQuantum download location:
```bash
export CUQUANTUM_ROOT=/path/to/cuquantum-folder
```

Compilation is then simple; we specify `ENABLE_CUQUANTUM` in addition to the above GPU CMake variables. 
For example
```bash
# configure
cmake .. -D ENABLE_CUDA=ON -D CMAKE_CUDA_ARCHITECTURES=80 -D ENABLE_CUQUANTUM=ON

# build
cmake --build . --parallel
```

<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
No other changes are necessary, nor does cuQuantum affect <a href="#compile_multi-gpu">hybridising</a> GPU acceleration and distribution. Launching the executable is the same as in the above section. See [`launch.md`](launch.md#gpu-acceleration).




------------------


<!-- permit doxygen to reference section -->
<a id="compile_distribution"></a>

## Distribution

Because statevectors grow exponentially with the number of simulated qubits, it is easy to run out of memory. In such settings, we may seek to use _distribution_ whereby multiple cooperating machines on a network each store a tractable partition of the state. Distribution can also be useful to speed up our simulations, when the benefit of additional parallelisation outweighs the inter-machine communication penalties.


<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
Enabling distribution requires compiling QuEST with an MPI compiler, such as those listed in [`compilers.md`](compilers.md#comm). Test your compiler is working via
```bash
mpicxx --version
```
> This command may differ on Intel MPI compilers; try `mpiicpc`.


Compiling QuEST's distributed backend is as simple as

```bash
# configure
cmake .. -D ENABLE_DISTRIBUTION=ON

# build
cmake --build . --parallel
```

<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
Note that distributed executables are launched in a distinct way to the other deployment mods, as explained in [`launch.md`](launch.md#distribution),



------------------


<!-- permit doxygen to reference section -->
<a id="compile_multi-gpu"></a>

## Multi-GPU

> TODO!
> - CUDA-aware MPI
> - UCX
> - launch flags
> - checking via reportenv
