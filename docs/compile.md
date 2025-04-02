<!--
  Instructions for compiling QuEST with CMake

  @author Tyson Jones

  Some notes about this guide:
  - we will always use a build directory called 'build'
  - we will use spaces around cmake argnames and values for clarity, e.g.
    cmake -B build -D ENABLE_CUDA=ON
  - we will demonstrate the simplest and visually clear (and likely sub-optimal) 
    use-cases before progressively more visually complicated examples
-->




# Compiling

QuEST can be compiled with [CMake](https://cmake.org/) to make a standalone executable, or an exported library, or a library installed on the system. 
Compiling is configured with variables supplied by the [`-D` flag](https://cmake.org/cmake/help/latest/command/add_definitions.html) to the [CMake CLI](https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#command-line-cmake-tool). This page details _how_ to compile QuEST for varying purposes and hardwares.

**TOC**:
- [Basic](#basic)
- [Optimisation](#optimisation)
- [Linking](#linking)
- [Examples](#examples)
- [Tests](#tests)
- [Multithreading](#multithreading)
- [GPU-acceleration](#gpu-acceleration)
- [cuQuantum](#cuquantum)
- [Distribution](#distribution)
- [GPUDirect](#gpudirect)

> **See also**:
> - [`cmake.md`](cmake.md) for the full list of passable compiler variables.
> - [`compilers.md`](compilers.md) for a list of compatible and necessary compilers.
> - [`qtechtheory.org`](https://quest.qtechtheory.org/download/) for help downloading the necessary compilers.


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


## Optimisation

QuEST's source code is careful to enable a myriad of optimisations such as [inlining](https://en.wikipedia.org/wiki/Inline_expansion), [loop unrolling](https://en.wikipedia.org/wiki/Loop_unrolling), [auto-vectorisation](https://en.wikipedia.org/wiki/Automatic_vectorization) and [cache optimisations](https://en.wikipedia.org/wiki/Cache_replacement_policies). To utilise them fully, we must instruct our compilers to enable them; like we might do with the [`-O3`](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html) flag when invoking a compiler like `gcc` directly.

On most platforms (with the exception of Windows), this is automatic with the commands above, but can otherwise be forced by specifying [`CMAKE_BUILD_TYPE`](https://cmake.org/cmake/help/latest/guide/tutorial/Packaging%20Debug%20and%20Release.html) at _configure time_:

```bash
# configure
cmake .. -D CMAKE_BUILD_TYPE=Release
```

When compiling on **Windows** however (using [Visual Studio](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#visual-studio-generators)), or otherwise using a "[_multi-config generator_](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html#other-generators)", we must always supply the build type at **_build time_** via [`config`](https://cmake.org/cmake/help/latest/manual/cmake.1.html#cmdoption-cmake-build-config):
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


## Linking

QuEST can be pre-compiled and later linked to other binaries, _or_ compiled directly alongside the user's source code. 
We focus on the latter use-case, common among scientists when writing simulation scripts. Users seeking to integrate QuEST into larger stacks are likely already familiar with linking libraries through CMake and should check out [`cmake.md`](/docs/cmake.md) directly.

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
```C++
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
```C++
/* otherfile.cpp */

#include <stdio.h>

void myfunc() {
    printf("hello world!\n");
}
```
simply separate them by `;` in `USER_SOURCE`, wrapped in `"`:
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


Additional flags needed by your files can be passed to the `C` and `C++` compilers, and the linker (respectively), at _configuration time_ via
- [`CMAKE_C_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)
- [`CMAKE_CXX_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_FLAGS.html)
- [`CMAKE_EXE_LINKER_FLAGS`](https://cmake.org/cmake/help/latest/variable/CMAKE_EXE_LINKER_FLAGS.html)

For example,
```bash
# configure
cmake .. -D CMAKE_C_FLAGS="-D MYMACRO=5" -D CMAKE_EXE_LINKER_FLAGS="-lm"
```

However, if your configuration is any more complicated or your source code requires different `C`/`C++` 
standards than the QuEST source code, you should consider separately compiling QuEST then linking it
to your source code as a library!


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


## Tests

### v4

To compile QuEST's latest unit and integration tests, use

```bash
# configure
cmake .. -D ENABLE_TESTING=ON

# build
cmake --build .
```
This will compile an executable which can be run directly (from within the `build` folder) via
```bash
./tests/tests
```
which accepts all of the [Catch2 CLI arguments](https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md)
```bash
./tests/tests applyHadamard
```
and which can be distributed (when compiled for such):
```bash
mpirun -np 8 ./tests/tests applyHadamard
```

Alternatively, the tests can be run through [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html) via either
```bash
ctest
```
```bash
make test
```
which will log each passing test alive, but alas cannot be deployed with distribution.

### v3

QuEST's deprecated v3 API has its own unit tests which can be additionally compiled (_except_ on Windows) via
```bash
# configure
cmake .. -D ENABLE_TESTING=ON -D ENABLE_DEPRECATED_API=ON

# build
cmake --build .
```
and run directly as a [Catch2 process](https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md) as
```bash
./tests/deprecated/dep_tests
```
or run through [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html) via
```bash
cd tests/deprecated
ctest
```

> [!CAUTION]
> The deprecated unit tests are non-comprehensive and the deprecated API should not be relied upon, for it may introduce
> undetectable corner-case bugs. Please only use the deprecated API and tests for assistance porting your application
> from QuEST v3 to v4.


## Multithreading

QuEST uses [OpenMP](https://www.openmp.org/) to perform multithreading, so accelerating QuEST over multiple CPUs or cores requires a compiler integrated with OpenMP. This is true of almost all major compilers.

> [!IMPORTANT]  
> Using [`Clang`](https://clang.llvm.org/) on MacOS requires use of the `libomp` library, obtainable via [Homebrew](https://brew.sh/):
> ```bash
> brew install libomp
> ```
> and which must be exposed to Clang prior to compilation via
> ```
> export OpenMP_ROOT=$(brew --prefix)/opt/libomp
> ```

Compiling with multithreading is as simple as enabling it during configuration
```bash
# configure
cmake .. -D ENABLE_MULTITHREADING=ON

# build
cmake --build .
```
which is in fact enabled by default!

The number of threads over which to parallelise QuEST's execution can be chosen _before launching_ the compiled executable via the [`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html) environment variable, e.g. either
```bash
OMP_NUM_THREADS=16 ./myexec
```
```bash
export OMP_NUM_THREADS=32
./myexec
```

<!-- the doxygen-doc hyperlink below includes a hash of the function name which should be unchanging! -->
One can verify the number of utilised threads by compiling and running a file which calls [`reportQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496), which outputs a subsection such as
```
  [cpu]
    numCpuCores.......10 per machine
    numOmpProcs.......10 per machine
    numOmpThrds.......32 per node
```


> TODO:
> - binding sockets



## GPU-acceleration

QuEST supports both NVIDIA GPUs (using CUDA) and AMD GPUs (using HIP). Using either requires obtaining a specialised compiler and passing some GPU-specific compiler flags.


### NVIDIA

> TODO!
> - CUDA and nvidia-smi and drivers eh
> - CUDA GPU
> - min CC
> - ENABLE_CUDA
> - CMAKE_CUDA_ARCHITECTURES

### AMD

> TODO!
> - ROCm
> - ENABLE_HIP
> - CMAKE_HIP_ARCHITECTURES

## cuQuantum

> TODO!
> - OS requirement
> - CC requirement
> - downloading/installing
> - resolving path (curoot?)

## Distribution

> TODO!
> - compiling
> - launching
> - power-of-2 nodes
> - distributing over sockets
> - oversubscribing

## GPUDirect

> TODO!
> - CUDA-aware MPI
> - UCX
> - launch flags
> - checking via reportenv
