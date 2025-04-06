# ⚙️  CMake

<!--
  Instructions for compiling QuEST with CMake
  (this comment must be under the title for valid doxygen rendering)

  @author Oliver Thomson Brown
  @author Tyson Jones (test variables)
-->

Version 4 of QuEST includes reworked CMake to support library builds, CMake export, and installation. Here we detail useful variables to configure the compilation of QuEST. If using a Unix-like operating system any of these variables can be set using the `-D` flag when invoking CMake, for example:

```
cmake -Bbuild -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/QuEST -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DENABLE_MULTITHREADING=ON -DENABLE_DISTRIBUTION=OFF ./
```

Then one need only move to the build directory, and invoke make:

```
cd build
make
```

> [!NOTE]
> Windows or MSVC users should additionally supply `--config Release` during compilation to enable optimisations.



------------------------

## QuEST variables

| Variable | (Default) Values | Notes |
| -------- | ---------------- | ----- |
| `LIB_NAME` | (`QuEST`), String | The QuEST library will be named `lib${LIB_NAME}.so`. Can be used to differentiate multiple versions of QuEST which have been compiled. |
| `VERBOSE_LIB_NAME` | (`OFF`), `ON` | When turned on `LIB_NAME` will be modified according to the other configuration options chosen. For example compiling QuEST with multithreading, distribution, and double precision with `VERBOSE_LIB_NAME` turned on creates `libQuEST-fp2+mt+mpi.so`. |
| `FLOAT_PRECISION` | (`2`), `1`, `4` | Determines which floating-point precision QuEST will use: double, single, or quad. *Note: Quad precision is not supported when also compiling for GPU.* |
| `BUILD_EXAMPLES` | (`OFF`), `ON` | Determines whether the example programs will be built alongside QuEST. Note that `min_example` is always built. |
| `ENABLE_MULTITHREADING` | (`ON`), OFF | Determines whether QuEST will be built with support for parallelisation with OpenMP. |
| `ENABLE_DISTRIBUTION` | (`OFF`), ON | Determines whether QuEST will be built with support for parallelisation with MPI. |
| `ENABLE_CUDA` | (`OFF`), `ON` | Determines whether QuEST will be built with support for NVIDIA GPU acceleration. If turned on, `CMAKE_CUDA_ARCHITECTURES` should probably also be set. |
| `ENABLE_CUQUANTUM` | (`OFF`), `ON` | Determines whether QuEST will make use of the NVIDIA CuQuantum library. Cannot be turned on if `ENABLE_CUDA` is off. |
| `ENABLE_HIP` | (`OFF`), `ON` | Determines whether QuEST will be built with support for AMD GPU acceleration. If turned on, `CMAKE_HIP_ARCHITECTURES` should probably also be set. |
| `ENABLE_DEPRECATED_API` | (`OFF`), `ON` | Determines whether QuEST will be built with support for the deprecated (v3) API. ***Note**: this will generate compiler warnings and is not supported by MSVC.* |
| `USER_SOURCE` | (Undefined), String | The source file for a user program which will be compiled alongside QuEST. `OUTPUT_EXE` *must* also be defined. |
| `OUTPUT_EXE` | (Undefined), String | The name of the executable which will be created from the provided `USER_SOURCE`. `USER_SOURCE` *must* also be defined. |


--------------------------

## Test variables

| Variable | (Default) Values | Notes |
| -------- | ---------------- | ----- |
| `ENABLE_TESTING` | (`OFF`), `ON` | Determines whether to additionally build QuEST's unit and integration tests. If built, tests can be run from the `build` directory with `make test`, or `ctest`, or manually launched with `./tests/tests` which enables distribution (i.e. `mpirun -np 8 ./tests/tests`) |
| `ENABLE_DEPRECATED_API` | (`OFF`), `ON` | As described above. When enabled alongside testing, the `v3 deprecated` unit tests will additionally be compiled and can be run from within `build` via `cd tests/deprecated; ctest`, or manually launched with `./tests/deprecated/dep_tests` (enabling distribution, as above).
| `DOWNLOAD_CATCH2` | (`ON`), `OFF` | QuEST's tests require Catch2. By default, if you don't have Catch2 installed (or CMake doesn't find it) it will be downloaded from Github and built for you. If you don't want that to happen, for example because you _do_ have Catch2 installed, set this to `OFF`. |
| `TEST_MAX_NUM_QUBIT_PERMUTATIONS` | (`0`), Integer | Determines the maximum number of control and target qubit permutations with which to test each API function. Set to `0` to test all permutations, or to a positive integer (e.g. `50`) to accelerate the unit tests. |
| `TEST_MAX_NUM_SUPEROP_TARGETS` | (`4`), Integer | Determines the maximum number of superoperator targets in the unit tests (for functions `mixKrausMap` and `mixSuperOp`). Set to `0` to impose no maximum (which is extraordinarily slow), or a positive integer (e.g. `3`) to accelerate the unit tests. |
| `TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS` | (`10`), Integer | Determines the number of times (minimum of `1`) to repeat the randomised unit tests of functions which accept multiple distinctly-deployed `Qureg`s. Set to a small, positive integer to accelerate mixed-deployment unit tests. |
| `TEST_ALL_DEPLOYMENTS` | (`ON`), `OFF` | Determines whether unit tests will be repeatedly run using all possible combinations of available `Qureg` deployments (i.e. `OpenMP` and/or `CUDA` and/or `MPI`), else only once using all available deployments simultaneously. Set to `OFF` to accelerate unit tests. |

---------------------------

## Standard CMake variables

| Variable | Description | CMake Doc Page |
| -------- | ----------- | ----- |
| `CMAKE_BUILD_TYPE` | Whether QuEST will be built with or without optimisations and debugging info. QuEST defaults to a `Release` build which is with optimisation and without debugging info. | [CMAKE_BUILD_TYPE](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html) |
| `CMAKE_CXX_COMPILER` | The C++ compiler that will be used to compile QuEST. | [CMAKE_\<LANG\>_COMPILER](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html) |
| `CMAKE_C_COMPILER` | The C compiler that will be used to compile QuEST. | [CMAKE_\<LANG\>_COMPILER](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html) |
| `CMAKE_INSTALL_PREFIX` | The directory to which QuEST will be installed when `make install` is invoked. A standard GNU directory structure (lib, bin, include) will be used inside the prefix directory. | [CMAKE_INSTALL_PREFIX](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html) <br> [GNUInstallDirs](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html) |
| `CMAKE_CUDA_ARCHITECTURES` | Used to set the value of `arch` when compiling for NVIDIA GPU. This is also known as the target GPU's "compute capability" and can be discovered [here](https://developer.nvidia.com/cuda-gpus). | [CMAKE_CUDA_ARCHITECTURES](https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html) |
| `CMAKE_HIP_ARCHITECTURES` | Used to set the HIP platform which QuEST is compiled for when compiling for AMD GPU. | [CMAKE_HIP_ARCHITECTURES](https://cmake.org/cmake/help/latest/variable/CMAKE_HIP_ARCHITECTURES.html) |