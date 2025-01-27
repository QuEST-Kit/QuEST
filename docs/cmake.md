# CMake Configuration Options in QuEST

Version 4 of QuEST includes reworked CMake to support library builds, CMake export, and installation. Here we detail useful variables to configure the compilation of QuEST. If using a Unix-like operating system any of these variables can be set using the `-D` flag when invoking CMake, for example:

```
cmake -Bbuild -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/QuEST -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DENABLE_MULTITHREADING=ON -DENABLE_DISTRIBUTION=OFF ./
```

Then one need only move to the build directory, and invoke make:

```
cd build
make
```

## QuEST CMake variables

| Variable | (Default) Values | Notes |
| -------- | ---------------- | ----- |
| `LIB_NAME` | (`QuEST`), String | The QuEST library will be named `lib${LIB_NAME}.so`. Can be used to differentiate multiple versions of QuEST which have been compiled. |
| `VERBOSE_LIB_NAME` | (`OFF`), `ON` | When turned on `LIB_NAME` will be modified according to the other configuration options chosen. For example compiling QuEST with multithreading, distribution, and double precision with `VERBOSE_LIB_NAME` turned on creates `libQuEST-fp2+mt+mpi.so`. |
| `FLOAT_PRECISION` | (`2`), `1`, `4` | Determines which floating-point precision QuEST will use: double, single, or quad. *Note: Quad precision is not supported when also compiling for GPU.* |
| `BUILD_EXAMPLES` | (`ON`), `OFF` | Determines whether the example programs will be built alongside QuEST. |
| `ENABLE_TESTING` | (`ON`), `OFF` | Determines whether Catch2 tests will be built alongisde QuEST. If built, tests can be run from the build directory with `make test`. |
| `ENABLE_MULTITHREADING` | (`ON`), OFF | Determines whether QuEST will be built with support for parallelisation with OpenMP. |
| `ENABLE_DISTRIBUTION` | (`OFF`), ON | Determines whether QuEST will be built with support for parallelisation with MPI. |
| `ENABLE_CUDA` | (`OFF`), `ON` | Determines whether QuEST will be built with support for NVIDIA GPU acceleration. If turned on, `CMAKE_CUDA_ARCHITECTURES` should probably also be set. |
| `ENABLE_CUQUANTUM` | (`OFF`), `ON` | Determines whether QuEST will make use of the NVIDIA CuQuantum library. Cannot be turned on if `ENABLE_CUDA` is off. |
| `ENABLE_HIP` | (`OFF`), `ON` | Determines whether QuEST will be built with support for AMD GPU acceleration. If turned on, `CMAKE_HIP_ARCHITECTURES` should probably also be set. |
| `ENABLE_DEPRECATION` | (`OFF`), `ON` | Determines whether QuEST will be built with support for the deprecated (v3) API. *Note: will generate compiler warnings, and not supported by GCC.` |
| `USER_SOURCE` | (Undefined), String | The source file for a user program which will be compiled alongside QuEST. `OUTPUT_EXE` *must* also be defined. |
| `OUTPUT_EXE` | (Undefined), String | The name of the executable which will be created from the provided `USER_SOURCE`. `USER_SOURCE` *must* also be defined. |

## Standard CMake variables

| Variable | Description | CMake Doc Page |
| -------- | ----------- | ----- |
| `CMAKE_BUILD_TYPE` | Whether QuEST will be built with or without optimisations and debugging info. QuEST defaults to a `Release` build which is with optimisation and without debugging info. | [CMAKE_BUILD_TYPE](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html) |
| `CMAKE_CXX_COMPILER` | The C++ compiler that will be used to compile QuEST. | [CMAKE_\<LANG\>_COMPILER](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html) |
| `CMAKE_C_COMPILER` | The C compiler that will be used to compile QuEST. | [CMAKE_\<LANG\>_COMPILER](https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER.html) |
| `CMAKE_INSTALL_PREFIX` | The directory to which QuEST will be installed when `make install` is invoked. A standard GNU directory structure (lib, bin, include) will be used inside the prefix directory. | [CMAKE_INSTALL_PREFIX](https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html) <br> [GNUInstallDirs](https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html) |
| `CMAKE_CUDA_ARCHITECTURES` | Used to set the value of `arch` when compiling for NVIDIA GPU. | [CMAKE_CUDA_ARCHITECTURES](https://cmake.org/cmake/help/latest/variable/CMAKE_CUDA_ARCHITECTURES.html) |
| `CMAKE_HIP_ARCHITECTURES` | Used to set the HIP platform which QuEST is compiled for when compiling for AMD GPU. | [CMAKE_HIP_ARCHITECTURES](https://cmake.org/cmake/help/latest/variable/CMAKE_HIP_ARCHITECTURES.html) |