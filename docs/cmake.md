# ⚙️  CMake

<!--
  Instructions for compiling QuEST with CMake
  (this comment must be under the title for valid doxygen rendering)

  @author Oliver Thomson Brown
  @author Tyson Jones (test variables)
-->

QuEST requires CMake 3.28 or newer. Version 4 includes CMake support for library builds, installation, exported targets, and binary and source packages. Here we detail useful variables to configure the compilation of QuEST. Set any of these cache variables with the `-D` flag when invoking CMake, for example:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/QuEST -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DQUEST_ENABLE_OMP=ON -DQUEST_ENABLE_MPI=OFF
```

Then, as detailed in [`compile.md`](compile.md), compile through CMake:

```
cmake --build build
```

Install an install-enabled build with CMake's portable install command:

```bash
cmake --install build --config Release
```

> [!NOTE]
> Windows or MSVC users should additionally supply `--config Release` during compilation to enable optimisations.



------------------------

## QuEST variables

| Variable | (Default) Values | Notes |
| -------- | ---------------- | ----- |
| `QUEST_ENABLE_INSTALL` | (`ON` standalone, `OFF` as a subproject), `ON`, `OFF` | Enables QuEST installation and exported CMake package files. Installable ADIOS2 builds require an externally installed ADIOS2 package. |
| `QUEST_ENABLE_PACKAGING` | (`ON` for a standalone installable build, otherwise `OFF`), `ON`, `OFF` | Enables CPack configuration. Packaging requires installation. |
| `QUEST_OUTPUT_LIB_NAME` | (`QuEST`), String | Changes the library artifact name. The installed CMake package and target remain `QuESTConfig.cmake` and `QuEST::QuEST`. |
| `QUEST_APPEND_CONFIG_TO_LIB_NAME` | (`OFF`), `ON` | When turned on `QUEST_OUTPUT_LIB_NAME` will be modified according to the other configuration options chosen. For example compiling QuEST with multithreading, distribution, and double precision with `QUEST_APPEND_CONFIG_TO_LIB_NAME` turned on creates `libQuEST-fp2+mt+mpi.so`. |
| `QUEST_FLOAT_PRECISION` | (`2`), `1`, `4` | Determines which floating-point precision QuEST will use: double, single, or quad. *Note: Quad precision is not supported when also compiling for GPU.* |
| `QUEST_BUILD_MIN_EXAMPLE` | (`ON` standalone, `OFF` as a subproject), `ON`, `OFF` | Determines whether the minimum example is built. |
| `QUEST_BUILD_EXAMPLES` | (`OFF`), `ON` | Determines whether the other example programs are built alongside QuEST. |
| `QUEST_INSTALL_BINARIES` | (`OFF`), `ON` | Determines whether compiled binaries such as the examples will be installed as well as the QuEST library. |
| `QUEST_ENABLE_OMP` | (`ON`), `OFF` | Determines whether QuEST will be built with support for parallelisation with OpenMP. |
| `QUEST_ENABLE_NUMA` | (`ON`), `OFF` | Determines whether QuEST will attempt to build with NUMA awareness when OpenMP is also enabled. |
| `QUEST_ENABLE_MPI` | (`OFF`), `ON` | Determines whether QuEST will be built with support for parallelisation with MPI. |
| `QUEST_ENABLE_SUBCOMM` | (`OFF`), `ON` | Determines whether QuEST will be built with support for custom MPI communicators. _**Note**: This has the unfortunate side-effect of requiring the MPI header in the public header for QuEST, meaning MPI will become a dependency of any application or library which includes the QuEST header whether it uses MPI or not._ |
| `QUEST_ENABLE_CUDA` | (`OFF`), `ON` | Determines whether QuEST will be built with support for NVIDIA GPU acceleration. If turned on, `CMAKE_CUDA_ARCHITECTURES` should probably also be set. |
| `QUEST_ENABLE_CUQUANTUM` | (`OFF`), `ON` | Determines whether QuEST will make use of the NVIDIA CuQuantum library. Cannot be turned on if `QUEST_ENABLE_CUDA` is off. |
| `QUEST_ENABLE_HIP` | (`OFF`), `ON` | Determines whether QuEST will be built with support for AMD GPU acceleration. If turned on, `CMAKE_HIP_ARCHITECTURES` should probably also be set. |
| `QUEST_ENABLE_BMI2` | (`OFF`), `ON` | Determines whether QuEST will be built with BMI2 intrinsics to accelerate CPU simulation of few-qubit Quregs. This is not compatible with all compilers and CPUs. **Beware** that if enabled, and the compiled QuEST executable is later run upon a different machine which lacks the BMI2 instructions, execution will crash. |
| `QUEST_ENABLE_ADIOS2` | (`OFF`), `ON` | Determines whether QuEST will be built with ADIOS2 to enable checkpointing, via functions `saveQuregToFile()` and `createQuregFromFile()`. |
| `QUEST_DOWNLOAD_ADIOS2` | (`ON`), `OFF` | Determines whether to download ADIOS2 from GitHub when ADIOS2 is enabled but not found. Downloading is available only when `QUEST_ENABLE_INSTALL=OFF`; installable builds must use an external compatible ADIOS2 package. |
| `QUEST_ENABLE_DEPRECATED_API` | (`OFF`), `ON` | Determines whether QuEST will be built with support for the deprecated (v3) API. ***Note**: this will generate compiler warnings and is not supported by MSVC.* |
| `QUEST_DISABLE_DEPRECATION_WARNINGS` | (`OFF`), `ON` | Whether to disable the compile-time deprecation warnings when using the deprecated (v3) API. |
| `USER_SOURCE_NAMES` | (Undefined), String | The source file for a user program which will be compiled alongside QuEST. `USER_OUTPUT_EXE_NAME` *must* also be defined. |
| `USER_OUTPUT_EXE_NAME` | (Undefined), String | The name of the executable which will be created from the provided `USER_SOURCE_NAMES`. `USER_SOURCE_NAMES` *must* also be defined. |
| `QUEST_DEFAULT_NUM_GPU_THREADS_PER_BLOCK` | (128), Number | The default number of threads per block QuEST will use when offloading to a GPU. *Must* be a multiple of 32 (on NVIDIA GPUs) or 64 (on AMD GPUs). This CMake variable sets the default if not later overridden. The number can be overridden at process launch time using an [environment variable](https://quest-kit.github.io/QuEST/group__modes.html#gaf1b71f54d270d3353fe072c66827339b) of the same name, or during runtime using [`setQuESTNumGpuThreadsPerBlock()`](https://quest-kit.github.io/QuEST/group__experimental.html#gae35a55c6d9366ce677e6aaaf4c1ff5ef). |



--------------------------

## Test variables

| Variable | (Default) Values | Notes |
| -------- | ---------------- | ----- |
| `QUEST_BUILD_TESTS` | (`OFF`), `ON` | Determines whether to additionally build QuEST's unit and integration tests. If built, tests can be run from the `build` directory with `make test`, or `ctest`, or manually launched with `./tests/tests` which enables distribution (i.e. `mpirun -np 8 ./tests/tests`) |
| `QUEST_BUILD_PACKAGING_TESTS` | (`OFF`), `ON` | Builds installation, relocation, exported-target, finder, and packaging checks. These tests do not require Catch2. Run them with `ctest --test-dir build -L packaging --output-on-failure`. |
| `QUEST_ENABLE_DEPRECATED_API` | (`OFF`), `ON` | As described above. When enabled alongside testing, the `v3 deprecated` unit tests will additionally be compiled and can be run from within `build` via `cd tests/deprecated; ctest`, or manually launched with `./tests/deprecated/dep_tests` (enabling distribution, as above). |
| `QUEST_TESTS_DOWNLOAD_CATCH2` | (`ON`), `OFF` | QuEST's tests require Catch2. By default, if you don't have Catch2 installed (or CMake doesn't find it) it will be downloaded from Github and built for you. If you don't want that to happen, for example because you _do_ have Catch2 installed, set this to `OFF`. |

> As of `v4.2`, macros which configure the unit tests such as `QUEST_TEST_MAX_NUM_QUBIT_PERMUTATIONS` have become environment variables specified before launch. See [`launch.md`](launch.md)

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
| `CMAKE_RUNTIME_OUTPUT_DIRECTORY` | The output directory to which to save compiled executables, overriding the default `build` folder | [`CMAKE_RUNTIME_OUTPUT_DIRECTORY`](https://cmake.org/cmake/help/latest/variable/CMAKE_RUNTIME_OUTPUT_DIRECTORY.html). |


---------------------------

## Using an installed QuEST

An install always publishes the canonical `QuESTConfig.cmake`, `QuESTConfigVersion.cmake`, and `QuESTTargets.cmake` files. `QUEST_OUTPUT_LIB_NAME` and `QUEST_APPEND_CONFIG_TO_LIB_NAME` change only the library artifact name. One QuEST configuration is supported per installation prefix.

Downstream projects need only discover the package and link its exported target:

```cmake
cmake_minimum_required(VERSION 3.28)
project(my_quest_program LANGUAGES C CXX)

find_package(QuEST CONFIG REQUIRED)
add_executable(my_quest_program main.c)
target_link_libraries(my_quest_program PRIVATE QuEST::QuEST)
```

Enabling both C and CXX is supported for a C application and permits CMake to satisfy a static QuEST library's C++ linker requirements. The exported target requests C11 for C consumers and C++14 for C++ consumers. QuEST's C++17 implementation and GPU C++20 requirements remain private build requirements. Installed GPU packages can be consumed without enabling CUDA or HIP as project languages.

Configure the consumer with the QuEST prefix when it is outside CMake's normal search locations:

```bash
cmake -S consumer -B consumer-build -DCMAKE_PREFIX_PATH=/opt/QuEST
cmake --build consumer-build
```

The installed configuration rediscovers the dependencies required by the built library before loading `QuEST::QuEST`. It does not consult downstream `QUEST_*` settings or `BUILD_SHARED_LIBS` to reinterpret the installed binary. Static packages can therefore require development packages for enabled OpenMP, NUMA, MPI, CUDA, HIP, cuQuantum, or ADIOS2 backends. MPI is also a public requirement when QuEST was built with the subcommunicator API because its public header exposes `mpi.h`. Shared packages preserve their external runtime requirements while avoiding private SDK development requirements where the link interface does not need them.

ADIOS2 builds select the serial or MPI C++ target to match QuEST's MPI configuration. Current ADIOS2 packages normally provide `adios2::cxx` and `adios2::cxx_mpi`; ADIOS2 2.9 packages, including Ubuntu 24.04, use the compatible legacy names `adios2::cxx11` and `adios2::cxx11_mpi`. QuEST records the selected target and requires the same interface when a static installation is consumed.

QuEST installs its reusable cuQuantum and cuTENSOR find modules with the package. When a static built library uses cuQuantum, `QuESTConfig.cmake` makes this dependency request before loading the exported target:

```cmake
find_package(CUQUANTUM MODULE REQUIRED COMPONENTS cuStateVec)
```

Consumers still call only `find_package(QuEST CONFIG REQUIRED)`. The static package records the producer's cuStateVec component version and requires an equal or newer version with the same major ABI. A shared QuEST package retains its cuQuantum runtime requirement without resolving private SDK development files during consumer configuration.

The available cuQuantum imported targets are `CUQUANTUM::cuStateVec`, `CUQUANTUM::cuTensorNet`, and `CUQUANTUM::cuDensityMat`. Calling the finder without components requests all three. Each component's header and shared library are resolved independently; the finder does not substitute the SDK's `_static` archives. Set `CUQUANTUM_ROOT` or its environment variable to the SDK prefix and `CUDAToolkit_ROOT` to CUDA. cuTensorNet and cuDensityMat also require cuTENSOR and accept `CUTENSOR_ROOT`. An explicitly set `CUQUANTUM_DIR` remains a legacy prefix hint with precedence over these general roots. Component versions are reported separately as `CUQUANTUM_<component>_VERSION`; they are not the overall SDK release version.


---------------------------

## Creating packages

For a standalone installable build, CPack is enabled by default after the installation rules. TGZ and ZIP produce complete binary archives, and the source CPack configuration produces complete source archives:

```bash
cmake -S . -B build -DQUEST_ENABLE_PACKAGING=ON
cmake --build build --config Release
cpack --config build/CPackConfig.cmake -G TGZ
cpack --config build/CPackConfig.cmake -G ZIP
cpack --config build/CPackSourceConfig.cmake -G TGZ
cpack --config build/CPackSourceConfig.cmake -G ZIP
```

Binary archive names record the QuEST version, platform, architecture, configuration, shared or static linkage, precision, and enabled backends. Packages contain only QuEST-owned files and combine the QuEST components into one complete archive. The install components are `Runtime`, `Development`, and, when installable examples were built, `Examples`. `Development` contains headers, CMake exports and find modules, static or import libraries, and linker namelinks. A shared `Development` package depends on the exact `Runtime` version; a static build has no empty runtime dependency.

Native package profiles are available with `-DQUEST_NATIVE_PACKAGE_PROFILE=ubuntu24.04` for DEB and `-DQUEST_NATIVE_PACKAGE_PROFILE=fedora44` for RPM. Supported hosts are detected when the variable is empty; use `custom` for an explicitly described vendor environment. Native profiles install under `/usr` and use GNU installation directories. They produce `libquest4`, `libquest-dev`, and optional `quest-examples` packages on Debian, or `quest`, `quest-devel`, and optional `quest-examples` packages on Fedora.

Fedora's stock Open MPI installation is module-based. Load it when configuring an MPI-enabled Fedora package and when running or building consumers of that package:

```bash
source /etc/profile.d/modules.sh
module load mpi/openmpi-x86_64
```

The stock profiles describe GCC, OpenMP, NUMA, and Open MPI dependencies and enable the platform's shared-library dependency scanner. For ADIOS2 or GPU variants, set complete native dependency metadata explicitly with the standard per-component CPack variables, such as `CPACK_DEBIAN_DEVELOPMENT_PACKAGE_DEPENDS`, `CPACK_DEBIAN_RUNTIME_PACKAGE_DEPENDS`, and `CPACK_DEBIAN_EXAMPLES_PACKAGE_DEPENDS`, or the corresponding `CPACK_RPM_*_PACKAGE_REQUIRES` variables. The same requirement applies to custom profiles and non-GCC native builds.
