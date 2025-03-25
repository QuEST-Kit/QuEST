# @author Oliver Thomson Brown
# @author Erich Essmann (patches including MSVC support)
# @author Tyson Jones (patches including clang multithreading)
#
# Contributions to previous builds from:
#  - Ania Brown
#  - Jacob Wilkins
#  - Balint Koczor
#  - Richard Meister
#  - Gleb Struchalin
#  - Sachin Compton
#  - Christopher J. Anders
#  - Drew Silcock

cmake_minimum_required(VERSION 3.21)

project(QuEST
  VERSION 4.0.0
  DESCRIPTION "Quantum Exact Simulation Toolkit"
  LANGUAGES CXX C
)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

## Dependencies

# GNUInstallDirs to provide sensible default install directory names
cmake_path(SET ORG_INSTALL_PATH NORMALIZE "${CMAKE_INSTALL_PREFIX}")
cmake_path(APPEND CMAKE_INSTALL_PREFIX "quest")
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Maths
if (NOT WIN32)
  find_library(MATH_LIBRARY m REQUIRED)
endif()

## Configuration options

# Build type
# Default to "Release"
# Using recipe from Kitware Blog post
# https://www.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Library type
# Shared library by default
option(BUILD_SHARED_LIBS "Build shared library. Turned ON by default." ON)
message(STATUS "Shared library is turned ${BUILD_SHARED_LIBS}. Set BUILD_SHARED_LIBS to modify.")

# Library naming
set(LIB_NAME QuEST 
  CACHE 
  STRING
  "Change library name. LIB_NAME is QuEST by default."
)
message(STATUS "Library will be named lib${LIB_NAME}. Set LIB_NAME to modify.")

option(VERBOSE_LIB_NAME "Modify library name based on compilation configuration. Turned OFF by default." OFF)
message(STATUS "Verbose library naming is turned ${VERBOSE_LIB_NAME}. Set VERBOSE_LIB_NAME to modify.")

# Precision
set(FLOAT_PRECISION 2 
  CACHE 
  STRING 
  "Whether to use single, double, or quad floating point precision in the state vector. {1,2,4}"
)
set_property(CACHE FLOAT_PRECISION PROPERTY STRINGS
  1
  2
  4
)
message(STATUS "Precision set to ${FLOAT_PRECISION}. Set FLOAT_PRECISION to modify.")

if (VERBOSE_LIB_NAME)
  string(CONCAT LIB_NAME ${LIB_NAME} "-fp${FLOAT_PRECISION}")
endif()

# Examples
option(
  BUILD_EXAMPLES
  "Whether the example programs will be built alongside the QuEST library. Turned OFF by default."
  OFF
)
message(STATUS "Examples are turned ${BUILD_EXAMPLES}. Set BUILD_EXAMPLES to modify.")

# Testing
option(
  ENABLE_TESTING
  "Whether the test suite will be built alongside the QuEST library. Turned ON by default."
  OFF
)
message(STATUS "Testing is turned ${ENABLE_TESTING}. Set ENABLE_TESTING to modify.")

option(
  DOWNLOAD_CATCH2
  "Whether Catch2 v3 will be downloaded if it is not found. Turned ON by default."
  ON
)

string(CONCAT DESCRIPTON_OF_TEST_MAX_NUM_QUBIT_PERMUTATIONS  # for multiline description
  "The maximum number of control and target qubit permutations under which to unit test each function. "
  "Set to 0 (default) to test all permutations, or to a positive integer (e.g. 50) to accelerate the unit tests. "
  "This is used to accelerate the v4 unit tests, and has no effect when ENABLE_TESTING is OFF, nor on the v3 deprecated tests."
)
set(TEST_MAX_NUM_QUBIT_PERMUTATIONS 0
  CACHE 
  STRING 
  "${DESCRIPTON_OF_TEST_MAX_NUM_QUBIT_PERMUTATIONS}"
)
if (ENABLE_TESTING)
  set(PERM_STRING "${TEST_MAX_NUM_QUBIT_PERMUTATIONS} random")
  if (TEST_MAX_NUM_QUBIT_PERMUTATIONS EQUAL 0)
    set(PERM_STRING "all")
  endif()
  message(STATUS "Tests will use ${PERM_STRING} qubit permutations. Set TEST_MAX_NUM_QUBIT_PERMUTATIONS to modify.")
endif()

string(CONCAT DESCRIPTON_OF_TEST_MAX_NUM_SUPEROP_TARGETS  # for multiline description
  "The maximum number of superoperator targets for which to unit test functions mixKrausMap and mixSuperOp. "
  "This is computationally equivalent to simulating unitaries with double the number of targets upon a density matrix. "
  "Set to 0 to test all sizes which is likely prohibitively slow, or to a positive integer (e.g. the default of 4) to accelerate the unit tests. "
  "This is used to accelerate the v4 unit tests, and has no effect when ENABLE_TESTING is OFF, nor on the v3 deprecated tests."
)
set(TEST_MAX_NUM_SUPEROP_TARGETS 4
  CACHE 
  STRING 
  "${DESCRIPTON_OF_TEST_MAX_NUM_SUPEROP_TARGETS}"
)
if (ENABLE_TESTING)
  set(PERM_STRING "${TEST_MAX_NUM_SUPEROP_TARGETS} random")
  if (TEST_MAX_NUM_SUPEROP_TARGETS EQUAL 0)
    set(PERM_STRING "all")
  endif()
  message(STATUS "Tests will use superoperators of up to ${PERM_STRING} qubits. Set TEST_MAX_NUM_SUPEROP_TARGETS to modify.")
endif()

set(TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS 10
  CACHE 
  STRING 
  "The number of times (minimum=1) to repeat each mixed-deployment unit test for each deployment combination."
)
if (ENABLE_TESTING)
  if (TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS LESS_EQUAL 0)
    message(FATAL_ERROR "Cannot set TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS <= 0. Must be at least 1.")
  endif()
  message(STATUS "Tests will repeat each mixed-deployment unit test ${TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS} times. Set TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS to modify.")
endif()

string(CONCAT DESCRIPTION_OF_TEST_ALL_DEPLOYMENTS  # for multiline description
  "Whether unit tests will be run using all possible deployment combinations (i.e. OpenMP, CUDA, MPI) in-turn (ON), "
  "or only once using all available deployments simultaneously (OFF). Turned ON by default. "
  "This is used to accelerate the v4 unit tests, and has no effect when ENABLE_TESTING is OFF, nor on the v3 deprecated tests."
)
option(TEST_ALL_DEPLOYMENTS "${DESCRIPTION_OF_TEST_ALL_DEPLOYMENTS}" ON)
if (ENABLE_TESTING)
  message(STATUS "Testing all deployments is turned ${TEST_ALL_DEPLOYMENTS}. Set TEST_ALL_DEPLOYMENTS to modify.")
endif()

# Multithreading
option(
  ENABLE_MULTITHREADING 
  "Whether QuEST will be built with shared-memory parallelism support using OpenMP. Turned ON by default." 
  ON
)
message(STATUS "Multithreading is turned ${ENABLE_MULTITHREADING}. Set ENABLE_MULTITHREADING to modify.")

# Distribution
option(
  ENABLE_DISTRIBUTION 
  "Whether QuEST will be built with distributed parallelism support using MPI. Turned OFF by default." 
  OFF
)
message(STATUS "Distribution is turned ${ENABLE_DISTRIBUTION}. Set ENABLE_DISTRIBUTION to modify.")

# GPU Acceleration
option(
  ENABLE_CUDA
  "Whether QuEST will be built with support for NVIDIA GPU acceleration. Turned OFF by default."
  OFF
)
message(STATUS "NVIDIA GPU acceleration is turned ${ENABLE_CUDA}. Set ENABLE_CUDA to modify.")

if (ENABLE_CUDA)
  option(
    ENABLE_CUQUANTUM
    "Whether QuEST will be built with support for NVIDIA CuQuantum. Turned OFF by default."
    OFF
  )
  message(STATUS "CuQuantum support is turned ${ENABLE_CUQUANTUM}. Set ENABLE_CUQUANTUM to modify.")
endif()

option(
  ENABLE_HIP
  "Whether QuEST will be built with support for AMD GPU acceleration. Turned OFF by default."
  OFF
)
message(STATUS "AMD GPU acceleration is turned ${ENABLE_HIP}. Set ENABLE_HIP to modify.")

# Throw on disallowed combinations
if (ENABLE_CUDA AND ENABLE_HIP)
  message(FATAL_ERROR "QuEST cannot support CUDA and HIP simultaneously.")
endif()

if ((ENABLE_CUDA OR ENABLE_HIP) AND FLOAT_PRECISION STREQUAL 4)
  message(FATAL_ERROR "Quad precision is not supported on GPU. Please disable GPU acceleration or lower precision.")
endif()

option(
  PERMIT_NODES_TO_SHARE_GPU
  "Whether to permit multiple distributed nodes to share a single GPU at the detriment of performance. Turned OFF by default."
  OFF
)
if (ENABLE_DISTRIBUTION AND (ENABLE_CUDA OR ENABLE_HIP))
  message(STATUS "Permitting nodes to share GPUs is turned ${PERMIT_NODES_TO_SHARE_GPU}. Set PERMIT_NODES_TO_SHARE_GPU to modify.")
endif()

# Deprecated API
option(
  ENABLE_DEPRECATED_API
  "Whether QuEST will be built with deprecated API support. Turned OFF by default."
  OFF
)
message(STATUS "Deprecated API support is turned ${ENABLE_DEPRECATED_API}. Set ENABLE_DEPRECATED_API to modify.")

# Windows Specific Options
if(WIN32)
  
  # Force MSVC to export all symbols in a shared library, like GCC and clang
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  if (ENABLE_TESTING)
    set(BUILD_SHARED_LIBS OFF)
  endif()
endif()

## Library

add_library(QuEST)

# Add namespaced alias to support inclusion of QuEST as a subproject
add_library(QuEST::QuEST ALIAS QuEST)

# Set include directories
target_include_directories(QuEST
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quest/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:../quest/include>
  $<INSTALL_INTERFACE:..>
)

# Add required C and C++ standards
target_compile_features(QuEST
  PUBLIC
  c_std_11
  cxx_std_20
)

# Turn on all compiler warnings
if (MSVC)
  set(WARNING_FLAG /W4)
else()
  set(WARNING_FLAG -Wall)
endif()

target_compile_options(QuEST 
  PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:${WARNING_FLAG}>
  $<$<COMPILE_LANGUAGE:C>:${WARNING_FLAG}>
)


# Set user options
target_compile_definitions(QuEST PUBLIC FLOAT_PRECISION=${FLOAT_PRECISION})

if (ENABLE_MULTITHREADING)

  # find OpenMP, but fail gracefully...
  find_package(OpenMP QUIET)

  # so that we can customise the error message on MacOS
  if (NOT OpenMP_FOUND)
    set(ErrorMsg "Could not find OpenMP, necessary for enabling multithreading.")
    if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      string(APPEND ErrorMsg " Try first calling `brew install libomp` then `export OpenMP_ROOT=$(brew --prefix)/opt/libomp`")
    endif()
    message(FATAL_ERROR ${ErrorMsg})
  endif()

  target_compile_definitions(QuEST PUBLIC COMPILE_OPENMP=1)
  target_link_libraries(QuEST
    PRIVATE
    OpenMP::OpenMP_CXX
    OpenMP::OpenMP_C
  )

  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+mt")
  endif()

else()

  # suppress GCC "unknown pragma" warning when OpenMP disabled
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(QuEST PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>)
  endif()
  
  target_compile_definitions(QuEST PUBLIC COMPILE_OPENMP=0)
endif()

if (ENABLE_DISTRIBUTION)
  find_package(MPI REQUIRED
    COMPONENTS CXX
  )
  target_compile_definitions(QuEST PUBLIC COMPILE_MPI=1)
  target_link_libraries(QuEST
    PRIVATE
    MPI::MPI_CXX
  )
  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+mpi")
  endif()
else()
  target_compile_definitions(QuEST PUBLIC COMPILE_MPI=0)
endif()

if (ENABLE_CUDA)

  # make nvcc use user cxx-compiler as default host (before cuda-host is set below)
  if (NOT DEFINED CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  endif()

  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  
  set_property(TARGET QuEST PROPERTY CUDA_STANDARD 20)
  
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  
  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+cuda")
  endif()
endif()

if (ENABLE_CUQUANTUM)
  find_package(CUQUANTUM REQUIRED)
  target_compile_definitions(QuEST PUBLIC COMPILE_CUQUANTUM=1)
  target_link_libraries(QuEST PRIVATE CUQUANTUM::cuStateVec)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+cuquantum")
  endif()
else()
  target_compile_definitions(QuEST PUBLIC COMPILE_CUQUANTUM=0)
endif()

if (ENABLE_HIP)

  # if generation fails (hip::amdhip64 not found), users can try setting
  # CMAKE_MODULE_PATH to '/opt/rocm/cmake' or '/opt/rocm/hip/lib/cmake/hip'
  # (suitable when shared library libamdhip64.so is located in /opt/rocm/lib/
  #  or /opt/rocm/hip/lib/ respectively). We avoid setting CMAKE_MODULE_PATH
  # pre-emptively since it made successful generation less likely in our tests!
  # example: list(APPEND CMAKE_MODULE_PATH "/opt/rocm/cmake"). Users should
  # also add '/opt/rocm/bin' or '/opt/rocm/hip/bin' to their $PATH env-var.

  enable_language(HIP)
  set(CMAKE_HIP_STANDARD_REQUIRED ON)
  set_property(TARGET QuEST PROPERTY HIP_STANDARD 20)

  find_package(HIP REQUIRED)
  message(STATUS "Found HIP: " ${HIP_VERSION})

  target_compile_definitions(QuEST PUBLIC COMPILE_CUQUANTUM=0)
  target_link_libraries(QuEST PRIVATE hip::host)

  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+hip")
  endif()
endif()

if (ENABLE_CUDA OR ENABLE_HIP)
  target_compile_definitions(QuEST PUBLIC COMPILE_CUDA=1)
else()
  target_compile_definitions(QuEST PUBLIC COMPILE_CUDA=0)
endif()

if (ENABLE_DEPRECATED_API)
  target_compile_definitions(QuEST PRIVATE INCLUDE_DEPRECATED_FUNCTIONS=1)

  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+depr")
  endif()

else()
  target_compile_definitions(QuEST PRIVATE INCLUDE_DEPRECATED_FUNCTIONS=0)
endif()


if (ENABLE_DISTRIBUTION AND (ENABLE_CUDA OR ENABLE_HIP))
  target_compile_definitions(
    QuEST PRIVATE 
    PERMIT_NODES_TO_SHARE_GPU=$<IF:$<BOOL:${PERMIT_NODES_TO_SHARE_GPU}>,1,0>
  )
endif()


# add math library
if (NOT MSVC)
  target_link_libraries(QuEST PRIVATE ${MATH_LIBRARY})
endif()

# Set output name
set_target_properties(QuEST PROPERTIES OUTPUT_NAME ${LIB_NAME})

# Add source files
add_subdirectory(quest)

## Examples

add_executable(min_example
  examples/tutorials/min_example.cpp
)
target_link_libraries(min_example PRIVATE QuEST::QuEST)

install(TARGETS min_example
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)
set_target_properties(min_example
  PROPERTIES 
    INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
)

if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

## User Source

if (USER_SOURCE AND NOT OUTPUT_EXE)
    message(SEND_ERROR "USER_SOURCE specified, but not OUTPUT_EXE.")
endif()
if (OUTPUT_EXE AND NOT USER_SOURCE)
    message(SEND_ERROR "OUTPUT_EXE specified, but not USER_SOURCE.")
endif()
if (USER_SOURCE AND OUTPUT_EXE)
  message(STATUS "Compiling ${USER_SOURCE} to executable ${OUTPUT_EXE}.")

  add_executable(${OUTPUT_EXE} ${USER_SOURCE})
  target_link_libraries(${OUTPUT_EXE} PUBLIC QuEST)
  install(TARGETS ${OUTPUT_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  set_target_properties(${OUTPUT_EXE} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif()

## Tests

if (ENABLE_TESTING)

  set(CatchVersion 3.8.0)

  find_package(Catch2 ${CatchVersion} QUIET)

  if (NOT TARGET Catch2::Catch2 AND DOWNLOAD_CATCH2)
    message(STATUS "Catch2 not found, it will be downloaded and built in the build directory.")
    Include(FetchContent)

    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v${CatchVersion}
    )
    
    FetchContent_MakeAvailable(Catch2)
  
  else()
    # We won't magically find it here, but this is the easiest way to
    # a) Force the build to fail, and
    # b) Print out all the useful information for helping CMake find Catch2
    find_package(Catch2 ${CatchVersion} REQUIRED)
  endif()

  include(Catch)
  enable_testing()
  add_subdirectory(tests)
endif()

## Install

install(TARGETS QuEST
        EXPORT QuESTTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

## Export
# Write CMake version file for QuEST
set(QuEST_INSTALL_CONFIGDIR "${CMAKE_INSTALL_LIBDIR}/cmake/QuEST")

# Write QuESTConfigVersion.cmake
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

# Configure QuESTConfig.cmake (from template)
configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/QuESTConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfig.cmake"
        INSTALL_DESTINATION "${QuEST_INSTALL_CONFIGDIR}"
)

# Install them
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/QuESTConfigVersion.cmake"
        DESTINATION "${QuEST_INSTALL_CONFIGDIR}"
)

install(EXPORT QuESTTargets
        FILE QuESTTargets.cmake
        NAMESPACE QuEST::
        DESTINATION "${QuEST_INSTALL_CONFIGDIR}"
)
