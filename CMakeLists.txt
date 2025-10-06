# @author Oliver Thomson Brown
# @author Erich Essmann (patches including MSVC support)
# @author Tyson Jones (tidying + patches including clang multithreading)
# @author Luc Jaulmes (NUMA awareness, patching install)
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



# ============================
# Project
# ============================


cmake_minimum_required(VERSION 3.21)


project(QuEST
  VERSION 4.1.0
  DESCRIPTION "Quantum Exact Simulation Toolkit"
  LANGUAGES CXX C
)


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})



# ============================
# Dependencies
# ============================


# GNUInstallDirs to provide sensible default install directory names
cmake_path(SET ORG_INSTALL_PATH NORMALIZE "${CMAKE_INSTALL_PREFIX}")
cmake_path(APPEND CMAKE_INSTALL_PREFIX "quest")
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


# Maths
if (NOT WIN32)
  find_library(MATH_LIBRARY m REQUIRED)
endif()



# ============================
# Declare options
# ============================


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

if (VERBOSE_LIB_NAME)
  # Same headers will be used for several verbosely-named libraries
  set(MULTI_LIB_HEADERS 1)
  function(compile_option VAR VALUE)
    target_compile_definitions(QuEST PUBLIC ${VAR}=${VALUE})
  endfunction()
else()
  # Headers will be used for a single library with a single valid configuration
  set(MULTI_LIB_HEADERS 0)
  function(compile_option VAR VALUE)
    target_compile_definitions(QuEST PRIVATE ${VAR}=${VALUE})
    set(${VAR} ${VALUE} PARENT_SCOPE)
  endfunction()
endif()


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

# Throw on disallowed GPU combinations
if (ENABLE_CUDA AND ENABLE_HIP)
  message(FATAL_ERROR "QuEST cannot support CUDA and HIP simultaneously.")
endif()

if ((ENABLE_CUDA OR ENABLE_HIP) AND FLOAT_PRECISION STREQUAL 4)
  message(FATAL_ERROR "Quad precision is not supported on GPU. Please disable GPU acceleration or lower precision.")
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



# ============================
# Library
# ============================


add_library(QuEST)

# Add namespaced alias to support inclusion of QuEST as a subproject
add_library(QuEST::QuEST ALIAS QuEST)


# Set include directories
target_include_directories(QuEST
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quest/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set_target_properties(QuEST PROPERTIES
        VERSION     ${PROJECT_VERSION}
        SOVERSION   ${PROJECT_VERSION_MAJOR}
)


# Add required C and C++ standards.
# Note the QuEST interface(s) require only C11 and C++14, 
# while the source code is entirely C++ and requires C++17,
# and the tests further require C++20 (handled in tests/).
# Yet, we here specify C++17 for the source, and C11 as only
# applies to the C interface when users specify USER_SOURCE,
# to attemptedly minimise user confusion. Users wishing to
# link QuEST with C++14 should separate compilation.
target_compile_features(QuEST
  PUBLIC
  c_std_11
  cxx_std_17
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



# ============================
# Pass options to library
# ============================


compile_option(FLOAT_PRECISION ${FLOAT_PRECISION})


# OpenMP
if (ENABLE_MULTITHREADING)

  # find OpenMP, but fail gracefully...
  find_package(OpenMP QUIET)

  # so that we can customise the error message on MacOS
  if (NOT OpenMP_FOUND)
    set(ErrorMsg "Could not find OpenMP, necessary for enabling multithreading.")
    if (APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      string(APPEND ErrorMsg " Try first calling \n\tbrew install libomp\nthen\n\texport OpenMP_ROOT=$(brew --prefix)/opt/libomp")
    endif()
    message(FATAL_ERROR ${ErrorMsg})
  endif()

  compile_option(COMPILE_OPENMP 1)
  target_link_libraries(QuEST
    PRIVATE
    OpenMP::OpenMP_CXX
    OpenMP::OpenMP_C
  )

  # Find NUMA - location of NUMA headers
  if (WIN32)
    compile_option(NUMA_AWARE 0)
    message(WARNING "Building on Windows, QuEST will not be aware of numa locality")
  else()
    include(FindPkgConfig)
    pkg_search_module(NUMA numa IMPORTED_TARGET GLOBAL)
    if (${NUMA_FOUND})
      compile_option(NUMA_AWARE ${NUMA_FOUND})
      target_link_libraries(QuEST PRIVATE PkgConfig::NUMA)
      message(STATUS "NUMA awareness is enabled.")
    else()
      compile_option(NUMA_AWARE 0)
      message(WARNING "libnuma not found, QuEST will not be aware of numa locality")
    endif()
  endif()

  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+mt")
  endif()

else()

  # suppress GCC "unknown pragma" warning when OpenMP disabled
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(QuEST PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>)
  endif()
  
  compile_option(COMPILE_OPENMP 0)
endif()


# MPI
if (ENABLE_DISTRIBUTION)
  find_package(MPI REQUIRED
    COMPONENTS CXX
  )
  compile_option(COMPILE_MPI 1)
  target_link_libraries(QuEST
    PRIVATE
    MPI::MPI_CXX
  )
  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+mpi")
  endif()
else()
  compile_option(COMPILE_MPI 0)
endif()


# CUDA
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

  # beware that compile_option(COMPILE_CUDA) is deferred to below because
  # it is triggered by both/either ENABLE_CUDA and ENABLE_HIP

endif()


# cuQuantum
if (ENABLE_CUQUANTUM)
  find_package(CUQUANTUM REQUIRED)
  compile_option(COMPILE_CUQUANTUM 1)
  target_link_libraries(QuEST PRIVATE CUQUANTUM::cuStateVec)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+cuquantum")
  endif()
else()
  compile_option(COMPILE_CUQUANTUM 0)
endif()


# HIP
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

  compile_option(COMPILE_CUQUANTUM 0)
  target_link_libraries(QuEST PRIVATE hip::host)

  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+hip")
  endif()
endif()


# set COMPILE_CUDA
if (ENABLE_CUDA OR ENABLE_HIP)
  compile_option(COMPILE_CUDA 1)
else()
  compile_option(COMPILE_CUDA 0)
endif()


# v3 API
if (ENABLE_DEPRECATED_API)
  target_compile_definitions(QuEST PRIVATE INCLUDE_DEPRECATED_FUNCTIONS=1)

  if (VERBOSE_LIB_NAME)
    string(CONCAT LIB_NAME ${LIB_NAME} "+depr")
  endif()

else()
  target_compile_definitions(QuEST PRIVATE INCLUDE_DEPRECATED_FUNCTIONS=0)
endif()



# ============================
# Patch CPU performance
# ============================


# Patch performance of CPU std::complex arithmetic operator overloads.
# The cpu_subroutines.cpp file makes extensive use of std::complex operator
# overloads, and alas these are significantly slower than hand-rolled 
# arithmetic, due to their NaN and inf checks, and interference with SIMD.
# It is crucial to pass additional optimisation flags to this file to restore
# hand-rolled performance (else QuEST v3 is faster than v4 eep). In theory,
# we can achieve this with specific, relatively 'safe' flags such as LLVM's:
#     -ffinite-math-only -fno-signed-zeros -ffp-contract=fast
# However, it is a nuisance to find equivalent flags for different compilers
# and monitor their performance vs accuracy trade-offs. So instead, we use the
# much more aggressive and ubiquitous -Ofast flag to guarantee performance. 
# This introduces many potentially dangerous optimisations, such as asserting
# associativity of flops, which would break techniques like Kahan summation.
# The cpu_subroutines.cpp must ergo be very conscious of these optimisations.
# We here also explicitly inform the file cpu_subroutines.cpp whether or not
# we are passing the flags, so it can detect/error when flags are forgotten.

if (CMAKE_BUILD_TYPE STREQUAL "Release")

  # Release build will pass -Ofast when known for the given compiler, and
  # fallback to giving a performance warning and proceeding with compilation

  if (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|Cray|CrayClang|GNU|HP|Intel|IntelLLVM|NVHPC|NVIDIA|XL|XLClang")
    set(patch_flags "-Ofast")
    set(patch_macro "-DCOMPLEX_OVERLOADS_PATCHED=1")
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "HP")
    set(patch_flags "+Ofast")
    set(patch_macro "-DCOMPLEX_OVERLOADS_PATCHED=1")
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(patch_flags "/fp:fast")
    set(patch_macro "-DCOMPLEX_OVERLOADS_PATCHED=1")
  else()
    message(WARNING 
      "The compiler (${CMAKE_CXX_COMPILER_ID}) is unrecognised and so crucial optimisation flags have not been "
      "passed to the CPU backend. These flags are necessary for full performance when performing complex algebra, "
      "otherwise a slowdown of 3-50x may be observed. Please edit the root CMakeLists.txt to include flags which are "
      "equivalent to GNU's -Ofast flag for your compiler (search this warning), or contact the QuEST developers for help."
    )
    set(patch_flags "")
    set(patch_macro "-DCOMPLEX_OVERLOADS_PATCHED=0")
  endif()
  
else()

  # Non-release builds (e.g. Debug) will pass no optimisation flags, and will
  # communicate to cpu_subroutines.cpp that this is intentional via a macro

  set(patch_flags "")
  set(patch_macro "-DCOMPLEX_OVERLOADS_PATCHED=0")

endif()

set_source_files_properties(
  quest/src/cpu/cpu_subroutines.cpp
  PROPERTIES
  COMPILE_FLAGS "${patch_flags} ${patch_macro}"
)



# ============================
# Pass files to library
# ============================


# add math library
if (NOT MSVC)
  target_link_libraries(QuEST PRIVATE ${MATH_LIBRARY})
endif()


# Set output name
set_target_properties(QuEST PROPERTIES OUTPUT_NAME ${LIB_NAME})


# Add source files
add_subdirectory(quest)



# ============================
# Examples
# ============================


# min example is always built
add_executable(min_example
  examples/tutorials/min_example.c
)
target_link_libraries(min_example PRIVATE QuEST::QuEST)

install(TARGETS min_example
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)


# all examples optionally built
if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()


## RATH
set(BUILD_RPATH_USE_ORIGIN ON)
if(APPLE)
  set(_RPATH_ORIGIN "@loader_path")
else()
  set(_RPATH_ORIGIN "$ORIGIN")
endif()

set(_INSTALL_RPATH  "${_RPATH_ORIGIN}/../${CMAKE_INSTALL_LIBDIR}")
set(_BUILD_RPATH    "${_RPATH_ORIGIN};$<TARGET_FILE_DIR:QuEST>")

# A tiny helper function so you can call it for every target
function(setup_quest_rpath tgt)
  set_target_properties(${tgt} PROPERTIES
          BUILD_RPATH   "${_BUILD_RPATH}"
          INSTALL_RPATH "${_INSTALL_RPATH}"
          # keeps RPATH from being stripped when installing
          INSTALL_RPATH_USE_LINK_PATH TRUE
  )
endfunction()

setup_quest_rpath(QuEST)
setup_quest_rpath(min_example)



# ============================
# User source
# ============================


# validate
if (USER_SOURCE AND NOT OUTPUT_EXE)
    message(SEND_ERROR "USER_SOURCE specified, but not OUTPUT_EXE.")
endif()
if (OUTPUT_EXE AND NOT USER_SOURCE)
    message(SEND_ERROR "OUTPUT_EXE specified, but not USER_SOURCE.")
endif()


# compile user source
if (USER_SOURCE AND OUTPUT_EXE)
  message(STATUS "Compiling ${USER_SOURCE} to executable ${OUTPUT_EXE}.")

  add_executable(${OUTPUT_EXE} ${USER_SOURCE})
  target_link_libraries(${OUTPUT_EXE} PUBLIC QuEST)
  install(TARGETS ${OUTPUT_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  setup_quest_rpath(${OUTPUT_EXE})
endif()



# ============================
# Tests
# ============================


if (ENABLE_TESTING)

  # try find Catch2
  set(CatchVersion 3.8.0)
  find_package(Catch2 ${CatchVersion} QUIET)

  # else try download Catch2
  if (NOT TARGET Catch2::Catch2 AND DOWNLOAD_CATCH2)
    message(STATUS "Catch2 not found, it will be downloaded and built in the build directory.")
    Include(FetchContent)

    FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v${CatchVersion}
    )
    
    FetchContent_MakeAvailable(Catch2)
  
  # otherwise fail
  else()
    # We won't magically find it here, but this is the easiest way to
    # a) Force the build to fail, and
    # b) Print out all the useful information for helping CMake find Catch2
    find_package(Catch2 ${CatchVersion} REQUIRED)
  endif()

  # compile tests
  include(Catch)
  enable_testing()
  add_subdirectory(tests)
endif()



# ============================
# Installation
# ============================


install(TARGETS QuEST
  EXPORT QuESTTargets
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
)


# Write CMake version file for QuEST
set(QuEST_INSTALL_CONFIGDIR "${CMAKE_INSTALL_LIBDIR}/cmake/QuEST")


# Write QuESTConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)


# Configure QuESTConfig.cmake (from template)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/QuESTConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
  INSTALL_DESTINATION "${QuEST_INSTALL_CONFIGDIR}"
)


# Install them
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}Config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}ConfigVersion.cmake"
  DESTINATION "${QuEST_INSTALL_CONFIGDIR}"
)

install(FILES 
  "${CMAKE_CURRENT_SOURCE_DIR}/quest/include/quest.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/include/quest/include/config.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/quest/include"
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/quest/include"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/quest"
  FILES_MATCHING PATTERN "*.h"
  PATTERN "quest.h" EXCLUDE
)

install(
  EXPORT QuESTTargets
  FILE "${LIB_NAME}Targets.cmake"
  NAMESPACE QuEST::
  DESTINATION "${QuEST_INSTALL_CONFIGDIR}"
)
