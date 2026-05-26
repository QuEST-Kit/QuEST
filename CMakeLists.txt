# @author Oliver Thomson Brown
# @author Erich Essmann (patches including MSVC support)
# @author Tyson Jones (centralising macros, doc, patches including clang multithreading)
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
  VERSION 4.2.0
  DESCRIPTION "Quantum Exact Simulation Toolkit"
  LANGUAGES CXX C
)


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})



# ============================
# Dependencies
# ============================


# GNUInstallDirs to provide sensible default install directory names
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)


# Maths
if (NOT WIN32)
  find_library(MATH_LIBRARY m REQUIRED)
endif()



# ============================
# Obtain options
# ============================


# Build type
# Default to "Release"
# Using recipe from Kitware Blog post
# https://www.kitware.com/cmake-and-the-default-build-type/
set(quest_default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${quest_default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${quest_default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


# Library type
# Shared library by default when we are the top level cmake project
if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_SHARED_LIBS "Build shared library. Turned ON by default." ON)
  message(STATUS "Shared library is turned ${BUILD_SHARED_LIBS}. Set BUILD_SHARED_LIBS to modify.")
endif ()

# Library naming
set(QUEST_OUTPUT_LIB_NAME QuEST
  CACHE
  STRING
  "Change library name. QUEST_OUTPUT_LIB_NAME is QuEST by default."
)
message(STATUS "Library will be named lib${QUEST_OUTPUT_LIB_NAME}. Set QUEST_OUTPUT_LIB_NAME to modify.")

option(QUEST_APPEND_CONFIG_TO_LIB_NAME "Modify library name based on compilation configuration. Turned OFF by default." OFF)
message(STATUS "Verbose library naming is turned ${QUEST_APPEND_CONFIG_TO_LIB_NAME}. Set QUEST_APPEND_CONFIG_TO_LIB_NAME to modify.")


# Precision
set(QUEST_FLOAT_PRECISION 2
  CACHE
  STRING
  "Whether to use single, double, or quad floating point precision in the state vector. {1,2,4}"
)
set_property(CACHE QUEST_FLOAT_PRECISION PROPERTY STRINGS
  1
  2
  4
)
message(STATUS "Precision set to ${QUEST_FLOAT_PRECISION}. Set QUEST_FLOAT_PRECISION to modify.")


# Examples
option(
  QUEST_BUILD_EXAMPLES
  "Whether the example programs will be built alongside the QuEST library. Turned OFF by default."
  OFF
)
message(STATUS "Examples are turned ${QUEST_BUILD_EXAMPLES}. Set QUEST_BUILD_EXAMPLES to modify.")


# Testing
option(
  QUEST_BUILD_TESTS
  "Whether the test suite will be built alongside the QuEST library. Turned ON by default."
  OFF
)
message(STATUS "Testing is turned ${QUEST_BUILD_TESTS}. Set QUEST_BUILD_TESTS to modify.")

option(
  QUEST_TESTS_DOWNLOAD_CATCH2
  "Whether Catch2 v3 will be downloaded if it is not found. Turned ON by default."
  ON
)


# Multithreading
option(
  QUEST_ENABLE_OMP
  "Whether QuEST will be built with shared-memory parallelism support using OpenMP. Turned ON by default."
  ON
)
message(STATUS "Multithreading is turned ${QUEST_ENABLE_OMP}. Set QUEST_ENABLE_OMP to modify.")


# NUMA
option(
  QUEST_ENABLE_NUMA
  "Whether QuEST will be built with NUMA awareness, when also using OpenMP. Turned ON by default."
  ON
)
message(STATUS "NUMA awareness is turned ${QUEST_ENABLE_NUMA}. Set QUEST_ENABLE_NUMA to modify.")


# Distribution
option(
  QUEST_ENABLE_MPI
  "Whether QuEST will be built with distributed parallelism support using MPI. Turned OFF by default."
  OFF
)
message(STATUS "Distribution is turned ${QUEST_ENABLE_MPI}. Set QUEST_ENABLE_MPI to modify.")


# GPU Acceleration
option(
  QUEST_ENABLE_CUDA
  "Whether QuEST will be built with support for NVIDIA GPU acceleration. Turned OFF by default."
  OFF
)
message(STATUS "NVIDIA GPU acceleration is turned ${QUEST_ENABLE_CUDA}. Set QUEST_ENABLE_CUDA to modify.")

option(
  QUEST_ENABLE_CUQUANTUM
  "Whether QuEST will be built with support for NVIDIA cuQuantum. Turned OFF by default."
  OFF
)
message(STATUS "CuQuantum support is turned ${QUEST_ENABLE_CUQUANTUM}. Set QUEST_ENABLE_CUQUANTUM to modify.")

option(
  QUEST_ENABLE_HIP
  "Whether QuEST will be built with support for AMD GPU acceleration. Turned OFF by default."
  OFF
)
message(STATUS "AMD GPU acceleration is turned ${QUEST_ENABLE_HIP}. Set QUEST_ENABLE_HIP to modify.")


# Deprecated API
option(
  QUEST_ENABLE_DEPRECATED_API
  "Whether QuEST will be built with deprecated API support. Turned OFF by default."
  OFF
)
message(STATUS "Deprecated API support is turned ${QUEST_ENABLE_DEPRECATED_API}. Set QUEST_ENABLE_DEPRECATED_API to modify.")

option(
  QUEST_DISABLE_DEPRECATION_WARNINGS
  "Whether to disable compile-time warnings ordinarily triggered by use of the deprecated API. Turned OFF by default."
  OFF
)
message(STATUS "Disabling of deprecated API warnings is turned ${QUEST_DISABLE_DEPRECATION_WARNINGS}. Set QUEST_DISABLE_DEPRECATION_WARNINGS to modify.")

option(QUEST_INSTALL_BINARIES "Whether to include example and user binaries in the install." OFF)



# ============================
# Validate options
# ============================


if (QUEST_ENABLE_CUDA AND QUEST_ENABLE_HIP)
  message(FATAL_ERROR "QuEST cannot support CUDA and HIP simultaneously.")
endif()


if ((QUEST_ENABLE_CUDA OR QUEST_ENABLE_HIP) AND QUEST_FLOAT_PRECISION STREQUAL 4)
  message(FATAL_ERROR "Quad precision is not supported on GPU. Please disable GPU acceleration or lower precision.")
endif()


if (QUEST_ENABLE_CUQUANTUM AND NOT QUEST_ENABLE_CUDA)
  message(FATAL_ERROR "Use of cuQuantum requires CUDA.")
endif()


if(WIN32)
  
  # Force MSVC to export all symbols in a shared library, like GCC and clang
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

  if (QUEST_BUILD_TESTS AND BUILD_SHARED_LIBS)
    message(WARNING "Compiling the tests on Windows requires BUILD_SHARED_LIBS=OFF which we now force.")
    set(BUILD_SHARED_LIBS OFF)
  endif()

  if (QUEST_ENABLE_DEPRECATED_API)
    message(FATAL_ERROR "The deprecated API is not compatible with MSVC.")
  endif()

endif()


# Encourage high-performance Release build

# Taken from Kitware's exmaple of problematic code at
# https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#build-configurations
# We have to use it here as we want to message the user! We are handling multiconfig
# separately anyway.
string(TOLOWER "${CMAKE_BUILD_TYPE}" build_type)

if(CMAKE_CONFIGURATION_TYPES)
  message(WARNING
    "You are using a multi-config generator, so make sure to subsequently "
    "build with '--config Release' for best performance.")
elseif(NOT (build_type STREQUAL "release" OR build_type STREQUAL "relwithdebinfo"))
  message(WARNING
    "Using a non-release build (CMAKE_BUILD_TYPE='${CMAKE_BUILD_TYPE}') "
    "can significantly degrade performance. "
    "Consider using -DCMAKE_BUILD_TYPE=Release.")
endif()



# ============================
# Extend verbose library name
# ============================


if (QUEST_APPEND_CONFIG_TO_LIB_NAME)

  string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "-fp${QUEST_FLOAT_PRECISION}")

  if (QUEST_ENABLE_OMP)
    string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "+mt")
  endif()

  if (QUEST_ENABLE_MPI)
    string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "+mpi")
  endif()

  if (QUEST_ENABLE_CUDA)
    string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "+cuda")
  endif()

  if (QUEST_ENABLE_HIP)
    string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "+hip")
  endif()

  if (QUEST_ENABLE_CUQUANTUM)
    string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "+cuquantum")
  endif()

  if (QUEST_ENABLE_DEPRECATED_API)
    string(CONCAT QUEST_OUTPUT_LIB_NAME ${QUEST_OUTPUT_LIB_NAME} "+depr")
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
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
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
# applies to the C interface when users specify USER_SOURCE_NAMES,
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
# Link optional dependencies
# ============================


# OpenMP
if (QUEST_ENABLE_OMP)

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

  target_link_libraries(QuEST
    PRIVATE
    OpenMP::OpenMP_CXX
    OpenMP::OpenMP_C
  )

else()

  # suppress GCC "unknown pragma" warning when OpenMP disabled
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(QuEST PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wno-unknown-pragmas>)
  endif()
  
endif()


# NUMA (only relevant when multithreading)
if (QUEST_ENABLE_OMP AND QUEST_ENABLE_NUMA)

  # Find NUMA - location of NUMA headers
  if (WIN32)
    set(QUEST_ENABLE_NUMA 0)
    message(WARNING "Building on Windows, QuEST will not be aware of numa locality")
  else()
    include(FindPkgConfig)
    pkg_search_module(NUMA numa IMPORTED_TARGET GLOBAL)
    if (${NUMA_FOUND})
      set(QUEST_ENABLE_NUMA ${NUMA_FOUND})
      target_link_libraries(QuEST PRIVATE PkgConfig::NUMA)
      message(STATUS "NUMA awareness is enabled.")
    else()
      set(QUEST_ENABLE_NUMA 0)
      message(WARNING "libnuma not found, QuEST will not be aware of numa locality")
    endif()
  endif()

else()
  set(QUEST_ENABLE_NUMA 0)
endif()


# MPI
if (QUEST_ENABLE_MPI)
  find_package(MPI REQUIRED
    COMPONENTS CXX
  )
  target_link_libraries(QuEST
    PRIVATE
    MPI::MPI_CXX
  )
endif()


# CUDA
if (QUEST_ENABLE_CUDA)

  # make nvcc use user cxx-compiler as default host (before cuda-host is set below)
  if (NOT DEFINED CMAKE_CUDA_HOST_COMPILER)
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  endif()

  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  
  set_property(TARGET QuEST PROPERTY CUDA_STANDARD 20)

  # force MSVC to use the modern preprocessor
  if (MSVC)
    target_compile_options(QuEST PRIVATE
      $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    )
  endif()

endif()


# HIP
if (QUEST_ENABLE_HIP)

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

  target_link_libraries(QuEST PRIVATE hip::host)

endif()


# cuQuantum
if (QUEST_ENABLE_CUQUANTUM)
  find_package(CUQUANTUM REQUIRED)
  target_link_libraries(QuEST PRIVATE CUQUANTUM::cuStateVec)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()



# ===============================
# Set options to save in config.h
# ===============================


# set vars which will be written to config.h.in (auto-converted to 0 or 1)
set(QUEST_COMPILE_OPENMP ${QUEST_ENABLE_OMP})
set(QUEST_COMPILE_MPI ${QUEST_ENABLE_MPI})
set(QUEST_COMPILE_CUQUANTUM ${QUEST_ENABLE_CUQUANTUM})
set(QUEST_INCLUDE_DEPRECATED_FUNCTIONS ${QUEST_ENABLE_DEPRECATED_API})


# (for the love of God cmake, create a concise syntax for this)
if (QUEST_ENABLE_CUDA OR QUEST_ENABLE_HIP)
  set(QUEST_COMPILE_CUDA 1)
else()
  set(QUEST_COMPILE_CUDA 0)
endif()


# these vars are already set, but repeated here for clarity
set(QUEST_FLOAT_PRECISION ${QUEST_FLOAT_PRECISION})
set(QUEST_ENABLE_NUMA ${QUEST_ENABLE_NUMA})
set(QUEST_DISABLE_DEPRECATION_WARNINGS ${QUEST_DISABLE_DEPRECATION_WARNINGS})


# these do not appear in src but are saved for record-keeping in config.h.in
set(QUEST_COMPILE_HIP ${QUEST_ENABLE_HIP})



# ============================
# Pass files to library
# ============================


# add math library
if (NOT MSVC)
  target_link_libraries(QuEST PRIVATE ${MATH_LIBRARY})
endif()


# Set output name
set_target_properties(QuEST PROPERTIES OUTPUT_NAME ${QUEST_OUTPUT_LIB_NAME})


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

if (QUEST_INSTALL_BINARIES)
  install(TARGETS min_example
    RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif ()


# all examples optionally built
if (QUEST_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()


## RPATH
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
if (USER_SOURCE_NAMES AND NOT USER_OUTPUT_EXE_NAME)
    message(SEND_ERROR "USER_SOURCE_NAMES specified, but not USER_OUTPUT_EXE_NAME.")
endif()
if (USER_OUTPUT_EXE_NAME AND NOT USER_SOURCE_NAMES)
    message(SEND_ERROR "USER_OUTPUT_EXE_NAME specified, but not USER_SOURCE_NAMES.")
endif()


# compile user source
if (USER_SOURCE_NAMES AND USER_OUTPUT_EXE_NAME)
  message(STATUS "Compiling ${USER_SOURCE_NAMES} to executable ${USER_OUTPUT_EXE_NAME}.")

  add_executable(${USER_OUTPUT_EXE_NAME} ${USER_SOURCE_NAMES})
  target_link_libraries(${USER_OUTPUT_EXE_NAME} PUBLIC QuEST)
   
  if (QUEST_INSTALL_BINARIES)
    install(TARGETS ${USER_OUTPUT_EXE_NAME} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
  
  setup_quest_rpath(${USER_OUTPUT_EXE_NAME})
endif()



# ============================
# Tests
# ============================


if (QUEST_BUILD_TESTS)

  # try find Catch2
  set(CatchVersion 3.8.0)
  find_package(Catch2 ${CatchVersion} QUIET)

  # else try download Catch2
  if (NOT TARGET Catch2::Catch2 AND QUEST_TESTS_DOWNLOAD_CATCH2)
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
set(quest_install_config_dir "${CMAKE_INSTALL_LIBDIR}/cmake/QuEST")


# Write QuESTConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${QUEST_OUTPUT_LIB_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)


# Configure QuESTConfig.cmake (from template)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/QuESTConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${QUEST_OUTPUT_LIB_NAME}Config.cmake"
  INSTALL_DESTINATION "${quest_install_config_dir}"
)


# Install them
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${QUEST_OUTPUT_LIB_NAME}Config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${QUEST_OUTPUT_LIB_NAME}ConfigVersion.cmake"
  DESTINATION "${quest_install_config_dir}"
)

install(FILES 
  "${CMAKE_CURRENT_SOURCE_DIR}/quest/include/quest.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/quest/include/config.h"
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
  FILE "${QUEST_OUTPUT_LIB_NAME}Targets.cmake"
  NAMESPACE QuEST::
  DESTINATION "${quest_install_config_dir}"
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif ()

