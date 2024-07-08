# author : Ania Brown
# author : Jacob Wilkins
# author : Balint Koczor (Windows compatibility)
# author : Tyson Jones (testing)
# author : Oliver Thomson Brown (v4)

cmake_minimum_required(VERSION 3.13)

project(QuEST
  VERSION 4.0.0
  DESCRIPTION "Quantum Exact Simulation Toolkit"
  LANGUAGES CXX C
)

## Dependencies

# Math
find_library(MATH_LIBRARY m REQUIRED)

## Configuration options

# Build type
# Default to "Release"
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

# Precision
set(PRECISION 2 
  CACHE 
  STRING 
  "Whether to use single, double, or quad floating point precision in the state vector. {1,2,4}"
)
set_property(CACHE PRECISION PROPERTY STRINGS
  1
  2
  4
)
message(STATUS "Precision set to ${PRECISION}.")

# Testing
option(
  TESTING
  "whether the test suite will be built alongside the QuEST library. Turned ON by default."
  ON
)
message(STATUS "Testing is turned ${TESTING}. Set TESTING to modify.")

# Multithreading
option(
  MULTITHREADED 
  "Whether QuEST will be built with shared-memory parallelism support using OpenMP. Turned ON by default." 
  ON
)
message(STATUS "Multithreading is turned ${MULTITHREADED}. Set MULTITHREADED to modify.")

# Distribution
option(
  DISTRIBUTED 
  "Whether QuEST will be built with distributed parallelism support using MPI. Turned OFF by default." 
  OFF
)
message(STATUS "Distribution is turned ${DISTRIBUTED}. Set DISTRIBUTED to modify.")

# GPU Acceleration
option(
  GPUACCELERATED 
  "Whether QuEST will be built with support for GPU acceleration. Turned OFF by default."
  OFF
)
message(STATUS "GPU Acceleration is turned ${GPUACCELERATED}. Set GPUACCELERATED to modify.")

# CuQuantum
option(
  CUQUANTUM
  "Whether QuEST will be built with support for NVIDIA CuQuantum. Turned OFF by default."
  OFF
)
message(STATUS "CuQuantum support is turned ${CUQUANTUM}. Set CUQUANTUM to modify.")

## Library

add_library(QuEST)

# Add required C and C++ standards
target_compile_features(QuEST
  PUBLIC
  c_std_11
  cxx_std_14
)

# Turn on all compiler warnings
if (MSVC)
  target_compile_options(QuEST PRIVATE /W4)
else()
  target_compile_options(QuEST PRIVATE -Wall)
endif()

# Set user options
target_compile_definitions(QuEST PUBLIC FLOAT_PRECISION=${PRECISION})

if (MULTITHREADED)
  find_package(OpenMP REQUIRED)
  target_compile_definitions(QuEST PUBLIC ENABLE_MULTITHREADING=1)
  target_link_libraries(QuEST
    PUBLIC
    OpenMP::OpenMP_CXX
    OpenMP::OpenMP_C
  )
else()
  target_compile_definitions(QuEST PUBLIC ENABLE_MULTITHREADING=0)
endif()

if (DISTRIBUTED)
else()
  target_compile_definitions(QuEST PUBLIC ENABLE_DISTRIBUTION=0)
endif()

if (GPUACCELERATED)
else()
  target_compile_definitions(QuEST PUBLIC ENABLE_GPU_ACCELERATION=0)
endif()

if (CUQUANTUM)
else()
  target_compile_definitions(QuEST PUBLIC ENABLE_CUQUANTUM=0)
endif()

# add math library
if (NOT MSVC)
  target_link_libraries(QuEST PUBLIC ${MATH_LIBRARY})
endif()

target_include_directories(QuEST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(quest)

## Examples

#add_subdirectory(examples)

## Tests
if (TESTING)
  add_subdirectory(tests)
endif()