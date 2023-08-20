# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----- QuEST LIBRARY BUILD SYSTEM --------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Builds QuEST as a shared library libQuEST.so


# -----------------------------------------------------------------------------
# ----- LIBRARY USER SETTINGS -------------------------------------------------
# -----------------------------------------------------------------------------

# These variables are cached across multiple calls to cmake. Use cmake -LH for 
# current cache values

# To override, use:
#   cmake -D[VAR]=[VALUE] [Path to root CMakeLists.txt]
# eg
#   cmake -DDISTRIBUTED=1 ..

cmake_minimum_required(VERSION 3.7)
project(QuEST_Library)

option(MULTITHREADED "Switches on OpenMP Threading if an OpenMP implementation can be found. Set to 0 to disable." 1)

option(DISTRIBUTED "Whether the program will be run in distributed mode using MPI. Set to 1 to enable" 0)

set(PRECISION 2 CACHE STRING
    "Whether to use single, double or quad floating point precision in the state vector. {1,2,4}")

option(GPUACCELERATED "Whether to program will run on GPU. Set to 1 to enable" 0)

set(GPU_COMPUTE_CAPABILITY 30 CACHE STRING "GPU hardware dependent, lookup at https://developer.nvidia.com/cuda-gpus. Write without fullstop")

option(USE_HIP "Whether to use HIP for GPU code compilation for AMD GPUs. Set to 1 to enable" 0)

set(GPU_ARCH gfx90 CACHE STRING "GPU hardware dependent, used for AMD GPUs when USE_HIP=1. Lookup at https://llvm.org/docs/AMDGPUUsage.html#amdgpu-processor-table. Write without fullstop")

option(USE_CUQUANTUM "Whether to use NVIDIA's cuQuantum library (requires prior installation) in lieu of QuEST's bespoke GPU kernel. Set to 1 to enable." 0)


# *****************************************************************************
# ***** NO CHANGES SHOULD BE REQUIRED FROM THE USER BEYOND THIS POINT *********
# *****************************************************************************

# Show the user their settings
message(STATUS "Precision is ${PRECISION}")
message(STATUS "GPU acceleration is ${GPUACCELERATED}")
message(STATUS "OMP acceleration is ${MULTITHREADED}")
message(STATUS "MPI distribution is ${DISTRIBUTED}")
if (${GPUACCELERATED})
  message(STATUS "HIP compilation is ${USE_HIP}")
  message(STATUS "cuQuantum compilation is ${USE_CUQUANTUM}")
endif()


# -----------------------------------------------------------------------------
# ----- CHECK FOR VALID CONFIG SETTINGS ---------------------------------------
# -----------------------------------------------------------------------------

# ----- FATAL -----------------------------------------------------------------

if (${DISTRIBUTED} AND ${GPUACCELERATED})
    message(FATAL_ERROR "DISTRIBUTED=${DISTRIBUTED} and \
        GPUACCELERATED=${GPUACCELERATED} set but \
        distributed GPU acceleration not supported. Aborting")
endif()

if ( NOT(${PRECISION} EQUAL 1) AND
     NOT(${PRECISION} EQUAL 2) AND 
     NOT(${PRECISION} EQUAL 4) )
    message(FATAL_ERROR "PRECISION=${PRECISION} set but valid options for \
        precision are only 1 for single precision, 2 for double precision \
        or 4 for quad precision. Aborting")
endif()

if ( (${PRECISION} EQUAL 4) AND 
        ${GPUACCELERATED} )
    message(FATAL_ERROR "PRECISION=${PRECISION} but quad precision is not \
        supported on GPU. Aborting")
endif()

# ----- WARNINGS --------------------------------------------------------------

if (${GPUACCELERATED} AND ${MULTITHREADED})
    message(WARNING "MULTITHREADED=${MULTITHREADED} and \
        GPUACCELERATED=${GPUACCELERATED} set but GPU \
        version makes no use of multithreading. Ignoring multithreading settings")
endif()

#TODO Add other supported Clang versions if found
if ( ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") AND 
        ${GPUACCELERATED} AND
        NOT("${CMAKE_C_COMPILER_VERSION}" STREQUAL "3.7.0") )
    message(WARNING "Some versions of Clang are not NVIDIA-GPU compatible. \
        If compilation fails, try Clang 3.7.")
endif()

if ( ${GPUACCELERATED} AND 
    ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") AND
    ("${CMAKE_SYSTEM_NAME}" STREQUAL "DARWIN") ) # DARWIN means Mac OS X
    message(WARNING "On some platforms (e.g. OSX), NVIDIA-GPUs are not \
        compatible with GNU compilers. If compilation fails, try an \
        alternative compiler, like Clang 3.7")
endif()


# -----------------------------------------------------------------------------
# ----- FIND PACKAGES ---------------------------------------------------------
# -----------------------------------------------------------------------------

if (DISTRIBUTED)
  find_package(MPI REQUIRED)
  # Backwards compatibility
  if (DEFINED MPIEXEC)
    set (MPIEXEC_EXECUTABLE ${MPIEXEC} CACHE STRING "MPI Executable")
  endif()
  include_directories(${MPI_INCLUDE_PATH})
endif()

if (GPUACCELERATED)
  if (USE_HIP)

    if(NOT DEFINED HIP_PATH)
      if(NOT DEFINED ENV{HIP_PATH})
        message(WARNING "WARNING: HIP_PATH is not defiend. Using default HIP_PATH=/opt/rocm/hip    " ${HIP_VERSION})
        set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
      else()
        set(HIP_PATH $ENV{HIP_PATH}	CACHE PATH "Path to which HIP has been installed")
      endif()
    endif()

    if(EXISTS "${HIP_PATH}")
      set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
      find_package(HIP REQUIRED)
      message(STATUS "Found HIP: " ${HIP_VERSION})
      message(STATUS "HIP PATH: " ${HIP_PATH})
    endif()
      
    ADD_DEFINITIONS( -DUSE_HIP )
    ADD_DEFINITIONS( -D__HIP_PLATFORM_AMD__ )

  elseif (USE_CUQUANTUM)
    find_package(CUDA REQUIRED)
    ADD_DEFINITIONS( -DUSE_CUQUANTUM )
  else()
    find_package(CUDA REQUIRED)
  endif()  
endif()



# -----------------------------------------------------------------------------
# ----- SET COMPILER FLAGS ----------------------------------------------------
# -----------------------------------------------------------------------------


# ----- OPENMP ----------------------------------------------------------------

if (${MULTITHREADED} AND NOT ${GPUACCELERATED})
  find_package(OpenMP)

  # If found, we must also check the version
  if (OPENMP_FOUND)
    # Version number is not cached, and we need to test for old versions
    set(OpenMP_C_VERSION ${OpenMP_C_VERSION} CACHE STRING "OpenMP C version")

    # MSVC, for instance, only implements OpenMP 2.0 (as of 2019)
    if (OpenMP_C_VERSION VERSION_LESS "2.0")  # todo find the real minimum required
      set(MULTITHREADED 0)
      message(WARNING "Found OpenMP ${OpenMP_C_VERSION} but this is too \
             old. Turning OpenMP support OFF.")
    else ()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif ()
  endif ()

endif ()

# ----- CUDA FLAGS ------------------------------------------------------------

if (GPUACCELERATED)
  if (USE_HIP)
    if( NOT DEFINED GPU_ARCH AND DEFINED GPU_COMPUTE_CAPABILITY )
      set ( GPU_ARCH ${GPU_COMPUTE_CAPABILITY})
    endif()

    set(HIP_HIPCC_FLAGS "${HIP_HIPCC_FLAGS} -fPIC  --offload-arch=${GPU_ARCH}" )
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
      set (HIP_HIPCC_FLAGS "${HIP_HIPCC_FLAGS} \
          -O2"
      )
    elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      set (HIP_HIPCC_FLAGS "${HIP_HIPCC_FLAGS} \
          -G -g -lineinfo"
      )
    else()
      set (HIP_HIPCC_FLAGS "${HIP_HIPCC_FLAGS} \
          -O2"
      )
    endif()
  else()
    set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
        -arch=compute_${GPU_COMPUTE_CAPABILITY} -code=sm_${GPU_COMPUTE_CAPABILITY}"
    )
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
            -O2"
        )
    elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        if ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
            set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
                -G"
            )
        else()
            set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
                -G -g -lineinfo"
            )
        endif()
    else()
        set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
            -O2"
        )
    endif()
  endif()  
endif()


# ----- C COMPILER FLAGS --------------------------------------------------

# set C flags that are common between compilers and build types
set (CMAKE_C_STANDARD 99)

# Use -O2 for all but debug mode by default 
if (NOT("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
        -O2"
    )
endif()

# Set c flags to use in debug mode

if (NOT("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC"))
  set(CMAKE_C_FLAGS_DEBUG 
      "-g"
  )
else()
  set(CMAKE_C_FLAGS_DEBUG 
      "-Zi"
  )
endif()

# Set c flags for release
set(CMAKE_C_FLAGS_RELEASE 
    "-O2"
)

# TODO standardize
# set C compiler flags based on compiler type
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -Wall"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -Wall"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable cpu-dispatch"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
  string(REGEX REPLACE "/W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  string(REGEX REPLACE "-W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -w"
  )
endif()

# ----- C++ COMPILER FLAGS --------------------------------------------------

# set C++ flags that are common between compilers and build types
if (USE_CUQUANTUM)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
else ()
  set (CMAKE_CXX_STANDARD 98)
endif ()

# Use -O2 for all but debug mode by default 
if (NOT("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} \
        -O2"
    )
endif()

# Set c++ flags for release
set(CMAKE_CXX_FLAGS_RELEASE
    "-O2"
)

# Set c++ flags to use in debug mode
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS_DEBUG 
    "-Zi -MDd"
  )
else()
  set(CMAKE_CXX_FLAGS_DEBUG 
    "-g"
  )
endif()

# set C++ compiler flags based on compiler type
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx -Wall"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx -Wall"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch -Wall"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
  string(REGEX REPLACE "/W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  string(REGEX REPLACE "-W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -w"
  )
endif()

if (VERBOSE_CMAKE)
    message("-- Compiler flags:")
    message("   C Compiler ID: ${CMAKE_C_COMPILER_ID}")
    message("   C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
    message("   C compilation flags: ${CMAKE_C_FLAGS}")
    message("   CUDA compilation flags: ${CUDA_NVCC_FLAGS}")
    message("   HIP compilation flags: ${HIP_HIPCC_FLAGS}")
    message("   C++ compilation flags: ${CMAKE_CXX_FLAGS}")

    message("-- Build type: ${CMAKE_BUILD_TYPE}")
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        message("   C debug flags: ${CMAKE_C_FLAGS_DEBUG}")
        message("   C++ debug flags: ${CMAKE_CXX_FLAGS_DEBUG}")
    endif()
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        message("   C release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
        message("   C++ release flags: ${CMAKE_CXX_FLAGS_RELEASE}")
    endif()
endif()

# -----------------------------------------------------------------------------
# ----- BUILD LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------


add_subdirectory(src)

if (NOT WIN32)
    set(BUILD_SHARED_LIBS TRUE CACHE BOOL "Whether to build a dynamic library")
else ()
    set(BUILD_SHARED_LIBS FALSE CACHE BOOL "Whether to build a dynamic library")
endif ()


if (GPUACCELERATED)
  if (USE_HIP)
    hip_add_library(QuEST
        ${QuEST_SRC} 
    )
  else()
    cuda_add_library(QuEST
        ${QuEST_SRC}
    )
  endif()  
else ()
    add_library(QuEST
        ${QuEST_SRC}
    )
endif()

# ----- Location of header files ----------------------------------------------

target_include_directories(QuEST 
    PRIVATE src
    PUBLIC include
)


# ----- Definitions -----------------------------------------------------------

target_compile_definitions(QuEST
    PUBLIC
    QuEST_PREC=${PRECISION}
)

# -----------------------------------------------------------------------------
# ----- LINK LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

# ----- OMP -------------------------------------------------------------------

if (MULTITHREADED AND OPENMP_FOUND)
  target_link_libraries(QuEST PUBLIC OpenMP::OpenMP_C)
endif ()

# ----- MPI -------------------------------------------------------------------

target_link_libraries(QuEST PUBLIC ${MPI_C_LIBRARIES})

# ----- GPU -------------------------------------------------------------------
if (USE_HIP)
  target_link_libraries(QuEST PUBLIC  ${HIP_PATH}/lib/libamdhip64.so )
elseif (USE_CUQUANTUM)
  find_library(CUQUANTUM_LIBRARIES custatevec)
  if (NOT CUQUANTUM_LIBRARIES)
    message(FATAL_ERROR "cuQuantum library (specifically custatevec) not found")
  endif ()

  target_link_libraries(QuEST ${CUDA_LIBRARIES} ${CUQUANTUM_LIBRARIES})
  target_include_directories(QuEST PUBLIC "/usr/local/cuda/include")
else()
  target_link_libraries(QuEST ${CUDA_LIBRARIES})
endif()

# ----- Coverage testing with GCC or Clang ------------------------------------
option(QUEST_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" FALSE)
if (QUEST_ENABLE_COVERAGE)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    message(STATUS "Configuring with coverage")
    target_compile_options(QuEST PUBLIC --coverage -O0)
    target_link_libraries(QuEST PUBLIC --coverage)
  else()
    message(FATAL_ERROR "GCC or Clang required with QUEST_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
  endif()
endif()

# ----- Memory testing with LLVM Clang ----------------------------------------
option(QUEST_MEMCHECK "Use LLVM AddressSanitizer to detecting memory errors" OFF)
if (QUEST_MEMCHECK)
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Configuring with LLVM AddressSanitizer")
    set(QUEST_MEMCHECK_FLAGS -fno-optimize-sibling-calls
        -fsanitize=address
        -fsanitize-address-use-after-scope
        )
    target_compile_options(QuEST PUBLIC -O1 -g -fno-omit-frame-pointer ${QUEST_MEMCHECK_FLAGS})
    target_link_libraries(QuEST PUBLIC -g ${QUEST_MEMCHECK_FLAGS})
  else()
    message(FATAL_ERROR "clang compiler required with QUEST_MEMCHECK: found ${CMAKE_CXX_COMPILER_ID}")
  endif()
endif()
