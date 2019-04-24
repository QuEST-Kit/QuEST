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


option(MULTITHREADED "Switches on OpenMP Threading if an OpenMP implementation can be found. Set to 0 to disable." 1)

option(DISTRIBUTED "Whether the program will be run in distributed mode using MPI. Set to 1 to enable" 0)

set(PRECISION 2 CACHE STRING
    "Whether to use single, double or quad floating point precision in the state vector. {1,2,4}")

option(GPUACCELERATED "Whether to program will run on GPU. Set to 1 to enable" 0)

set(GPU_COMPUTE_CAPABILITY 30 CACHE STRING "GPU hardware dependent, lookup at https://developer.nvidia.com/cuda-gpus. Write without fullstop")



# *****************************************************************************
# ***** NO CHANGES SHOULD BE REQUIRED FROM THE USER BEYOND THIS POINT *********
# *****************************************************************************

# Show the user their settings
message(STATUS "GPU acceleration is ${GPUACCELERATED}")
message(STATUS "OMP acceleration is ${MULTITHREADED}")
message(STATUS "MPI distribution is ${DISTRIBUTED}")


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
    find_package(CUDA REQUIRED)
endif()



# -----------------------------------------------------------------------------
# ----- SET COMPILER FLAGS ----------------------------------------------------
# -----------------------------------------------------------------------------


# ----- OPENMP ----------------------------------------------------------------

if (${MULTITHREADED} AND NOT ${GPUACCELERATED})
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# ----- CUDA FLAGS ------------------------------------------------------------

if (GPUACCELERATED)
    set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
        -arch=compute_${GPU_COMPUTE_CAPABILITY} -code=sm_${GPU_COMPUTE_CAPABILITY}"
    )
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
            -O2"
        )
    endif()
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
            -G -g -lineinfo"
        )
    else()
        set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
            -O2"
        )
    endif()
endif()


# ----- C COMPILER FLAGS --------------------------------------------------

# set C flags that are common between compilers and build types
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -std=c99"
)

# Use -O2 for all but debug mode by default 
if (NOT("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
        -O2"
    )
endif()

# Set c flags for release
set(CMAKE_C_FLAGS_RELEASE 
    "-O2"
)

# Set c flags to use in debug mode
set(CMAKE_C_FLAGS_DEBUG 
    "-g"
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
endif()

# ----- C++ COMPILER FLAGS --------------------------------------------------

# set C++ flags that are common between compilers and build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -std=c++98 -Wall"
)

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
set(CMAKE_CXX_FLAGS_DEBUG 
    "-g"
)

# set C++ compiler flags based on compiler type
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch"
  )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
endif()

if (VERBOSE_CMAKE)
    message("-- Compiler flags:")
    message("   C Compiler ID: ${CMAKE_C_COMPILER_ID}")
    message("   C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
    message("   C compilation flags: ${CMAKE_C_FLAGS}")
    message("   CUDA compilation flags: ${CUDA_NVCC_FLAGS}")
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

if (GPUACCELERATED)
    cuda_add_library(QuEST SHARED
        ${QuEST_SRC}
    )
else ()
    add_library(QuEST SHARED
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
    PRIVATE
    QuEST_PREC=${PRECISION}
)

# -----------------------------------------------------------------------------
# ----- LINK LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

# ----- MPI -------------------------------------------------------------------

target_link_libraries(QuEST ${MPI_C_LIBRARIES})

# ----- GPU -------------------------------------------------------------------

target_link_libraries(QuEST ${CUDA_LIBRARIES})



