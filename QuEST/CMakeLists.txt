# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----- QuEST LIBRARY BUILD SYSTEM --------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Builds QuEST as a shared library libQuEST.so


# -----------------------------------------------------------------------------
# ----- LIBRARY USER SETTINGS -------------------------------------------------
# -----------------------------------------------------------------------------

# These defaults can also be overridden using:
#   cmake -D[VAR]=[VALUE] [Path to root CMakeLists.txt]
# eg
#   cmake -DDISTRIBUTED=1 ..

# TODO -- this should only be part of debug mode. User shouldn't care.
# Switches on OpenMP threading if an OpenMP implementation can be found.
# Set to 0 to disable
if (NOT DEFINED MULTITHREADED)
    #TODO change back
    set(MULTITHREADED 1)
endif()

# Whether the program will be run in distributed mode using MPI
# Set to 1 to enable
if (NOT DEFINED DISTRIBUTED)
    set(DISTRIBUTED 0)
endif()

# Whether to use single, double or quad floating point precision in the 
# state-vector {1,2,4}
if (NOT DEFINED PRECISION)
    set(PRECISION 2)
endif()

# GPU config. Set to 1 to enable
if (NOT DEFINED GPUACCELERATED)
    set(GPUACCELERATED 0)
endif()
# GPU hardware dependent, lookup at https://developer.nvidia.com/cuda-gpus, write without fullstop
if (NOT DEFINED GPU_COMPUTE_CAPABILITY)
    set(GPU_COMPUTE_CAPABILITY 52)
endif()

# -----------------------------------------------------------------------------
# ----- FIND PACKAGES ---------------------------------------------------------
# -----------------------------------------------------------------------------

if (DISTRIBUTED)
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
endif()

if (GPUACCELERATED)
    find_package(CUDA REQUIRED)
endif()


# -----------------------------------------------------------------------------
# ----- SET COMPILER FLAGS ----------------------------------------------------
# -----------------------------------------------------------------------------


# ----- OPENMP ----------------------------------------------------------------

if (${MULTITHREADED})
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# ----- CUDA FLAGS ------------------------------------------------------------

if (GPUACCELERATED)
    set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \
        -O2 -arch=compute_${GPU_COMPUTE_CAPABILITY} -code=sm_${GPU_COMPUTE_CAPABILITY}"
    )
endif()


# ----- C COMPILER FLAGS --------------------------------------------------

# set C flags that are common between compilers 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -O2 -std=c99"
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

# set C++ flags that are common between compilers 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -O2 -std=c++11 -Wall"
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
    message("\n")
    message("-- Compiler flags:")
    message("   C Compiler ID: ${CMAKE_C_COMPILER_ID}")
    message("   C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
    message("   C compilation flags: ${CMAKE_C_FLAGS}")
    message("   CUDA compilation flags: ${CUDA_NVCC_FLAGS}")
    message("   C++ compilation flags: ${CMAKE_CXX_FLAGS}")
    message("\n")
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



