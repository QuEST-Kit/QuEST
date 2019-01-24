# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ----- QuEST LIBRARY BUILD SYSTEM --------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


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
    set(MULTITHREADED 1)
endif()

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
    set(GPU_COMPUTE_CAPABILITY 30)
endif()

message("\n\nmultithreaded\n\n ${MULTITHREADED}")


# -----------------------------------------------------------------------------
# ----- BUILD LIBRARY ---------------------------------------------------------
# -----------------------------------------------------------------------------

add_subdirectory(src)

add_library(QuEST_lib SHARED
    ${QuEST_SRC}
)

# -----------------------------------------------------------------------------
# ----- CONFIGURE LIBRARY -----------------------------------------------------
# -----------------------------------------------------------------------------

# ----- Location of header files ----------------------------------------------

target_include_directories(QuEST_lib 
    PRIVATE src
    PUBLIC include
)


# ----- Macros -----------------------------------------------------

target_compile_definitions(QuEST_lib
    PRIVATE
    QuEST_PREC=${PRECISION}
)


# ----- OPENMP ----------------------------------------------------------------

if (${MULTITHREADED})
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# ----- COMPILER FLAGS --------------------------------------------------------

# set C/C++ compiler flags that are common between compilers 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -O2 -std=c99"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -O2 -std=c++11 -Wall"
)

# TODO standardize
# set C/C++ compiler flags based on compiler type
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -Wall"
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} \
    -mavx"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -Wall"
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} \
    -mavx"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable cpu-dispatch"
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} \
    -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
endif()



