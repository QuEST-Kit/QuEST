add_subdirectory(src)

add_library(QuEST_lib
    ${QuEST_SRC}
)

target_include_directories(QuEST_lib 
    PRIVATE src
    PUBLIC include
)

target_compile_definitions(QuEST_lib
    PRIVATE
    QuEST_PREC=2
)

# set C/C++ compiler flags that are common between compilers 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} 
    -02 -std=c99"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
    -O2 -std=c++11 -Wall"
)

# set C/C++ compiler flags based on compiler type
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}
      -mavx -Wall"
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}
      -mavx"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}
      -mavx -Wall"
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}
      -mavx"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}
      -fprotect-parens -Wall -xAVX -axCORE-AVX2 -diag-disable cpu-dispatch"
  )
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}
      -xAVX -axCORE-AVX2 -diag-disable -cpu-dispatch"
  )
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio
endif()



