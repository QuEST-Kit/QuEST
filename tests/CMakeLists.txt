# author: Tyson Jones

# locate Catch2 header for compiler
include_directories(catch)

# compile Catch2 unit tests in C++14
add_executable(tests
    test_main.cpp
    test_utilities.cpp
    test_calculations.cpp
    test_data_structures.cpp
    test_decoherence.cpp
    test_gates.cpp
    test_operators.cpp
    test_state_initialisations.cpp
    test_unitaries.cpp
)

# link QuEST
target_link_libraries(tests QuEST)

if (ENABLE_DISTRIBUTION)
    # distributed unit tests must supply an MPI-clue to (modified) catch
    add_definitions(-DDISTRIBUTED_MODE)
endif ()

# locate Catch2 scripts for CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/catch)
include(CTest)
include(Catch)

# register Catch2 tests as individual CTest tests
catch_discover_tests(tests)