# author: Tyson Jones

# locate Catch2 header for compiler
include_directories(catch)

# compile Catch2 unit tests in C++14
add_executable(tests
    main.cpp
    utilities.cpp
    test_calculations.cpp
    test_data_structures.cpp
    test_decoherence.cpp
    test_gates.cpp
    test_operators.cpp
    test_state_initialisations.cpp
    test_unitaries.cpp
)
set_property(TARGET tests PROPERTY CXX_STANDARD 14)

# link QuEST
if (WIN32)
    target_link_libraries(tests QuEST)
else ()
    target_link_libraries(tests QuEST m)
endif ()

if (${DISTRIBUTED})
    # distributed unit tests must supply an MPI-clue to (modified) catch
    add_definitions(-DDISTRIBUTED_MODE)
    
    # compile and link C++ MPI
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(tests ${MPI_CXX_LIBRARIES})
endif ()

# locate Catch2 scripts for CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/catch)
include(CTest)
include(Catch)

# register Catch2 tests as individual CTest tests
catch_discover_tests(tests)