# author: Tyson Jones

# locate Catch2 header for compiler
include_directories(catch)

# compile Catch2 unit tests in C++14
add_executable(tests
    catch_main.cpp
    QuEST_test_utils.cpp
    QuEST_test_calculations.cpp
    QuEST_test_data_structures.cpp
    QuEST_test_decoherence.cpp
    QuEST_test_gates.cpp
    QuEST_test_operators.cpp
    QuEST_test_state_initialisations.cpp
    QuEST_test_unitaries.cpp
)
set_property(TARGET tests PROPERTY CXX_STANDARD 14)
if (WIN32)
    target_link_libraries(tests QuEST)
else ()
    target_link_libraries(tests QuEST m)
endif ()

# locate Catch2 scripts for CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/catch)
include(CTest)
include(Catch)

# register Catch2 tests as individual CTest tests
catch_discover_tests(tests)