# @author Oliver Thomson Brown
# @author Erich Essmann (patched MPI)

add_executable(dep_tests
    test_main.cpp
    test_calculations.cpp
    test_data_structures.cpp
    test_decoherence.cpp
    test_gates.cpp
    test_operators.cpp
    test_state_initialisations.cpp
    test_unitaries.cpp
    test_utilities.cpp
)
target_link_libraries(dep_tests PUBLIC QuEST::QuEST Catch2::Catch2)
target_compile_definitions(dep_tests PRIVATE INCLUDE_DEPRECATED_FUNCTIONS=1)

if (ENABLE_DISTRIBUTION)
    target_link_libraries(dep_tests PRIVATE MPI::MPI_CXX)
endif()

catch_discover_tests(dep_tests)