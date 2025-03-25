# @author Oliver Thomson Brown
# @author Tyson Jones (macros and MPI ctests)

add_executable(tests
  main.cpp
)
target_link_libraries(tests PRIVATE QuEST::QuEST Catch2::Catch2)

add_subdirectory(unit)
add_subdirectory(utils)
add_subdirectory(integration)

if (ENABLE_DEPRECATED_API)
  add_subdirectory(deprecated)
endif()


# map test-related cmake vars to preprocessors
target_compile_definitions(tests PRIVATE TEST_MAX_NUM_QUBIT_PERMUTATIONS=${TEST_MAX_NUM_QUBIT_PERMUTATIONS})
target_compile_definitions(tests PRIVATE TEST_MAX_NUM_SUPEROP_TARGETS=${TEST_MAX_NUM_SUPEROP_TARGETS})
target_compile_definitions(tests PRIVATE TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS=${TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS})
target_compile_definitions(tests PRIVATE TEST_ALL_DEPLOYMENTS=$<BOOL:${TEST_ALL_DEPLOYMENTS}>)


# let Catch2 register all tests with CTest
catch_discover_tests(tests)
