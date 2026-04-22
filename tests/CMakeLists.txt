# @author Oliver Thomson Brown

add_executable(tests
  main.cpp
)
target_link_libraries(tests PRIVATE QuEST::QuEST Catch2::Catch2)
target_compile_features(tests PUBLIC cxx_std_20)

add_subdirectory(unit)
add_subdirectory(utils)
add_subdirectory(integration)

if (ENABLE_DEPRECATED_API)
  add_subdirectory(deprecated)
endif()


# let Catch2 register all tests with CTest
catch_discover_tests(tests)
