# @author Oliver Thomson Brown

add_executable(tests
  main.cpp
)
target_link_libraries(tests PRIVATE QuEST::QuEST Catch2::Catch2)
target_compile_features(tests PUBLIC cxx_std_20)

# DEBUG
if (MSVC)
  # add_compile_options(/openmp:llvm /EHsc /bigobj /Zm500)
  # add_compile_definitions(CATCH_CONFIG_DISABLE_STRINGIFICATION)
  # add_compile_definitions(CATCH_CONFIG_FAST_COMPILE)

  # see if openmp is reaching tests or Catch2
  get_target_property(test_opts tests COMPILE_OPTIONS)
  message(STATUS "tests COMPILE_OPTIONS: ${test_opts}")
  get_target_property(catch_opts Catch2 COMPILE_OPTIONS)
  message(STATUS "Catch2 COMPILE_OPTIONS: ${catch_opts}")

  # try force off
  target_compile_options(tests PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/openmp->
  )
  target_compile_options(Catch2 PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/openmp->
  )
endif()

add_subdirectory(unit)
add_subdirectory(utils)
add_subdirectory(integration)

if (ENABLE_DEPRECATED_API)
  add_subdirectory(deprecated)
endif()


# let Catch2 register all tests with CTest
catch_discover_tests(tests)
