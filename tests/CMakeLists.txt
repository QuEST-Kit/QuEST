add_subdirectory(catch)

add_executable(tests
  test_main.cpp
  test_utilities.cpp
)
target_link_libraries(tests PUBLIC QuEST::QuEST Catch2)

add_subdirectory(unit)
add_subdirectory(utils)

catch_discover_tests(tests)