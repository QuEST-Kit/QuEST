add_executable(tests
  main.cpp
)
target_link_libraries(tests PUBLIC QuEST::QuEST Catch2::Catch2)

add_subdirectory(unit)
add_subdirectory(utils)

if (ENABLE_DEPRECATED_API)
  add_subdirectory(deprecated)
endif()

catch_discover_tests(tests)