# @author Oliver Thomson Brown
# @author Tyson Jones (patched MSVC and test discovery)

add_executable(tests
  main.cpp
)
target_link_libraries(tests PRIVATE QuEST::QuEST Catch2::Catch2)
target_compile_features(tests PUBLIC cxx_std_20)

if (QUEST_ENABLE_MPI AND QUEST_ENABLE_SUBCOMM)
  target_link_libraries(tests PRIVATE MPI::MPI_CXX)
endif()

# extend the MSVC max object size
if (MSVC)
  target_compile_options(tests PRIVATE /bigobj)
endif()

add_subdirectory(unit)
add_subdirectory(utils)
add_subdirectory(integration)

if (QUEST_ENABLE_DEPRECATED_API)
  add_subdirectory(deprecated)
endif()

# defer test discovery, so that (e.g.) MPI libs aren't loaded during build
catch_discover_tests(tests DISCOVERY_MODE PRE_TEST)
