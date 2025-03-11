# @author Oliver Thomson Brown
# @author Tyson Jones (MPI ctests)

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
if (DEFINED FAST_UNIT_TESTS)
  target_compile_definitions(tests PRIVATE FAST_UNIT_TESTS=${FAST_UNIT_TESTS})
endif()


# when QuEST is MPI-compiled, launch CTest distributed,
# with as many nodes as set in var NUM_CTEST_MPI_PROCS
if (ENABLE_DISTRIBUTION)

  if (DEFINED NUM_CTEST_MPI_PROCS)
    set_property(
      TARGET tests 
      PROPERTY CROSSCOMPILING_EMULATOR 
      '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_CTEST_MPI_PROCS}')
    message(STATUS "CTests will be launched using ${NUM_CTEST_MPI_PROCS} MPI processes.")
  else()
    message(WARNING "Distributed tests were compiled but NUM_CTEST_MPI_PROCS was not passed, so CTest will subsequently run serially.")
  endif()

endif()


catch_discover_tests(tests)
