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
target_compile_definitions(tests PRIVATE TEST_ALL_DEPLOYMENTS=$<BOOL:${TEST_ALL_DEPLOYMENTS}>)


# when QuEST is MPI-compiled, launch CTest distributed,
# with as many nodes as set in var NUM_CTEST_MPI_PROCS.
# note that building/compilation will fail if this exceeds
# the number of available system slots which is a nuisance
# if compiling on a smaller system with the intention of
# running ctest on a larger one. Alas the --oversubscribe
# flag is OpenMPI specific (and not necessary at all with
# e.g. MPICH) so we forego solving this for now
if (ENABLE_DISTRIBUTION)

  if (DEFINED NUM_CTEST_MPI_PROCS)

    ## works on clang, fails to compile on Fritz and ARCHER2 (missing linker symbols)
    # set_property(
    #   TARGET tests 
    #   PROPERTY CROSSCOMPILING_EMULATOR 
    #   '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_CTEST_MPI_PROCS}')



    # didn't do anything on Fritz
    #set_property(TARGET tests PROPERTY TEST_LAUNCHER ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_CTEST_MPI_PROCS})
    #set_property(TARGET tests PROPERTY TEST_EXECUTOR ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_CTEST_MPI_PROCS})

    # when paired with set_property() below, it also did nothing on MacOS


    # but potentially because it is missing the below property which needs to be applied PER TEST
    # set_property(TEST tests PROPERTY PROCESSORS ${NUM_CTEST_MPI_PROCS})



    message(STATUS "CTests will be launched using ${NUM_CTEST_MPI_PROCS} MPI processes.")
  else()
    message(WARNING "Distributed tests were compiled but NUM_CTEST_MPI_PROCS was not passed, so CTest will subsequently run serially.")
  endif()

endif()


# call after MPI setup above to distribute CTest
# catch_discover_tests(tests)



# catch_discover_tests(tests OUTPUT_VARIABLE TEST_NAMES)

# foreach(test ${TEST_NAMES})
#     set_property(TEST ${test} PROPERTY PROCESSORS ${NUM_CTEST_MPI_PROCS})
# endforeach()



# works sometimes, but discovers less than 1/4th of the tests, and
# strangely fails on Fritz and ARCHER2 when ctest is run (though passes
# on LLVM)
# include(CTest)
# include(ParseAndAddCatchTests)

# set(OptionalCatchTestLauncher ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4)
# ParseAndAddCatchTests(tests)

catch_discover_tests(tests OUTPUT_VARIABLE TEST_NAMES)