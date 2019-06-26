# -----------------------------------------------------------------------------
# ----- User Options ----------------------------------------------------------
# -----------------------------------------------------------------------------

option(ENABLE_ADVANCED_TEST "Enable specific testing for advanced users" OFF)

# -----------------------------------------------------------------------------
# ----- Test Generation Options -----------------------------------------------
# -----------------------------------------------------------------------------

set(TEST_DIR "${QuEST_Project_BINARY_DIR}/CustomTests" CACHE PATH "Testing files output directory")
set(TEST_NQUBITS 4 CACHE STRING "Number of qubits to run when running tests")
set(TEST_TYPES "PMS" CACHE STRING "Type of tests to generate")
set(TEST_QUREG "ZPDN" CACHE STRING "Types of qureg to generate in tests")
set(TEST_GENERATE "unit" CACHE STRING "Tests to generate on generate benchmarks")

# -----------------------------------------------------------------------------
# ----- Do not edit beyond here -----------------------------------------------
# -----------------------------------------------------------------------------

set(QuEST_TEST_DIR_NAME "tests" CACHE STRING
  "Name of the directory containing the QuEST library tests.")
set(QuEST_TEST_DIR "${QuEST_Project_SOURCE_DIR}/${QuEST_TEST_DIR_NAME}")

set(TESTCOMMAND ${PYTHON_EXECUTABLE} -m QuESTTest)

# Add target to generate benchmark results for this build
add_custom_target("benchmark"
  DEPENDS QuESTTest
  COMMAND mkdir -p ${TEST_DIR}
  COMMENT "Writing benchmark to ${TEST_DIR}
Test parameters: 
Testing ${TEST_GENERATE} with ${TEST_NQUBITS} qubits
Testing ${TEST_TYPES} on Quregs initialised in ${TEST_QUREG}"
  COMMAND ${TESTCOMMAND} -Q ${QuEST_LIB_PATH} ${INCLUDE_TESTS} -q -n ${TEST_NQUBITS} -T ${TEST_TYPES} -V ${TEST_QUREG} -p ${TEST_DIR} -g ${TEST_GENERATE}
  COMMAND ${CMAKE_COMMAND} -DBENCHMARKS_EXIST=ON ${QuEST_Project_SOURCE_DIR}
  WORKING_DIRECTORY ${UTIL_ROOT}
  )

if (${DISTRIBUTED})
  set(TESTCOMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${TESTCOMMAND})
endif()

# Do benchmarks exist? If so include test_dir
if (${BENCHMARKS_EXIST})
  set (INC_CUST "-p ${TEST_DIR}")
else()
  set (INC_CUST "")
endif()

# Add Test to build core QuEST
add_test(NAME build
  COMMAND ${CMAKE_MAKE_PROGRAM} "QuESTTest"
  WORKING_DIRECTORY ${QuEST_Project_BINARY_DIR}
  )

# Set up test for set "name"
macro(Testee name)
  add_test(
    NAME ${name}
    COMMAND ${TESTCOMMAND} ${INC_CUST} -Q ${QuEST_LIB_PATH} -l
    ${QuEST_Project_BINARY_DIR}/QuESTLog.log ${name}
    WORKING_DIRECTORY ${UTIL_ROOT}
    )
  set_tests_properties(${name}
    PROPERTIES
    DEPANDS build
    PASS_REGULAR_EXPRESSION " 0 failed"
    FAIL_REGULAR_EXPRESSION "Error")
  add_custom_target("test_${name}"
    DEPENDS QuESTTest
    DEPENDS QuEST
    COMMAND ctest -R ${name}
    WORKING_DIRECTORY ${QuEST_Project_BINARY_DIR}
    )
endmacro(Testee)
  
if(${ENABLE_ADVANCED_TEST})
  # Get list of test sets
  execute_process(
    COMMAND ${TESTCOMMAND} -Q ${QuEST_LIB_PATH} -L
    OUTPUT_VARIABLE AVAILABLE_TESTS
    WORKING_DIRECTORY ${UTIL_ROOT}
    )
  # Make semi-colon separated list as required by Cmake
  string(REPLACE "\n" "\;" AVAILABLE_TESTS "${AVAILABLE_TESTS}")
  # Cast string into list
  set(AVAILABLE_TESTS ${AVAILABLE_TESTS})

  # Add all tests into available test list
  message(STATUS "Advanced testing is ${ENABLE_ADVANCED_TEST}")
  foreach(TEST ${AVAILABLE_TESTS})
    Testee(${TEST})
    message(STATUS "Found test set: ${TEST}")
  endforeach(TEST)

else()
  # Add unit test into tests (Core)
  Testee(unit)
endif()
