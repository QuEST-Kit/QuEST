# Calculate possible tests

option(GENERATE_BENCHMARK "Generate set of unit tests for this build" OFF)
option(ENABLE_ADVANCED_TEST "Enable specific unit testing for advanced users" OFF)
set(TEST_NQUBITS 4 CACHE STRING "Number of qubits to run when running tests")
set(TEST_TYPES "PMS" CACHE STRING "Type of tests to generate")
set(TEST_QUREG "ZPDN" CACHE STRING "Types of qureg to generate in tests")
set(TEST_GENERATE "unit" CACHE STRING "Tests to generate on generate benchmarks")

set(TEST_DIR "${QuEST_Project_BINARY_DIR}/CustomTests" CACHE PATH "Testing files output directory")
set(TESTCOMMAND "${CMAKE_CURRENT_SOURCE_DIR}/QuESTTest.py")


add_custom_target("benchmark"
  DEPENDS QuEST
  COMMAND ${CMAKE_COMMAND} -P  ${QuEST_Project_SOURCE_DIR} -DGENERATE_BENCHMARK=ON
  )

if (GENERATE_BENCHMARK)
  message(STATUS "Generating benchmarking test results")
  message(STATUS "${TEST_NQUBITS} Qubits, ${TEST_TYPES}, ${TEST_QUREG}, ${TEST_GENERATE}")
  set(GENERATE_BENCHMARK OFF CACHE BOOL "Generate set of unit tests for this build" FORCE)
  execute_process(
    COMMAND mkdir -p ${TEST_DIR}
    COMMAND ${PYTHON_EXECUTABLE} ${TESTCOMMAND} -Q ${QuEST_LIB_PATH} -p ${TEST_DIR} -q -n ${TEST_NQUBITS} -T ${TEST_TYPES} -V ${TEST_QUREG} -g ${TEST_GENERATE}
    )
  message(STATUS "Generating benchmarking test results -- DONE") 
endif()
  
if(${ENABLE_ADVANCED_TEST})

  execute_process(
    COMMAND ${TESTCOMMAND} -L
    OUTPUT_VARIABLE AVAILABLE_TESTS
    )
  string(REPLACE "\n" "\;" AVAILABLE_TESTS "${AVAILABLE_TESTS}")
  set(AVAILABLE_TESTS ${AVAILABLE_TESTS})
  
  foreach(TEST ${AVAILABLE_TESTS})
    add_test(${TEST}
      ${PYTHON_EXECUTABLE} ${TESTCOMMAND} -Q ${QuEST_LIB_PATH} ${TEST}) # -p ${TEST_DIR}
    set_tests_properties(${TEST} PROPERTIES  PASS_REGULAR_EXPRESSION " 0 failed" FAIL_REGULAR_EXPRESSION "Error")
    add_custom_target("test_${TEST}"
      DEPENDS QuEST
      COMMAND ctest -R ${TEST}
      )
  endforeach(TEST)

else()
    add_test(all
      ${PYTHON_EXECUTABLE} ${TESTCOMMAND} -p ${TEST_DIR} -Q ${QuEST_LIB_PATH} unit) #
    set_tests_properties(all PROPERTIES  PASS_REGULAR_EXPRESSION " 0 failed" FAIL_REGULAR_EXPRESSION "Error") 
endif()
  


