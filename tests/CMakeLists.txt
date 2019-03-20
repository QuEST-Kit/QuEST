# Calculate possible tests

option(ENABLE_ADVANCED_TEST   "Enable specific unit testing for advanced users" OFF)

if(${ENABLE_ADVANCED_TEST})
  option(ENABLE_TEST_UNIT       "Enable basic unit testing" OFF)
  option(ENABLE_TEST_BENCHMARKS "Enable testing benchmarks" OFF)
  option(ENABLE_TEST_SCIENCE    "Enable longer science test cases" OFF)

  execute_process(
    COMMAND ${TESTCOMMAND} -L
    OUTPUT_VARIABLE AVAILABLE_TESTS
    )
  string(REPLACE "\n" "\;" AVAILABLE_TESTS "${AVAILABLE_TESTS}")
  set(AVAILABLE_TESTS ${AVAILABLE_TESTS})
  
  foreach(TEST ${AVAILABLE_TESTS})
    add_test(${TEST}
      ${TESTCOMMAND} -Q ${QuEST_LIB_PATH} -p ${CMAKE_CURRENT_SOURCE_DIR}/unitPy ${TEST})
    set_tests_properties(${TEST} PROPERTIES  PASS_REGULAR_EXPRESSION "0 failed" FAIL_REGULAR_EXPRESSION "Error")
  endforeach(TEST)

else()
    add_test(all
      ${TESTCOMMAND} -Q ${QuEST_LIB_PATH} -p ${CMAKE_CURRENT_SOURCE_DIR}/unitPy unit)
    set_tests_properties(all PROPERTIES  PASS_REGULAR_EXPRESSION "0 failed" FAIL_REGULAR_EXPRESSION "Error") 
endif()
  
set(TESTCOMMAND ${CMAKE_CURRENT_SOURCE_DIR}/QuESTTest.py)

