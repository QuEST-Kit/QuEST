# CMake generated Testfile for 
# Source directory: /home/izemize/postdoc/2019/quEST/utilities
# Build directory: /home/izemize/postdoc/2019/quEST/build/utilities
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(build "/usr/bin/make" "QuESTTest")
set_tests_properties(build PROPERTIES  WORKING_DIRECTORY "/home/izemize/postdoc/2019/quEST/build")
add_test(unit "/usr/bin/python3.6" "-m" "QuESTTest" "-Q" "/home/izemize/postdoc/2019/quEST/build/QuEST" "-l" "/home/izemize/postdoc/2019/quEST/build/QuESTLog.log" "unit")
set_tests_properties(unit PROPERTIES  DEPANDS "build" FAIL_REGULAR_EXPRESSION "Error" PASS_REGULAR_EXPRESSION " 0 failed" WORKING_DIRECTORY "/home/izemize/postdoc/2019/quEST/utilities")
