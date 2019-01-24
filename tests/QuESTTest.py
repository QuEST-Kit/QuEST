#!/usr/bin/env python3

from ctypes import *
from QuESTCore import *
from sys import argv

testResults = init_tests(unitTestPath = 'unitPy/', logFilePath = 'QuESTLog.log', tolerance = 1.e-10)


# If our argument is generate
if len(argv) > 1 and argv[1] == "generate":
    gen_tests()
    quit()
    
# Load the test sets
from testset import *

for test in testSets["essential"]:
    if test in testSets or test in tests:
        try:
            testResults.run_std_test(tests[test],test)
        except KeyError:
            print("Function '{}' does not exist, are you sure you wrote it correctly?".format(test))
    else:
        testResults.run_cust_test(test)
if testResults.fails > 0:
    raise ValueError("System failed essential qubit initalisation tests, impossible to continue!")


if (len(argv) == 1): # Default to "all"
    testsToRun = testSets['all']
else: # Run through command line args
    testsToRun = []
    for test in argv[1:]:
        testsToRun += testSets.get(test,[test])


for test in testsToRun:
    if test in testSets or test in tests:
        if test in specSets:
            for testFile in testSets[test]:
                testResults.run_cust_test(testFile)
        else:
            try:
                testResults.run_std_test(tests[test],test)
            except KeyError:
                print("Function '{}' does not exist, are you sure you wrote it correctly?".format(test))
    else:
        testResults.run_cust_test(test)

testResults.print_results()

finalise_tests()
