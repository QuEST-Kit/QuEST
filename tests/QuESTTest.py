#!/usr/bin/env python3

from ctypes import *
from QuESTCore import *
from testset import *
from sys import argv
import argparse

parser = argparse.ArgumentParser(description='Python Test Suite for the QUantum Exact Simulation Toolkit (QuEST). ' )
parser.add_argument('-g','--generate', help='Generate a new set of benchmark tests', action='store_true')
parser.add_argument('-q','--quiet', help='Do not print results to screen', action='store_true')
parser.add_argument('-l','--logfile', help='Redirect log. DEFAULT=QuESTLog.log', default='QuESTLog.log')
parser.add_argument('-p','--testpath', help='Set test directory search path as colon-separated list. DEFAULT: unitPy/', default='unitPy/')
parser.add_argument('-t','--tolerance', type=float, help='Set the test failure tolerance for float values. DEFAULT=1e-10', default=1.e-10)
parser.add_argument('tests', nargs=argparse.REMAINDER, 
                    help='Set of tests one wishes to run, available default sets:\n'+", ".join(printSets)+" or any exposed QuEST function. DEFAULT=all")

argList = parser.parse_args()
if not argList.tests: argList.tests = ["all"]
print(argList)

testResults = init_tests(unitTestPath = argList.testpath, logFilePath = argList.logfile, tolerance = argList.tolerance, quiet = argList.quiet)

# If our argument is generate
if argList.generate:
    gen_tests(argList.tests)
    quit()

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

testsToRun = []
for test in argList.tests:
    testsToRun += testSets.get(test,[test])

for test in testsToRun:
    if test in testSets or test in tests:
        try:
            testResults.run_std_test(tests[test],test)
        except KeyError:
            print("Function '{}' does not exist, are you sure you wrote it correctly?".format(test))
    else:
        testResults.run_cust_test(test)

testResults.print_results()

finalise_tests()
