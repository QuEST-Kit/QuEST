#!/usr/bin/env python3

# Import libraries needed for initialisation
import argparse
from QuESTBase import init_QuESTLib

parser = argparse.ArgumentParser(description='Python test suite for the Quantum Exact Simulation Toolkit (QuEST).',epilog='''NOTE: Tests can be specified as full filepaths or testnames, which will be searched for in the TESTPATH, with the earliest path taking priority. 

Tests with a full filepath can have ".test" or ".py" extensions. 

Custom .test files can be found as TESTPATH/TESTS.test or by a full filepath. 
Custom .py files must be specified by a full filepath. ''', add_help=False)

# Need to pull some trickery to allow QuESTLib redirection. Probably cleaner way to do this, but...
parser.add_argument('-h','--help', help="Show this help message and exit", action='store_true')
parser.add_argument('-Q','--questpath', help="Define alternative QuEST library location. The library must be named 'libQuEST.so' to be found. Default=%(default)s", default='../build/QuEST')

# Just parse -Q
QuESTPath = parser.parse_known_args()

init_QuESTLib(QuESTPath[0].questpath)

# Import remaining libraries
from QuESTCore import *
from testset import *

parser.add_argument('-q','--quiet', help='Do not print results to screen', action='store_true')
parser.add_argument('-l','--logfile', help='Redirect log. DEFAULT=%(default)s', default='QuESTLog.log')
parser.add_argument('-p','--testpath', help='Set test directory search path as colon-separated list. DEFAULT=%(default)s', default='unitPy/')
parser.add_argument('-t','--tolerance', type=float, help='Set the test failure tolerance for float values. DEFAULT=%(default)s', default=1.e-10)
parser.add_argument('-g','--generate', help='Generate a new set of benchmark tests for tests listed redirected to TESTPATH.', action='store_true')
parser.add_argument('-n','--numqubits', type=int, help='Specify the number of qubits to generate on generation. DEFAULT=%(default)s', default=3)
parser.add_argument('-f','--mpilog', help='Full MPI logging on a per-process basis, creates a new file for each process of "<LOGFILE>.<MPIRANK>" . Default=False', action='store_true')
parser.add_argument('tests', nargs=argparse.REMAINDER, metavar="TESTS",
                    help='Set of tests one wishes to run, available default sets:'+", ".join(printSets)+", any custom test (see NOTE) or any exposed QuEST function. DEFAULT=all")

argList = parser.parse_args()

if not argList.tests: argList.tests = ["all"]

if argList.help:
    if root: parser.print_help()
    quit()

# Set up Parallel environment and testing framework
init_tests(unitTestPath = argList.testpath, logFilePath = argList.logfile, tolerance = argList.tolerance, quiet = argList.quiet, fullLogging=argList.mpilog)

# If our argument is generate
if argList.generate:
    if root: gen_tests(argList.tests, argList.numqubits)
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

finalise_tests()
