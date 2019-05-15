#!/usr/bin/env python3

# Import libraries needed for initialisation
import os.path
import argparse
from QuESTPy.QuESTBase import init_QuESTLib

from QuESTPy.QuESTLibDir import defaultQuESTDir
# If set up
QuESTPath = defaultQuESTDir or "../build/QuEST"

parser = argparse.ArgumentParser(description='Python test suite for the Quantum Exact Simulation Toolkit (QuEST).',epilog='''NOTE: Tests can be specified as full filepaths or testnames, which will be searched for in the TESTPATH, with the earliest path taking priority. 

Tests with a full filepath can have ".test" or ".py" extensions. 

Custom .test files can be found as TESTPATH/TESTS.test or by a full filepath. 
Custom .py files must be specified by a full filepath. ''', add_help=False)

mutExc = parser.add_mutually_exclusive_group()

# Need to pull some trickery to allow QuESTLib redirection. Probably cleaner way to do this, but...
mutExc.add_argument('-h','--help', help="Show this help message and exit", action='store_true')
parser.add_argument('-Q','--questpath', help="Define alternative QuEST library location. The library must be named 'libQuEST.so' located in the specified directory to be found. Default=%(default)s", default=QuESTPath)

# Just parse -Q
QuESTPath = parser.parse_known_args()

# Dummy printsets

try:
    # Load QuEST Library if not printing help
    init_QuESTLib(QuESTPath[0].questpath)

    # Import remaining libraries
    from QuESTTest.QuESTCore import *
except FileNotFoundError:
    parser._actions[0].help = argparse.SUPPRESS
    parser.epilog = ""
    parser.description = ""
    print("Library not found, here's how to redirect QuESTTest to find library:\n")
    parser.print_help()
    print()
    raise 

mutExc.add_argument('-L','--list-tests', help="Print available test sets and exist", action='store_true')
mutExc.add_argument('-A','--list-all-tests', help="Print available tests and exist", action='store_true')
parser.add_argument('-q','--quiet', help='Do not print results to screen', action='store_true')
parser.add_argument('-l','--logfile', help='Redirect log. DEFAULT=%(default)s', default='QuESTLog.log')
parser.add_argument('-p','--testpath', help='Set test directory search path as colon-separated list. DEFAULT=essential:algor:benchmarks:unit', default = "")
parser.add_argument('-t','--tolerance', type=float, help='Set the test failure tolerance for float values. DEFAULT=%(default)s', default=1.e-10)
parser.add_argument('-f','--mpilog', help='Full MPI logging on a per-process basis, creates a new file for each process of "<LOGFILE>.<MPIRANK>" . Default=False', action='store_true')
parser.add_argument('-F','--full-log', help='Log all passes as well as failures. Default=False', action='store_true')
testArg = parser.add_argument('tests', nargs=argparse.REMAINDER, metavar="TESTS",
                    help="Set of tests one wishes to run, this can be any, any custom test (see NOTE) or any exposed QuEST function. DEFAULT=unit")
genGroup = parser.add_argument_group('Generation', 'Arguments related to the generation of tests')
genGroup.add_argument('-g','--generate', help='Generate a new set of benchmark tests for tests listed redirected to TESTPATH.', action='store_true')
genGroup.add_argument('-n','--numqubits', type=int, help='Specify the number of qubits to generate on generation. DEFAULT=%(default)s', default=3)
genGroup.add_argument('-c','--checks', help='Specify the checks to be generated. P: Total probability, M: Probability of each qubit being in 0 or 1 state, S: Full State Vector, as a single string. DEFAULT=%(default)s', default='PMS')
genGroup.add_argument('-V','--quregtypes', help='Specify which types of Quregs are generated in the tests. Z: Zero state, P: Plus state, D: Debug state, R: Random state, N: Normalised random state. States can be multiply defined, e.g. \'RRR\' will generate 3 different random configurations. DEFAULT=%(default)s', default='ZPDRzpdn')
genGroup.add_argument('-T','--targets', help='Specify a comma-separated list of targets for gates to operate on, or specify "E" for each qubit. DEFAULT=%(default)s', default="E")
genGroup.add_argument('-C','--controls', help='Specify a comma-separated list of controls for control gates to use, or specify "E" for each qubit. DEFAULT=%(default)s', default="E")

argList = parser.parse_args()

# Set up Parallel environment and testing framework
init_tests(unitTestPath = argList.testpath, logFilePath = argList.logfile, tolerance = argList.tolerance, quiet = argList.quiet, mpiLog=argList.mpilog, fullLog=argList.full_log)

# Now we manually handle the print with *all* potential arguments included

# Set default for the tests to run
if not argList.tests: argList.tests = ["unit"]

if any( [argList.generate, argList.list_tests, argList.list_all_tests, argList.help] ):
    if argList.help:
        testList = ", ".join(list_test_sets())
        testArg.help = "Set of tests one wishes to run, this can be any of: "+testList+", any custom test (see NOTE) or any exposed QuEST function. DEFAULT=all"
        if root: parser.print_help()
    
    if argList.list_tests:
        for test in list_test_sets():
            print(test)
    if argList.list_all_tests:
        for test in list_all_tests():
            print(test)
    
    # If our argument is generate
    if argList.generate:
        repEnv()
        if not argList.testpath:
            print('Will not overwrite core tests, specify output directory with -p')
            quit()

        if root:
            testResults.gen_tests(testsToGen = argList.tests, nQubits = argList.numqubits,
                                  qubitGen = argList.quregtypes, testGen = argList.checks,
                                  targScan = argList.targets, contScan = argList.controls)
    quit()

repEnv()
testResults.run_tests(["essential"])
if testResults.fails > 0:
    finalise_tests()
    raise ValueError("System failed essential qubit initalisation tests, impossible to continue!")

testResults.run_tests(argList.tests)
    
# Print final answer
finalise_tests()
