# Python Test Suite

The testing suite is designed to be a flexible, but intuitive testing framework for the testing of the QuEST package. Most users will just want to run the standard tests through the CMake-Ctest interface which runs all unit tests (see main build docs), this is a brief guide for advanced users who might wish to run more comprehensive tests.

## One-line instant guide
To run the test suite for most users:
```
./QuESTTest.py unit algor
```
This will run all the tests considered unit (regression) tests and algorithm (functionality) tests in the `tests` folder structure. 

## Available tests

By default the test suite will search the `tests` directory for any matching `.test` files. Directories and sub-directories comprise super- and subsets of tests which are available to run. Any folder or `.test` filename may be specified on the command line call to run any tests matching that definition:
```
./QuESTTest.py hadamard
```
will run any test (`.test` file) or test sets (directory inc. sub-directories) called hadamard.

A list of available test files may be found by passing flags via the command line interface.

It is possible to further specify subsets of tests by using the directories as tags like so:
```
./QuESTTest.py density_matrix:hadamard
```
Which will run any test or test-set called "hadamard" with "density_matrix" in the path.  

## Test folder structure

The tests in the `tests` folder are structured first and foremost according to the type of test they represent. These are considered to be: essential, unit, algor and benchmarks.

- Essential tests consists of the tests without which the test suite would fail to operate. These are run at the start of any instantiation of the test suite whether or not explicitally requested and any failure in these sets means that the test suite will not run as results are not guaranteed.
- Unit tests are tests of single gate operations on initial state vectors.
- Algor tests are more advanced scientific tests which may take longer to run, but will test the conservation of various key physical/mathematical properties of combinations of gates to prove their veracity.
- Benchmark tests are large tests designed to strain the system to get reliable measures of speed and memory usage of certain machines, builds or changes.

## Logging
The test suite writes a log (by default to QuESTLog.log) of its status as it goes so that issues can be easily identified. By default it writes a log only on failures and (for MPI) on the root node. These can all be switched by flags.

## Flags
A list of available flags and options can be found by running:
```
./QuESTTest.py -h
```
It should be noted that a full description of how to use the test generation options can be found in the full developer's documentation which may be made available on request.

