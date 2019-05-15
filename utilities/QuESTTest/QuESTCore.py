import os.path
import os
from QuESTPy.QuESTFunc import *
from functools import total_ordering
import importlib.util
import importlib.machinery

# Add .test to valid python names
importlib.machinery.SOURCE_SUFFIXES += ['.test']

@total_ordering
class TestSet:
    def __init__(self, name, depth, path, tests):
        self.name = name
        self.depth = depth
        self.path = path
        self.tests = tests

    def __eq__(self, other):
        return self.depth == other.depth

    def __lt__(self, other):
        return self.depth < other.depth

    def __repr__(self):
        return str(( self.name, self.depth, self.path, self.tests ))

    def names(self):
        return ( os.path.basename(os.path.splitext(test)[0]) for test in self.tests)

class SuperTestSet(list):
    def __call__(self, target):
        *specs, name  = target.split(":")
        if name == "all": return SuperTestSet([ testSet for testSet in self.depth(1)])
        return SuperTestSet([ testSet for testSet in self
                              if testSet.name == name and all( [path in testSet.path for path in specs] ) ])

    def depth(self, depth):
        if isinstance(depth, int) or depth is None: mindepth = maxdepth = depth
        elif len(depth) == 2: mindepth, maxdepth = depth
        else: raise TypeError('Unexpected type in SuperTestSet.depth')

        if mindepth is None: mindepth = 0
        if maxdepth is None: maxdepth = 999

        return SuperTestSet([ testSet for testSet in self
                              if testSet.depth >= mindepth and testSet.depth <= maxdepth ])

    def names(self, depth = None):
        for testSet in self.depth( depth ):
            yield testSet.name

    def tests(self):
        for testSet in self:
            for test in testSet.tests:
                yield test

    def paths(self):
        for testSet in self:
            yield testSet.path

    def dump(self, depth = None):
        for testSet in sorted(self.depth( depth )):
            print(testSet.name, testSet.depth, testSet.path)
            for test in testSet.tests:
                print("    -- ", os.path.relpath(test,testDir))

def init_tests(unitTestPath, logFilePath, tolerance=None, quiet=False, mpiLog=False, fullLog = False):
    """ Initialise the testing environment """
    global unitPath
    global testResults
    global tests
    global testSets

    unitPath = [path for path in unitTestPath.split(':') if path]
    unitPath += [testDir]

    for path in range(len(unitPath)):
        unitPath[path] = unitPath[path].rstrip('/ ')
    testResults.set_fulllog(fullLog)
    testResults.set_mpilog(mpiLog)
    testResults.open_log(logFilePath)
    testResults.set_tol(tolerance)
    testResults.set_quiet(quiet)
    testSets = construct_test_list(unitPath)

def construct_test_list(path):

    testSets = SuperTestSet()

    def find_tests(direc, depth):

        currDir = os.path.basename(direc)
        depth += 1
        tests = []

        # Search directory
        for obj in os.listdir(direc):

            loc = os.path.join(direc,obj) # Absolute path to current target file obj

            if obj.startswith('__') and obj.endswith('__'): continue # Skip dunders

            elif os.path.isfile(loc):

                filename, fileType = os.path.splitext(obj)
                if fileType in [".test",".py"]:
                    tests.append( loc )
                    testSets.append( TestSet(filename, depth + 1, direc, [ loc ]) )

                else: continue

            elif os.path.isdir(loc):

                inDir = os.path.basename(loc)
                tests += find_tests(loc, depth).tests

        testSet = TestSet( currDir, depth, direc, tests )
        testSets.append( testSet )

        return testSet

    depth = 0

    for direc in path:
        inDir = os.path.basename(direc)
        find_tests(direc, depth)

    return testSets

def list_test_sets():
    return set(testSets.names())

def list_all_tests():
    return testSets.paths()

def find_file(filename, returnPath=False, write=False, core=False, searchPath = None):
    if searchPath: pass
    elif core: searchPath = testDir
    else: searchPath = unitPath

    for path in searchPath:
        if write:
            return open(os.path.join(path,filename), 'w')
        for root, dirs, filenames in os.walk(path):
            if filename in filenames:
                if returnPath:
                    return os.path.join(root,filename)
                else:
                    return QuESTTestFile(os.path.join(root,filename), fullTest=False)
    else:
        testResults.log(fnfWarning.format(filename))
        return None

def finalise_tests():
    """ Finalise the testing environment """
    global unitPath
    global Env
    global testResults
    del unitPath
    destroyQuESTEnv(Env)
    del Env
    testResults.print_results()
    testResults.close_log()
    del testResults

class QuESTTestFile:
    """ Class containing test file information """

    def __init__(self,filename, fullTest=True):
        self.File = open(filename,'r')
        self.path = filename
        self.name = os.path.basename(filename) # Remove path
        self.nLine = 0
        self.nTests = 0
        if fullTest: self.testType = self._file_type()

    def init_testfile(self):
        temp = ''
        try:
            temp = self.readline()
            self.nTests = int(temp)
        except ValueError:
            raise IOError(fileWarning.format(message='Header of file :\n'+temp+"\n does not contain the number of tests",
                                             file=self.name, line=self.nLine))

    def __del__(self):
        self.File.close()

    def readline(self,retSkip=False):
        """ Reads a line from a test file """
        skip = []
        for line in self.File:
            self.nLine += 1
            if line.find('#') > -1:
                lineStrip = line[:line.find('#')].strip()
            else:
                lineStrip = line.strip()
            if lineStrip:
                if retSkip:
                    return lineStrip, skip
                else:
                    return lineStrip
            else:
                if line.lstrip('#').strip():
                    skip += [line.lstrip('#').strip()]
        raise IOError(fileWarning.format(message='Unexpected end of file',file=self.name, line=self.nLine))

    def parse_args(self, line):
        """ Split arguments, but maintain arrays and complex arrays as a block """
        line = self._remove_brackets(line)
        return line.split()

    def _remove_brackets(self, line):
        """ Remove all brackets and underscores from a given string (for parsing complex/arrays nicely) """
        remBrac = ''.maketrans('[{()}]_|><','          ','[{()}]_|><')
        return line.translate(remBrac)

    def read_state_vec(self, numQubits = 0, denMat=None):
        """ Read the expected state vector into a qubit state """

        stateVec = ""
        if denMat:
            nStates = 2**(int(numQubits)*2)
            loadStr = "C"
        else :
            nStates = 2**int(numQubits)
            loadStr = "c"

        for state in range(nStates): # Compare final with expected states
            try:
                stateVec += self.readline()+","
            except ValueError:
                raise IOError(fileWarning.format(message='Bad state line', file=self.name, line=self.nLine))
        stateVec = self._remove_brackets(stateVec).rstrip(',')
        QubitsOut = argQureg(numQubits, loadStr, testFile = self, initBits = stateVec, denMat = denMat)

        del stateVec
        return QubitsOut

    def _file_type(self):
        """ Determine testfile type """

        if self.name.endswith(".py"):
            return "Python"
        # Read first non-blank line
        line = self.File.readline().lstrip('# ').strip()
        while(not line):
            line = self.File.readline().lstrip('# ').strip()
        return line

class TestResults:
    """ Main class regarding testing framework stores results and comparisons """

    testLog = '{name}:{qureg}{string} {testing}'
    resultTab = '\nResult : {result}                    Expected: {expect}'

    def __init__(self, tolerance = 1.e-6, printToScreen = True, mpiLog=False, fullLog = False):
        self.passes, self.fails, self.numTests = [0]*3
        self.set_tol(tolerance)
        self.set_quiet(not printToScreen)
        self._logFilePath = None
        self._logFile = None
        self.set_fulllog(fullLog)
        self.set_mpilog(mpiLog)

    def _write_term(self, *out, **kwargs):
        """ Wrapper for print with enforced flushing and silencing """
        if self._printToScreen and root: print(*out, **kwargs, flush=True)

    def open_log(self, logFile):
        """ Open a new logFile associated with testResults """
        self.close_log()
        if Env.numRanks > 1 and self._mpiLog:
            self._logFile = open(logFile+".{}".format(Env.rank), 'w')
        else:
            self._logFile = open(logFile, 'w')
        self._logFilePath = logFile

    def close_log(self):
        """ Close the logFile associated with testResults """
        if self._logFile is not None:
            self._logFile.close()
        self._logFile = None

    def log(self, message = "", end = "\n"):
        """ Write a message to the log file (by default followed by new line) """
        if self._mpiLog:
            self._logFile.write(message+end)
        elif root:
            self._logFile.write(message+end)
        self._logFile.flush()

    def set_fulllog(self, fullLog):
        """ Set whether to log successes as well as failures """
        self._fullLog = fullLog

    def set_mpilog(self, mpiLog):
        """ Set whether to log for each process of the test results """
        self._mpiLog = mpiLog
                           
    def set_tol(self, tol = None):
        """ Set the tolerance of the test results """
        self._tolerance = tol
        self._tolOrder = -int(math.log10(tol))

    def set_quiet(self, quiet = True):
        """ Set the quiet status of the test results """
        self._printToScreen = not quiet

    def compareStates(self, a, b, tol = None):
        """ Check that all the elements of two state vectors are within tolerance of each other """

        if tol is None:
            tol = self._tolerance

        if a.isDensityMatrix and not b.isDensityMatrix or a.isDensityMatrix and not b.isDensityMatrix:
            raise TypeError('A and B are not both density matrices')

        if a.numQubitsRepresented != b.numQubitsRepresented:
            raise IndexError('A and B registers are not the same size')

        # Compare final with expected states
        if a.isDensityMatrix and b.isDensityMatrix:

            for row in range(a.numQubitsRepresented):
                for col in range(b.numQubitsRepresented):
                    aState = getDensityAmp(a,row,col)
                    bState = getDensityAmp(b,row,col)
                    if not self.compareComplex(aState,bState,tol): return False

        else:
            for state in range(getNumAmps(a)):
                aState = getAmp(a,state)
                bState = getAmp(b,state)
                if not self.compareComplex(aState,bState,tol): return False

        return True

    def compareReals(self, a, b, tol = None):
        """ Check that two real numbers are within tolerance of each other """
        if tol is None:
            tol = self._tolerance
        if abs(a - b) > tol: return False
        return True

    def compareComplex(self, a, b, tol = None):
        """ Check that two complex numbers are within tolerance of each other """
        if tol is None:
            tol = self._tolerance
        if abs(a.real - b.real) > tol or abs(a.imag - b.imag) > tol: return False
        return True

    def pass_test(self, testName=""):
        """ Force a pass to be logged """
        self._write_term('.',end='')
        if self._fullLog: self.log('Test {} Passed\n'.format(testName.strip()))
        self.numTests += 1
        self.passes += 1

    def fail_test(self, testName = "", message = ""):
        """ Force a fail test to be logged """
        self._write_term('F',end='')
        if testName or message:
            self.log('Test {} Failed: {}\n'.format(testName,message))
        self.numTests += 1
        self.fails += 1

    def validate(self, arg, test = "", message = ""):
        """ Call corresponding pass/fail function depending on value of arg """
        if arg:
            self.pass_test(test)
        else:
            self.fail_test(test, message)

    def print_results(self):
        """ Print number of passes and fails """
        self._write_term('\nPassed {} of {} tests, {} failed.\n'.format(self.passes,self.numTests,self.fails))

    def _run_test(self, testFunc, testFile):
        """ Read a test file and run corresponding test if in standard .test format """
        qubitTypeNames = {
            "Z":"Zero ", "C":"Custom ", "B":"BitState ", "P":"Plus ", "D":"Debug ",
            "z":"Zero ", "c":"Custom ", "b":"BitState ", "p":"Plus ", "d":"Debug "
        }

        for test in range(testFile.nTests):
            line, testComment = testFile.readline(True)
            testString,nBits,*args = testFile.parse_args(line)
            qubitType, *testType = testString.split('-')

            # Skip non-tests
            if int(nBits) == 0: continue
            
            denMat = qubitType.isupper()

            bitString = ""
            if qubitType in "CBcb":
                bitString = args[0]
                del args[0]
                Qubits = argQureg(nBits, qubitType, testFile, initBits = bitString, denMat = denMat)
            else:
                Qubits = argQureg(nBits, qubitType, testFile, denMat = denMat)

            args.insert(0,Qubits)

            retType = testFunc.thisFunc.restype
            if retType is None:
                testFunc(args)
                if not testType: testType = ["S"]
                for test in testType[0]:
                    if test in "Mm":
                        expectState = [None]*Qubits.numQubitsRepresented
                        success = True
                        for qubit in range(Qubits.numQubitsRepresented):
                            expectState[qubit] = list(map(float,testFile.readline().split()))
                            success = (success and
                                       self.compareReals(calcProbOfOutcome(Qubits, qubit, 0), expectState[qubit][0]) and
                                       self.compareReals(calcProbOfOutcome(Qubits, qubit, 1), expectState[qubit][1]))
                        self.validate(success,
                                      self.testLog.format(name = testFunc.funcname,
                                                          qureg = qubitTypeNames[qubitType],
                                                          string = bitString,
                                                          testing = "Measure"),
                                      self.resultTab.format(expect = "", result = "")
                        )
                        if not success:
                            for qubit in range(Qubits.numQubitsRepresented):
                                self.log('{} {}'.format([calcProbOfOutcome(Qubits, qubit, 0),calcProbOfOutcome(Qubits, qubit, 1)],expectState[qubit]))

                    elif test in "Ss":
                        expectState = testFile.read_state_vec(nBits,denMat = denMat)
                        success = self.compareStates(Qubits, expectState)
                        self.validate(success,
                                      self.testLog.format(name = testFunc.funcname,
                                                          qureg = qubitTypeNames[qubitType],
                                                          string = bitString,
                                                          testing = "State Vec"),
                                      self.resultTab.format(expect = "", result = "")
                        )
                        if not success: # Print resultant state vectors
                            if not Qubits.isDensityMatrix:
                                for state in range(getNumAmps(Qubits)):
                                    a = getAmp(Qubits, state)
                                    b = getAmp(expectState, state)
                                    self.log('{} {}'.format(a, b))
                            else:
                                for row in range(2**Qubits.numQubitsRepresented):
                                    for col in range(2**Qubits.numQubitsRepresented):
                                        a = getDensityAmp(Qubits, row, col)
                                        b = getDensityAmp(expectState, row, col)
                                        self.log('{} {}'.format(a, b))

                    elif test in "Pp":
                        expectState = float(testFile.readline())
                        success = self.compareReals(calcTotalProb(Qubits),expectState)

                        self.validate(success,
                                      self.testLog.format(name = testFunc.funcname,
                                                          qureg = qubitTypeNames[qubitType],
                                                          string = bitString,
                                                          testing = "Total Probability"),
                                      self.resultTab.format( expect = expectState, result = calcTotalProb(Qubits) )
                        )

                    else:
                        raise IOError('Unrecognised test type')
                    del expectState

            else: # Returning functions
                result = testFunc(args)

                if retType is Complex:
                    expect = argComplex(testFile.readline())
                    success = self.compareComplex(result,expect)
                elif retType is qreal:
                    expect = float(testFile.readline())
                    success = self.compareReals(result,expect)
                elif retType is c_int:
                    expect = int(testFile.readline())
                    success = expect == result

                else:
                    raise TypeError('Cannot test type {} currently'.format(retType.__name__))

                self.validate(success,
                              self.testLog.format(name = testFunc.funcname,
                                                  qureg = qubitTypeNames[qubitType],
                                                  string = bitString,
                                                  testing = "Return Val"),
                              self.resultTab.format( result = result, expect = expect )
                )

            destroyQureg(Qubits,Env)

        del testFile

    def run_test(self, testFileName, core=False):
        """ Run a test which is not listed in standard test sets """
        if os.path.isfile(testFileName):
            testPath = testFileName
        else:
            for path in unitPath:
                testPath = path+testFileName+'.test'
                if os.path.isfile(testPath) : break
            else:
                self.log(fnfWarning.format(testPath))
                self.fail_test()
                return

        if not core: self._write_term('Running test '+testPath+":", end=' ')

        testFunc = QuESTTestFile(testPath).testType

        if testFunc.capitalize() == "Python": # If file flagged as Python
            self._run_python_test(testPath)
            return

        elif testFunc not in list_funcnames():
            raise IOError(funWarning.format(testFunc))

        testFunc = QuESTTestee.get_func(testFunc)
        testFile = QuESTTestFile(testPath)
        testFile.init_testfile()

        self._run_test(testFunc, testFile)

        if not core: self._write_term()


    def _run_python_test(self, testPath, generate = False, nQubits = None):
        """ Internal function to import and run a Python style test """
        spec = importlib.util.spec_from_file_location("templib", testPath)
        templib = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(templib)
        if generate and hasattr(templib,'gen_tests'):
            templib.gen_tests(nQubits)
            self._write_term(".", end="")
        elif generate and not hasattr(templib,'gen_tests'):
            self._write_term("S", end="")
            self.log('Unable to generate test for Python test {} no gen_tests'.format(os.path.basename(testPath)))
        elif not generate: templib.run_tests()

        del templib


    def run_tests(self,testsToRun):
        """ Run corresponding tests """

        for test in testsToRun:
            for testSet in testSets(test):
                if not testSet.tests: continue
                self._write_term('Running tests '+testSet.name+":", end=' ')
                path = testSet.path
                core = testDir in path
                for test in testSet.tests:
                    # try:
                    self.run_test(test, core=core)
                    # except Exception as err:
                    #     self.fail_test(testName = test, message = err)

                self._write_term()

    def _write_gen_results(self, outputFile, testGen, qubitOp, result = None, Qubits = None):
        if qubitOp:
            for test in testGen:
                if test in "Pp":
                    outputFile.write("{}\n".format(round(calcTotalProb(Qubits),self._tolOrder+2)))
                elif test in "Ss":
                    for elem in Qubits._state_vec() : outputFile.write(str(elem))
                elif test in "Mm":
                    for qubit in range(Qubits.numQubitsRepresented):
                        outputFile.write("{} {}\n".format( round(calcProbOfOutcome(Qubits, qubit, 0),self._tolOrder+2),
                                                           round(calcProbOfOutcome(Qubits, qubit, 1),self._tolOrder+2) ))
                else:
                    raise IOError('Test type {} not recognised'.format(test))
        elif result:
            outputFile.write("{}\n".format(result))
        else:
            raise IOError('Unknown return in write_gen_results')


   
    def gen_std_test(self,testFunc, testFile, nQubits = 3, qubitGen = "zpdnZPDR", testGen = "PMS", targScan = "E", contScan = "E"):
        """ Generate individual test for a given function """

        for i in range(1,testFunc.nArgs):
            if testFunc.defArg[i] is None:
                self._write_term('Unable to generate test for function {} invalid default arguments'.format(testFunc.funcname))
                return


        def skip_gen(reason):
            outputFile.write("# {}\n".format(reason))
            outputFile.write("C- 0\n")

        qubitOp = testFunc.thisFunc.restype is None

        # Expand out
        while qubitGen.find('*') > -1:
            i = qubitGen.find('*')
            dup = qubitGen[i-1]
            mul = ""
            j = 1
            while i+j < len(qubitGen) and qubitGen[i+j].isdecimal():
                mul += qubitGen[i+j]
                j += 1
            if j == 1: raise IOError('Repeater not followed by repeat value')
            qubitGen = qubitGen.replace(qubitGen[i-1:i+j], dup*int(mul))

        if testFunc.target:
            if testFunc.targetType == "Qubit":
                if targScan == "E": targScan = list(range(nQubits))
                elif targScan == "D": targScan = testFunc.defArg[testFunc.target]
                else: targScan = map(int, targScan.split(","))
            elif testFunc.targetType == "Index":
                if targScan == "E": targScan = list(range(2**nQubits))
                elif targScan == "D": targScan = testFunc.defArg[testFunc.target]
                else: targScan = map(int, targScan.split(","))
        else:
            targScan = [0]


            
        nTests = len(qubitGen) * len(targScan)

        niceNames = {
            "z": "Zero State Vector", "p": "Plus State Vector", "d": "Debug State Vector",
            "r": "Random State Vector", "n": "Normalised Random State Vector",
            "Z": "Zero Density Matrix", "P": "Plus Density Matrix", "D": "Debug Density Matrix",
            "R": "Random Density Matrix"
        }

        with open(testFile,'w') as outputFile:

            outputFile.write('# {}\n'.format(testFunc.funcname))
            # Number of tests to run
            outputFile.write('{}\n'.format(nTests))

            for target in targScan:
                for qubitType in qubitGen:
                    denMat = qubitType.isupper()

                    # Write test comment header
                    outputFile.write("\n# {}\n".format(niceNames[qubitType]))

                    # If sensible otherwise skip test
                    if ((testFunc.denMat is False and denMat is True) or
                        (testFunc.denMat is True  and denMat is False) ):
                        skip_gen("Not valid for this function")
                        continue

                    args = [argQureg(nQubits, qubitType, denMat=testFunc.denMat)]

                    # Build argslist
                    for arg in range(1,testFunc.nArgs): args.append(testFunc.defArg[arg])

                    # Overwrite target qubit
                    if testFunc.target:  args[testFunc.target] = target

                    if testFunc.control is not None:
                        if testFunc.controlType == "Single" and target == args[testFunc.control]:
                            skip_gen("Target same as control")
                            continue
                        elif testFunc.controlType == "Multi" and target in args[testFunc.control]:
                            skip_gen("Target same as control")
                            continue

                    
                    
                    # Build random qureg
                    if qubitType in "RNrn":

                        # If too large and not happy
                        if not args[0]._size_warn(maxElem=2**25):
                            skip_gen("Output too large, not printing")
                            continue

                        if denMat: outputFile.write("C-{} {} [".format(testGen,nQubits))
                        else:      outputFile.write("c-{} {} [".format(testGen,nQubits))

                        # Write state vec inline
                        for elem in args[0]._state_vec() : outputFile.write(str(elem).rstrip()+",")
                        outputFile.write("] ")
                    else: outputFile.write("{}-{} {}".format(qubitType,testGen,nQubits))

                    # Write final argslist
                    for arg in range(1,testFunc.nArgs): outputFile.write(" "+str(args[arg]))

                    outputFile.write("\n")

                    # Run function
                    result = testFunc(*args)
                    retType = testFunc.thisFunc.restype
                    if retType is None:
                        for test in testGen:
                            if test in "Pp":
                                outputFile.write("{}\n".format(round(calcTotalProb(args[0]),self._tolOrder+2)))
                            elif test in "Ss":
                                for elem in args[0]._state_vec() : outputFile.write(str(elem))
                            elif test in "Mm":
                                for qubit in range(args[0].numQubitsRepresented):
                                    probOut0 = round(calcProbOfOutcome(args[0], qubit, 0),self._tolOrder+2)
                                    probOut1 = round(calcProbOfOutcome(args[0], qubit, 1),self._tolOrder+2)
                                    outputFile.write("{} {}\n".format(probOut0,probOut1))
                            else:
                                raise IOError('Test type {} not recognised'.format(test))
                    else:
                        outputFile.write("{}\n".format(result))

        self._write_term(".", end="")


    def gen_cust_test(self, testFileName, nQubits = 3):
        """ Generate a test which is not listed in standard test sets """

        if os.path.isfile(testFileName):
            testPath = testFileName
        else:
            for path in unitPath:
                testPath = path+testFileName+'.test'
                if os.path.isfile(testPath) : break
            else:
                self.log(fnfWarning.format(testPath))
                return

        testType = QuESTTestFile(testPath).testType

        if testType.capitalize() == "Python": # If file flagged as Python
            self._run_python_test(testPath, generate = True, nQubits = nQubits)
        elif testType in list_funcnames():
            self._write_term('\n Cannot generate non-Python tests in this way')
            quit()
        else:
            raise IOError('Unrecognised filetype in generation of file {}'.format(testPath))


    def gen_tests(self, testsToGen=["all"], nQubits=None, qubitGen = None, testGen = "PMS", targScan = "E", contScan = "E"):
        """ Generate sets of tests and skip if listed in don't_generate """

        protected = ["essential", "calcFidelity", # , "destroyQuESTEnv", "createQuESTEnv"
                  "calcInnerProduct", "measure", "measureWithStats"]
        coreTests = []
        for test in testsToGen:
            self._write_term("Generating tests " + test + ": ", end="")
            for testSet in testSets(test):
                for testFunc in testSet.names():
                    try:
                        if testFunc in protected: continue
                        elif testFunc in list_funcnames():
                            toGen = QuESTTestee.get_func(testFunc)
                            self.gen_std_test(toGen, os.path.join(unitPath[0],toGen.funcname+".test"), nQubits,
                                              qubitGen, testGen, targScan, contScan)
                        else:
                            self.gen_cust_test(next(testSets(testFunc).tests()))
                    except Exception as err:
                        self._write_term("F", end = "")
                        self.log("Error while generating {}: {}".format(testFunc, err))
        self._write_term()


def argQureg(nBits, qubitType, testFile=None, initBits = None, denMat = None):
    """ Create a qubit register from a standard .test style input """
    nBits = int(nBits)


    qubitTypes = {
        "z":initZeroState,"p":initPlusState,"d":initDebugState,"c":setAmps,"b":setAmps,"r":setAmps,"n":setAmps, # State Vector
        "Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps,"B":setAmps,"R":setAmps              # Density Matrix
    }

    if qubitType not in qubitTypes:
        if testFile is not None:
            raise IOError(fileWarning.format(message = 'Expected qubit state ({}), received {}'.format(",".join(qubitTypes.keys()), qubitType),
                                             file = testFile.name, line=testFile.nLine))
        else:
            raise IOError(fileWarning.format(message = 'Expected qubit state ({}), received {}'.format(",".join(qubitTypes.keys()), qubitType),
                                             file = "Unknown", line="Unknown"))


    if denMat is None: denMat = qubitType.isupper()

    qubitType = qubitType.upper()

    # Initialise Qubits
    if denMat :
        Qubits = createDensityQureg(nBits, Env)
    else :
        Qubits = createQureg(nBits, Env)


    if qubitType == "B":
        try:
            state = int(initBits, 2)
        except TypeError:
            if testFile is not None:
                raise IOError(fileWarning.format(message = 'Expected qubit state, received {}'.format(state),
                                                 file = testFile.name, line=testFile.nLine))
            else:
                raise IOError(fileWarning.format(message = 'Expected qubit state, received {}'.format(state),
                                                 file = "Unknown", line="Unknown"))

        nIn = len(initBits)
        if (nBits != nIn):
            if testFile is not None:
                raise IOError(
                    fileWarning.format(message = 'Bad number of states expected {}, received {}'.format(nBits, nIn)),
                    file = testFile.name, line=testFile.nLine)
            else:
                raise IOError(fileWarning.format(message = 'Expected qubit state, received {}'.format(state),
                                                 file = "Unknown", line="Unknown"))

        initClassicalState(Qubits, state)

        del nIn

    elif qubitType == "C": # Handle custom initialisation
        nStates = Qubits.numAmpsTotal
        nReqStates = nStates*2 # Account for complexes
        initBits = stringToList(initBits)
        nInStates = len(initBits)
        if nReqStates != nInStates:
            raise IOError(
                fileWarning.format(message = 'Bad number of states expected {}, received {}'.format(nStates, nInStates),
                file = testFile.path, line=testFile.nLine))

        qAmpsReal, qAmpsImag = initBits[0::2], initBits[1::2] # Split into real and imag
        if not denMat: setAmps(Qubits, 0, qAmpsReal, qAmpsImag, nStates)
        else: setDensityAmps(Qubits, qAmpsReal, qAmpsImag)
        del nStates, nReqStates, nInStates, qAmpsReal, qAmpsImag

    elif qubitType == "R":
        # Random
        nStates = Qubits.numAmpsTotal

        qAmpsReal = []
        qAmpsImag = []
        for i in range(nStates):
            qAmpsReal += [random.random()]
            qAmpsImag += [random.random()]

        if not denMat: setAmps(Qubits, 0, qAmpsReal, qAmpsImag, nStates)
        else: setDensityAmps(Qubits, qAmpsReal, qAmpsImag)
        del nStates, qAmpsReal, qAmpsImag

    elif qubitType == "N":
        # Normalised Random
        nStates = Qubits.numAmpsTotal

        qAmpsReal = []
        qAmpsImag = []
        for i in range(nStates):
            qAmpsReal += [random.random()]
            qAmpsImag += [random.random()]

        tProb = 0
        for i in range(nStates):
            tProb += (qAmpsReal[i] ** 2) + (qAmpsImag[i] ** 2)

        if tProb < 0: tProb *= -1
        tProb = math.sqrt(tProb)
        for i in range(nStates):
            qAmpsReal[i] /= tProb
            qAmpsImag[i] /= tProb

        if not denMat: setAmps(Qubits, 0, qAmpsReal, qAmpsImag, nStates)
        else: setDensityAmps(Qubits, qAmpsReal, qAmpsImag)
        del nStates, qAmpsReal, qAmpsImag

    else:
        qubitTypes[qubitType](Qubits)

    return Qubits


def repEnv():
    reportQuESTEnv(Env)

# Make Key variables publically accessible
Env = createQuESTEnv()
root = Env.rank == 0
unitPath = None
testResults = TestResults()
testSets = None
testDir = os.path.normpath(os.path.join(os.path.dirname(__file__),"../../tests"))
