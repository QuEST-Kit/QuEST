import os.path
from QuESTFunc import *
from testset import tests
import importlib.util
import importlib.machinery


# Add .test to valid python names
importlib.machinery.SOURCE_SUFFIXES += ['.test']

def init_tests(unitTestPath, logFilePath, tolerance=None, quiet=False, fullLogging = False):
    """ Initialise the testing environment """
    global unitPath
    global testResults

    unitPath = unitTestPath.split(':')
    for path in range(len(unitPath)):
        unitPath[path] = unitPath[path].rstrip('/ ') + '/'

    testResults.set_fulllog(fullLogging)
    testResults.open_log(logFilePath)
    testResults.set_tol(tolerance)
    testResults.set_quiet(quiet)

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

    def __init__(self,filename):
        self.File = open(filename,'r')
        self.name = filename[filename.rfind('/')+1:] # Remove path
        self.nLine = 0
        self.nTests = 0

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

    def read_state_vec(self, numQubits = 0, denMat=False):
        """ Read the expected state vector into a qubit state """

        stateVec = ""
        if denMat: nStates = 2**(int(numQubits)*2)
        else : nStates = 2**int(numQubits)
        
        for state in range(nStates): # Compare final with expected states
            try:
                stateVec += self.readline()+","
            except ValueError:
                raise IOError(fileWarning.format(message='Bad state line', file=self.name, line=self.nLine))
        stateVec = self._remove_brackets(stateVec).rstrip(',')
        QubitsOut = argQureg(numQubits, 'C', testFile = self, initBits = stateVec, denMat = denMat)
        
        del stateVec
        return QubitsOut


class TestResults:
    """ Main class regarding testing framework stores results and comparisons """

    def __init__(self, tolerance = 1.e-6, printToScreen = True, fullLogging = False):
        self.passes, self.fails, self.numTests = [0]*3
        self._tolerance = tolerance
        self._printToScreen = printToScreen
        self._logFilePath = None
        self._logFile = None
        self._fullLog = fullLogging

    def _write_term(self, *out, **kwargs):
        """ Wrapper for print with enforced flushing and silencing """
        if self._printToScreen and root: print(*out, **kwargs, flush=True)

    def open_log(self, logFile):
        """ Open a new logFile associated with testResults """
        self.close_log()
        if Env.numRanks > 1 and self._fullLog:
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
        if self._fullLog:
            self._logFile.write(message+end)
        elif root:
            self._logFile.write(message+end)
        self._logFile.flush()

    def set_fulllog(self, fullLog = False):
        """ Set whether to log for each process of the test results """
        self._fullLog = fullLog

    def set_tol(self, tol = None):
        """ Set the tolerance of the test results """
        self._tolerance = tol

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
        self.log('{} Passed\n'.format(testName.strip()))
        self.numTests += 1
        self.passes += 1

    def fail_test(self, testName = "", message = ""):
        """ Force a fail test to be logged """
        self._write_term('F',end='')
        if testName or message:
            self.log('Test {} failed: {}\n'.format(testName,message))
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

        qubitTypeNames = {"Z":"Zero ", "C":"Custom ", "B":"BitState ", "P":"Plus ", "D":"Debug "}
        for test in range(testFile.nTests):
            line, testComment = testFile.readline(True)
            qubitType,nBits,*args = testFile.parse_args(line)

            bitString = ""
            if qubitType in "CBcb":
                bitString = args[0]
                del args[0]
                Qubits = argQureg(nBits, qubitType, testFile, initBits = bitString, denMat = testFunc.denMat)
            else:
                Qubits = argQureg(nBits, qubitType, testFile, denMat = testFunc.denMat)

            args.insert(0,Qubits)

            retType = testFunc.thisFunc.restype
            if retType is None:
                testFunc(args)
                expectState = testFile.read_state_vec(nBits,denMat = testFunc.denMat)
                success = testResults.compareStates(Qubits, expectState)
            else:
                result = testFunc(args)

                if retType is Complex:
                    expect = argComplex(testFile.readline())
                    success = self.compareComplex(result,expect)
                elif retType is c_double:
                    expect = float(testFile.readline())
                    success = self.compareReals(result,expect)
                elif retType is c_int:
                    expect = int(testFile.readline())
                    success = expect == result

                else:
                    raise TypeError('Cannot test type {} currently'.format(retType.__name__))

            if success:
                self.pass_test("{}{}".format(qubitTypeNames[qubitType],bitString))
            else:
                if testComment:
                    self.log('Test {Func}:{Comm} failed in {File}\n'.format(
                        Func = testFunc.funcname,
                        Comm = "\n".join(testComment),
                        File = testFile.name))
                else:
                    self.log('Testing {Func}:{Comm} failed in {File}\n'.format(
                        Func = testFunc.funcname,
                        Comm = qubitTypeNames[qubitType],
                        File = testFile.name))
                if retType is None:
                    if not Qubits.isDensityMatrix:
                        self.log('Receieved                               Expected')
                        for state in range(getNumAmps(Qubits)):
                            a = getAmp(Qubits, state)
                            b = getAmp(expectState, state)
                            self.log('{} {}'.format(a, b))
                    else:
                        self.log('Receieved                               Expected')
                        for row in range(Qubits.numQubitsRepresented):
                            for col in range(Qubits.numQubitsRepresented):
                                a = getDensityAmp(Qubits, row, col)
                                b = getDensityAmp(expectState, row, col)
                                self.log('{} {}'.format(a, b))
                        
                else:
                    self.log('{} {}'.format( result, expect))
                self.fail_test()

            destroyQureg(Qubits,Env)

        del testFile 

    def run_std_test(self, testFuncsList,name=''):
        """ Run set of tests in testsets.py """
        self._write_term('Running tests '+name+":", end=' ')

        for testFunc in testFuncsList:
            self.log('\nRunning test {}\n'.format(testFunc.funcname))

            for path in unitPath:
                testPath = path+testFunc.funcname+'.test'
                if os.path.isfile(testPath) : break
            else:
                self.log(fnfWarning.format(testPath))
                self.fail_test()
                continue

            with open(testPath,'r') as testFile:
                testPyth = testFile.readline().lstrip('# ').strip()

                if testPyth == "Python": # If file flagged as Python
                    self._run_python_test(testPath)
                    continue
            testFile = QuESTTestFile(testPath)
            testFile.init_testfile()
            self._run_test(testFunc, testFile)

        self._write_term()

    def _run_python_test(self, testPath, generate = False, nQubits = None):
        """ Internal function to import and run a Python style test """
        spec = importlib.util.spec_from_file_location("templib", testPath)
        templib = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(templib)
        if generate and templib.gen_tests is not None: templib.gen_tests(nQubits)
        else: templib.run_tests()
        del templib

    def run_cust_test(self, testFileName):
        """ Run a test which is not listed in standard test sets """
        if os.path.isfile(testFileName):
            testPath = testFileName
        else:
            for path in unitPath:
                testPath = path+testFileName+'.test'
                if os.path.isfile(testPath) : break
            else:
                self.log(fnfWarning.format(testPath))
                self._write_term('Running test '+testFileName+":", end=' ')
                self.fail_test()
                return

        self._write_term('Running test '+testPath+":", end=' ')

        if testPath.endswith('.test'): # Run standard formatted tests
            with open(testPath,'r') as testFile:
                testFunc = testFile.readline().lstrip('# ').strip()

                if testFunc.capitalize() == "Python": # If file flagged as Python
                    self._run_python_test(testPath)
                    return

                if testFunc not in list_funcnames():
                    raise IOError(funWarning.format(testFunc.funcname))
            testFunc = tests[testFunc]
            testFile = QuESTTestFile(testPath)
            testFile.init_testfile()

            self._run_test(*testFunc, testFile)


        elif testPath.endswith('.py'): # Run custom test scripts
            self._run_python_test(testPath)
        else:
            raise IOError('Unrecognised filetype in test run of file {}'.format(testPath))

        self._write_term()

    def run_tests(self,testsToRun):
        """ Run corresponding tests """
        from testset import testSets, tests

        for test in testsToRun:
            if test in testSets or test in tests:
                # try:
                self.run_std_test(tests[test],test)
                # except KeyError:
                #     print("Function '{}' does not exist, are you sure you wrote it correctly?".format(test))
            else:
                self.run_cust_test(test)

    def gen_std_test(self,testFunc, testFile, nQubits = 3):
        """ Generate individual test for a given function """

        for i in range(1,testFunc.nArgs):
            if testFunc.defArg[i] is None:
                print('Unable to generate test for function {} invalid default arguments'.format(testFunc.funcname))
                return

        with open(testFile,'w') as outputFile:

            outputFile.write('# {}\n'.format(testFunc.funcname))
            # Standard run 3 tests
            outputFile.write('3\n')

            for qubitType in "ZPD":
                args = [argQureg(nQubits, qubitType, denMat=testFunc.denMat)]
                argString = f"{qubitType} {nQubits}"
                for arg in range(1,testFunc.nArgs):
                    args += [testFunc.defArg[arg]]
                    argString += " "+str(testFunc.defArg[arg])
                outputFile.write(argString+"\n")
                result = testFunc(*args)
                retType = testFunc.thisFunc.restype
                if retType is None:
                    for elem in args[0]._state_vec() : outputFile.write(str(elem))
                else:
                    outputFile.write(f"{result}\n")

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

        if testPath.endswith('.test'): # Run standard formatted tests
            with open(testPath,'r') as testFile:
                testFunc = testFile.readline().lstrip('# ').strip()

            if testFunc.capitalize() == "Python": # If file flagged as Python
                self._run_python_test(testPath, generate = True, nQubits = nQubits)
            else:
                self.log('Cannot generate not python tests in this way')

        elif testPath.endswith('.py'): # Run custom test scripts
            self._run_python_test(testPath, generate = True)
        else:
            raise IOError('Unrecognised filetype in test run of file {}'.format(testPath))

    def gen_tests(self, testsToGen=["all"], nQubits=None):
        """ Generate sets of tests and skip if listed in don't_generate """

        # Test if enough memory exists
        temp = createQureg(nQubits,Env)
        if not temp._size_warn():
            destroyQureg(temp,Env)
            return
        destroyQureg(temp,Env)

        from testset import tests
        for testSet in testsToGen:
            if testSet in tests:
                for testFunc in tests[testSet] :
                    if testFunc in tests["don't_generate"]: continue
                    if testFunc in tests["pyth_generate"]:
                        self.gen_cust_test(testSet)
                    else:
                        self.gen_std_test(testFunc, unitPath[0]+testFunc.funcname+".test", nQubits)
            else:
                self.gen_cust_test(testSet)


def argQureg(nBits, qubitType, testFile=None, initBits = None, denMat = False):
    """ Create a qubit register from a standard .test style input """
    nBits = int(nBits)

    #Upcase qubitType
    qubitType = qubitType.upper()
    # Initialise Qubits
    if denMat :
        Qubits = createDensityQureg(nBits, Env)
    else :
        Qubits = createQureg(nBits, Env)

    qubitTypes = {"Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps,"B":setAmps}

    if qubitType not in qubitTypes:
        if testFile is not None:
            raise IOError(fileWarning.format(message = 'Expected qubit state ({}), received {}'.format(",".join(qubitTypes.keys()), qubitType),
                                             file = testFile.name, line=testFile.nLine))
        else:
            raise IOError(fileWarning.format(message = 'Expected qubit state ({}), received {}'.format(",".join(qubitTypes.keys()), qubitType),
                                             file = "Unknown", line="Unknown"))

    elif qubitType == "B":
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
                file = testFile.name, line=testFile.nLine))

        qAmpsReal, qAmpsImag = initBits[0::2], initBits[1::2] # Split into real and imag
        if not denMat: setAmps(Qubits, 0, qAmpsReal, qAmpsImag, nStates)
        else: setDensityAmps(Qubits, qAmpsReal, qAmpsImag)
        del nStates, nReqStates, nInStates, qAmpsReal, qAmpsImag

    else:
        qubitTypes[qubitType](Qubits)

    return Qubits

# Make Key variables publically accessible
Env = createQuESTEnv()
root = Env.rank == 0
unitPath = None
testResults = TestResults()
