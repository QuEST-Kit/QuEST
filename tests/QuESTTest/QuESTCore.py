import os.path
import os
from QuESTPy.QuESTFunc import *
#from .testset import tests
import importlib.util
import importlib.machinery

# Add .test to valid python names
importlib.machinery.SOURCE_SUFFIXES += ['.test']

def init_tests(unitTestPath, logFilePath, tolerance=None, quiet=False, fullLogging = False):
    """ Initialise the testing environment """
    global unitPath
    global testResults
    global tests
    global testSets
    
    currDir = os.path.dirname(__file__)
    coreDirs = [f'{currDir}/{directory}' for directory in ['essential','algor','benchmarks','unit']]

    if unitTestPath: unitPath = unitTestPath.split(':') + coreDirs
    else: unitPath = coreDirs

    for path in range(len(unitPath)):
        unitPath[path] = unitPath[path].rstrip('/ ')
        
    testResults.set_fulllog(fullLogging)
    testResults.open_log(logFilePath)
    testResults.set_tol(tolerance)
    testResults.set_quiet(quiet)
    tests, testSets = construct_test_list(unitPath)

def construct_test_list(path):
    def find_tests(direc, depth):
        currDir = os.path.basename(direc)
        depth += 1
        testSets.add( (currDir, depth) )
        tests[currDir] = []
        if depth not in tests: tests[depth] = []
        for obj in os.listdir(direc):
            loc = os.path.join(direc,obj)
            target = (direc, obj)
            if obj.startswith('__') and obj.endswith('__'): continue # Skip dunders
            elif os.path.isfile(loc):
                filename, fileType = os.path.splitext(obj)
                if fileType in [".test",".py"]:
                    tests[currDir].append( target )
                    tests[filename] = [ target ]
                else:
                    continue
            elif os.path.isdir(loc):
                inDir = os.path.basename(loc)
                find_tests(loc, depth)
                tests[currDir] += tests[inDir]
                
    depth = 0
    tests = {}
    testSets = set( [("all", 0)] )
    tests["all"] = []
    for direc in path:
        inDir = os.path.basename(direc)
        find_tests(direc, depth)
        tests["all"] += tests[inDir]

    return tests, testSets

def list_test_sets():
    return sorted(testSets, key=lambda x: x[1])

def list_all_tests():
    return (key for key in tests.keys() if not isinstance(key,int))

def find_file(filename, write=False):
    for path in unitPath:
        if write:
            return open(os.path.join(path,filename), 'w')
        for root, dirs, filenames in os.walk(path):
            if filename in filenames:
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
        self.name = os.path.basename(filename) #filename[filename.rfind('/')+1:] # Remove path
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

    def __init__(self, tolerance = 1.e-6, printToScreen = True, fullLogging = False):
        self.passes, self.fails, self.numTests = [0]*3
        self.set_tol(tolerance)
        self.set_quiet(not printToScreen)
        self._logFilePath = None
        self._logFile = None
        self.set_fulllog(fullLogging)

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
            testString,nBits,*args = testFile.parse_args(line)
            qubitType, *testType = testString.split('-') 
            
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
                if not testType: testType = ["S"]
                for test in testType[0]:
                    if test in "Mm":
                        expectState = [None]*Qubits.numQubitsRepresented
                        success = True
                        for qubit in range(Qubits.numQubitsRepresented):
                            expectState[qubit] = list(map(float,testFile.readline().split()))
                            success = (success and
                                       testResults.compareReals(calcProbOfOutcome(Qubits, qubit, 0), expectState[qubit][0]) and
                                       testResults.compareReals(calcProbOfOutcome(Qubits, qubit, 1), expectState[qubit][1]))
                        testResults.validate(success,
                                             f'{testFunc.funcname}:{qubitTypeNames[qubitType]}{bitString} Measure',
                                             f'\nResult :                     Expected:')
                        for qubit in range(Qubits.numQubitsRepresented):
                            self.log(f'{[calcProbOfOutcome(Qubits, qubit, 0),calcProbOfOutcome(Qubits, qubit, 1)]} {expectState[qubit]}')
                            
                    elif test in "Ss":
                        expectState = testFile.read_state_vec(nBits,denMat = testFunc.denMat)
                        success = testResults.compareStates(Qubits, expectState)
                        testResults.validate(success,
                                             f"{testFunc.funcname}:{qubitTypeNames[qubitType]}{bitString}",
                                             "\nResult :                         Expected:  ")
                        if not success: # Print resultant state vectors
                            if not Qubits.isDensityMatrix:
                                for state in range(getNumAmps(Qubits)):
                                    a = getAmp(Qubits, state)
                                    b = getAmp(expectState, state)
                                    self.log('{} {}'.format(a, b))
                            else:
                                for row in range(Qubits.numQubitsRepresented):
                                    for col in range(Qubits.numQubitsRepresented):
                                        a = getDensityAmp(Qubits, row, col)
                                        b = getDensityAmp(expectState, row, col)
                                        self.log('{} {}'.format(a, b))

                    elif test in "Pp":
                        expectState = float(testFile.readline())
                        success = self.compareReals(calcTotalProb(Qubits),expectState)
                        
                        testResults.validate(success,
                                             f"{testFunc.funcname}:{qubitTypeNames[qubitType]}{bitString} Total Probability",
                                             f"Expected: {expectState}, Result: {calcTotalProb(Qubits)}")

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

                testResults.validate(success,
                                     f"{testFunc.funcname}:{qubitTypeNames[qubitType]}{bitString}",
                                     f"Expected: {expect}, Result: {result}")


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
                self._write_term('Running test '+testFileName+":", end=' ')
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
        if generate and hasattr(templib,'gen_tests'): templib.gen_tests(nQubits)
        elif generate and not hasattr(templib,'gen_tests'):
            print('Unable to generate test for Python test {} no gen_tests'.format(os.path.basename(testPath)))
        elif not generate: templib.run_tests()
        
        del templib


    def run_tests(self,testsToRun):
        """ Run corresponding tests """

        for test in testsToRun:
            if test in tests:
                self._write_term('Running test '+test+":", end=' ')
                for testee in tests[test]:
                    self.run_test( os.path.join(*testee), core = True )
                self._write_term()
            else:
                self.run_test(test)

    def gen_std_test(self,testFunc, testFile, nQubits = 3, qubitGen = "ZPDN", testGen = "PMS"):
        """ Generate individual test for a given function """

        for i in range(1,testFunc.nArgs):
            if testFunc.defArg[i] is None:
                print('Unable to generate test for function {} invalid default arguments'.format(testFunc.funcname))
                return

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
            
            
        niceNames = {"Z": "Zero State", "P": "Plus State", "D": "Debug State", "R": "Random State", "N": "Normalised Random State"}
            
        with open(testFile,'w') as outputFile:

            outputFile.write('# {}\n'.format(testFunc.funcname))
            # Standard run 5 tests
            outputFile.write(f'{len(qubitGen)}\n')

            for qubitType in qubitGen:
                outputFile.write(f"\n# {niceNames[qubitType]}\n")
                args = [argQureg(nQubits, qubitType, denMat=testFunc.denMat)]
                if qubitType in "RN":
                    if not args[0]._size_warn():
                        outputFile.write(f"C- 0\n")
                        continue
                    outputFile.write(f"C-{testGen} {nQubits} [")
                    for state in range(args[0].numAmpsTotal-1):
                        outputFile.write(str(args[0][state]) + ",")
                    outputFile.write(str(args[0][args[0].numAmpsTotal-1])+"] ")
                else: outputFile.write(f"{qubitType}-{testGen} {nQubits}")
                for arg in range(1,testFunc.nArgs):
                    args += [testFunc.defArg[arg]]
                    outputFile.write(" "+str(testFunc.defArg[arg]))
                outputFile.write("\n")
                result = testFunc(*args)
                retType = testFunc.thisFunc.restype
                if retType is None:
                    for test in testGen:
                        if test in "Pp":
                            outputFile.write(f"{round(calcTotalProb(args[0]),self._tolOrder+2)}\n")
                        elif test in "Ss":
                            for elem in args[0]._state_vec() : outputFile.write(str(elem))
                        elif test in "Mm":
                            for qubit in range(args[0].numQubitsRepresented):
                                outputFile.write(f"{round(calcProbOfOutcome(args[0], qubit, 0),self._tolOrder+2)} {round(calcProbOfOutcome(args[0], qubit, 1),self._tolOrder+2)}\n")
                        else:
                            raise IOError(f'Test type {test} not recognised')
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

        testType = QuESTTestFile(testPath).testType

        if testType.capitalize() == "Python": # If file flagged as Python
            self._run_python_test(testPath, generate = True, nQubits = nQubits)
        elif testType in list_funcnames():
            self.log('Cannot generate non-Python tests in this way')
        else:
            raise IOError('Unrecognised filetype in generation of file {}'.format(testPath))

    def gen_tests(self, testsToGen=["all"], nQubits=None, qubitGen = "ZPDN", testGen = "PMS"):
        """ Generate sets of tests and skip if listed in don't_generate """

        # Test if enough memory exists
        temp = createQureg(nQubits,Env)
        if not temp._size_warn():
            destroyQureg(temp,Env)
            return
        destroyQureg(temp,Env)

        if __package__ in unitPath[0].split(os.sep):
            self._write_term('Will not overwrite core tests, specify output directory with -p')
            return

        protected = ["essential", "destroyQuESTEnv", "createQuESTEnv", "calcFidelity",
                  "calcInnerProduct", "measure", "measureWithStats"]
        coreTests = []
        for i in protected:
            coreTests += [os.path.splitext(j[1])[0] for j in tests[i]] # Strip to function names
        for testSet in testsToGen:
            if testSet in tests:
                for path, testFunc in tests[testSet]:
                    toGen = os.path.splitext(testFunc)[0]
                    fullPath = os.path.join(path, testFunc)
                    if testFunc in coreTests: continue
                    if QuESTTestFile(fullPath).testType == "Python":
                        self.gen_cust_test(fullPath)
                    else:
                        toGen = QuESTTestee.get_func(toGen)
                        self.gen_std_test(toGen, os.path.join(unitPath[0],toGen.funcname+".test"), nQubits, qubitGen, testGen)
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

    qubitTypes = {"Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps,"B":setAmps,"R":setAmps,"N":setAmps}

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

# Make Key variables publically accessible
Env = createQuESTEnv()
root = Env.rank == 0
unitPath = None
testResults = TestResults()
tests = None
testSets = None
