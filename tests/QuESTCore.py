import os.path
from QuESTFunc import *
from QuESTTypes import *
import importlib


logFile = None
unitPath = None
Env = None
testResults = None
def init_tests(unitTestPath, logFilePath, tolerance=None):
    global logFile
    global unitPath
    global Env
    global testResults
    unitPath = unitTestPath
    logFile = open(logFilePath,'w')
    Env = createQuESTEnv()
    testResults = TestResults(tolerance)
    return testResults

def finalise_tests():
    global unitPath
    global logFile
    global Env
    global testResults
    del unitPath
    logFile.close()
    del logFile
    destroyQuESTEnv(Env)
    del Env
    del testResults
 
class QuESTTestFile:
    """ Class containing test file information """
    
    def __init__(self,filename):
        self.File = open(filename,'r')
        
        self.name = filename[filename.rfind('/')+1:] # Remove path
        self.nLine= 0
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
            lineStrip = line[:line.find('#')].strip()
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
        line = self.remove_brackets(line)
        return line.split()

    def remove_brackets(self, line):
        remBrac = ''.maketrans('[{()}]','      ','[{()}]')
        return line.translate(remBrac)
                
    def read_state_vec(self, numQubits = 0):
        """ Read the expected state vector into a qubit state """
        QubitsOut = createQureg(numQubits, Env)

        for state in range(getNumAmps(QubitsOut)): # Compare final with expected states
            try:
                stateElem = argComplex(self.readline())
            except ValueError:
                raise IOError(fileWarning.format(message='Bad state line', file=self.name, line=self.nLine))
            setAmps(QubitsOut, state, qreal(stateElem.real), qreal(stateElem.imag), 1)
        return QubitsOut
        
class TestResults:
    def __init__(self, tolerance = 1.e-6):
        self.passes, self.fails, self.numTests = [0]*3
        self.tolerance = tolerance
        
    def compareStates(self, a, b, tol = None):
        if tol is None:
            tol = self.tolerance

        if (a.numQubitsRepresented != b.numQubitsRepresented):
            raise IndexError('A and B registers are not the same size')
    
        for state in range(getNumAmps(a)): # Compare final with expected states
            aState = getAmp(a,state)
            bState = getAmp(b,state)
            if abs(aState.real - bState.real) > tol or abs(aState.imag - bState.imag) > tol: return False
        return True

    def compareReals(self, a, b, tol = None):
        if tol is None:
            tol = self.tolerance
        if abs(a - b) > tol: return False
        return True

    def compareComplex(self, a, b, tol = None):
        if tol is None:
            tol = self.tolerance
        if abs(aState.real - bState.real) > tol or abs(aState.imag - bState.imag) > tol: return False
        return True
                
    def pass_test(self):
        print('.',end='')
        self.numTests += 1
        self.passes += 1

    def fail_test(self):
        print('F',end='')
        self.numTests += 1
        self.fails += 1

    def print_results(self):
        print('\nPassed {} of {} tests, {} failed.\n'.format(self.passes,self.numTests,self.fails))
        
    def run_test(self, testFunc, testFile):
        for test in range(testFile.nTests):
            line, skip = testFile.readline(True)
            
            qubitType,nBits,*args = testFile.parse_args(line)
    
            if qubitType in "CBcb":
                Qubits = argQureg(nBits, qubitType, testFile, initBits = args[0])
                del args[0]
            else:
                Qubits = argQureg(nBits, qubitType, testFile)
                
            args.insert(0,Qubits)
    
            testFunc(args)
        
            if testFunc.thisFunc.restype is None:
                expectState = testFile.read_state_vec(nBits)
                success = testResults.compareStates(Qubits, expectState)
            else:
                pass
    
            if success:
                self.pass_test()
            else:
                for state in range(getNumAmps(Qubits)):
                    logFile.write('Test {} Failed \n'.format(''))
                    a = getAmp(Qubits, state)
                    b = getAmp(expectState, state)
                    logFile.write('{} {}\n'.format(a, b))
                self.fail_test()
                
            destroyQureg(Qubits,Env)

        del testFile

    def run_std_test(self, testFuncsList,name=''):
    
        print('Running tests '+name+":", end=' ')
        
        for testFunc in testFuncsList:
    
            logFile.write('\nRunning test {}\n'.format(testFunc.funcname))
            
            testPath = unitPath+testFunc.funcname+'.test'
            
            if os.path.isfile(testPath) :
                testFile = QuESTTestFile(testPath)
            else:
                logFile.write(fnfWarning.format(testPath))
                self.fail_test()
                continue
    
            self.run_test(testFunc, testFile)
    
        print()
    
    
    def run_cust_test(self, testFileName):
    
        if os.path.isfile(testFileName):
            testPath = testFileName
        else:
            testPath = unitPath+testFileName+'.test'
            
        print('Running test '+testPath+":", end=' ')
    
        if testPath.endswith('.test'): # Run standard formatted tests
            if os.path.isfile(testPath) :
                with open(testPath,'r') as testFile:
                    testFunc = testFile.readline().lstrip('# ').strip()
                    if testFunc not in list_funcnames():
                        raise IOError(funWarning.format(testFunc.funcname))
                testFile = QuESTTestFile(testPath)
    
                self.run_test(testFunc, testFile)
                  
            else:
                logFile.write(fnfWarning.format(testPath))
                self.fail_test()
    
    
        elif testPath.endswith('.py'): # Run custom test scripts
            testPath, ext = os.path.splitext(os.path.split(testPath)[-1])
            templib = importlib.import_module(testPath)
            templib.run_tests()
            del templib
        else:
            raise IOError('Unrecognised filetype in test run of file {}'.format(testPath))
                
        print()


def argQureg(nBits, qubitType, testFile=None, initBits = None):
    nBits = int(nBits)
        
    #Upcase qubitType
    qubitType = qubitType.upper()
    # Initialise Qubits
    Qubits = createQureg(nBits, Env)

    qubitTypes = {"Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps,"B":setAmps}
    
    if qubitType not in qubitTypes:
        raise IOError(fileWarning.format(message = 'Unrecognised qubit initialisation state "'+qubitType+'"',
                                  file = testFile.name, line=testFile.nLine))

    elif qubitType == "B":
        if any(bit not in "01" for bit in initBits ):
            raise IOError(fileWarning.format(message = 'Expected qubit state, received {}'.format(state)))

        try:
            state = int(initBits, 2)
        except TypeError:
            print(initBits)
            raise IOError(fileWarning.format(message = 'Expected qubit state, received {}'.format(state)))
        
        nIn = len(initBits)
        if (nBits != len(initBits)):
            raise IOError(
                fileWarning.format(message = 'Bad number of states expected {}, received {}'.format(nBits, nIn)),
                file = testFile.name, line=testFile.nLine)

        # Entirely zero state
        initZeroState(Qubits)
        setAmps(Qubits, 0, byref(qreal(0.0)), byref(qreal(0.0)), 1)

        # Set selected state
        setAmps(Qubits, state, byref(qreal(1.0)), byref(qreal(0.0)), 1)
        
        del nIn
        
    
    elif qubitType == "C": # Handle custom initialisation
        nStates = getNumAmps(Qubits)
        initBits = stringToList(initBits)
        nIn = len(initBits)
        if (nStates != len(initBits)):
            raise IOError(
                fileWarning.format(message = 'Bad number of states expected {}, received {}'.format(nStates, nIn)),
                file = testFile.name, line=testFile.nLine)
        
        qAmpsReal, qAmpsImag = initBits[0::2], initBits[1::2] # Split into real and imag
    
        setAmps(Qubits, 0, byref(qAmpsReal), byref(qAmpsImag), nStates)
        
        del nStates, nIn, qAmps, qAmpsReal, qAmpsImag
    
    else:
        qubitTypes[qubitType](Qubits)

    return Qubits


def gen_test(testFunc, testFile):
    for i in range(1,testFunc.nArgs):
        if testFunc.defArg[i] is None:
            print('Unable to generate test for function {} invalid default arguments'.format(testFunc.funcname))
            return

    with open(testFile,'w') as outputFile:
        
        outputFile.write('# {}\n'.format(testFunc.funcname))
        # Standard run 3 tests
        outputFile.write('3\n')
        
        for qubitType in "ZPD":
            nQubits = 2
            args = [argQureg(nQubits, qubitType)]
            argString = "{} {}".format(qubitType, nQubits)
            for arg in range(1,testFunc.nArgs):
                args += [testFunc.defArg[arg]]
                argString += " "+str(testFunc.defArg[arg])
            outputFile.write(argString+"\n")
            testFunc(*args)
            outputFile.write(str(args[0])+"\n")


def gen_tests():
    from testset import tests
    for testSet in ["stnd_operations","cont_operations"]: #["math_operations", "mcon_operations"]:
        for testFunc in tests[testSet] :
            gen_test(testFunc, unitPath+testFunc.funcname+".test")
