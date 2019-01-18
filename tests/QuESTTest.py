from ctypes import *
from QuESTTypes import *
from QuESTFunc import *
from sys import argv
import os.path
import importlib

logFile = open('QuESTTest.log','w')
tol = 1e-6
qubitTypes = {"Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps,"B":setAmps}

unitPath = "unitPy/"


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
        remBrac = ''.maketrans('[{()}]','      ','[{()}]')
        line = line.translate(remBrac)
        return line.split()

        
def run_test(testFunc, testFile):

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
    
        success = True
        if testFunc.thisFunc.restype is None:
            for state in range(getNumAmps(Qubits)): # Compare final with expected states
                resultState = getAmp(Qubits,state)
                expectState = Complex(*list(map(float,testFile.readline().split(','))))
                if abs(resultState.real - expectState.real) > tol or abs(resultState.imag - expectState.imag) > tol:
                    logFile.write('{:s}\n{:s}\ndiff {:s}\n'.format(
                        resultState, expectState, expectState - resultState))
                          
                    success = False
                    break
        else:
            pass

        if success:
            print('.',end='')
        else:
            print('F',end='')
            
        destroyQureg(Qubits,Env)
        
    del testFile

def argQureg(nBits, qubitType, testFile=None, initBits = None):
    nBits = int(nBits)
        
    #Upcase qubitType
    qubitType = qubitType.upper()
    # Initialise Qubits
    Qubits = createQureg(nBits, Env)

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


def run_std_test(testFuncsList,name=''):

    print('Running tests '+name+":", end=' ')
    
    for testFunc in testFuncsList:

        logFile.write('\nRunning test {}\n'.format(testFunc.funcname))
        
        testPath = unitPath+testFunc.funcname+'.test'
        
        if os.path.isfile(testPath) :
            testFile = QuESTTestFile(testPath)
        else:
            logFile.write(fnfWarning.format(testPath))
            print('F',end='')
            continue

        run_test(testFunc, testFile)

    print()


def run_cust_test(testFileName):

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

            run_test(testFunc, testFile)
              
        else:
            logFile.write(fnfWarning.format(testPath))
            print('F',end='')


    elif testPath.endswith('.py'): # Run custom test scripts
        testPath, ext = os.path.splitext(os.path.split(testPath)[-1])
        templib = importlib.import_module(testPath)
        templib.run_tests()
        del templib
    else:
        raise IOError('Unrecognised filetype in test run of file {}'.format(testPath))
            
    print()


def gen_test(testFunc, testFile):
    with open(testFile,'w') as outputFile:
        
        outputFile.write('# {}'.format(testFunc.funcname))
        # Standard run 3 tests
        outputFile.write('3')
        
        for qubitType in "ZPD":
            nQubits = 2
            Qubits = argQureg(nQubits, qubitType)
            
            
    print(testFunc.funcname)
        
def gen_tests():
    from testset import tests
    for testSet in ["math_operations", "stnd_operations","cont_operations", "mcon_operations"]:
        for testFunc in tests[testSet] :
            gen_test(testFunc, unitPath+testFunc.funcname+".test")
    
# If our argument is generate
if len(argv) > 1 and argv[1] == "generate":
    gen_tests()
    quit()
    
Env = createQuESTEnv()

# Load the test sets
from testset import *

if (len(argv) == 1): # Default to "all"
    testsToRun = testSets['all']
else: # Run through command line args
    testsToRun = []
    for test in argv[1:]:
        testsToRun += testSets.get(test,[test])
        
for test in testsToRun:
    if test in testSets or test in tests:
        try:
            run_std_test(tests[test],test)
        except KeyError:
            print("Error: function '{}' does not exist, are you sure you wrote it correctly?".format(test))
    else:
        run_cust_test(test)
        
logFile.close()

