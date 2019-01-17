from ctypes import *
from QuESTTypes import *
from QuESTFunc import *
from sys import argv
import os.path
import importlib

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
        
    def readline(self):
        """ Reads a line from a test file """
        for line in self.File:
            self.nLine += 1
            line = line[:line.find('#')].strip()
            if line:
                return line
        raise IOError(fileWarning.format(message='Unexpected end of file',file=self.name, line=self.nLine))

    def parse_args(self):
        """ Split arguments, but maintain arrays and complex arrays as a block """
        remBrac = ''.maketrans('[{()}]','      ','[{()}]')
        line = self.readline().translate(remBrac)
        return line.split()

        
def run_test(testFunc, testFile):

    for test in range(testFile.nTests):
        # line = testFile.readline()
    
        qubitType,nBits,*args = testFile.parse_args()
    
        nBits = int(nBits)
            
        #Upcase qubitType
        qubitType = qubitType.upper()
        # Initialise Qubits
        Qubits = createQureg(nBits, Env)
    
        args.insert(0,Qubits)
    
        nStates = getNumAmps(Qubits)
        
        if qubitType not in qubitTypes:
            raise IOError(file.format(message = 'Unrecognised qubit initialisation state "'+qubitType+'"',
                                      file = testFile.name, line=testFile.nLine))
        
        elif qubitType == "C": # Handle custom initialisation
            qAmps, args = args[0:nBits-1], args[nBits:] # Slice custom amplitudes off the front of args
            qAmpsReal, qAmpsImag = qAmps[0::2], qAmps[1::2] # Split into real and imag
        
            setAmps(Qubits, 0, byref(qAmpsReal), byref(qAmpsImag), nBits)
            
            del qAmps, qAmpsReal, qAmpsImag
        
        else:
            qubitTypes[qubitType](Qubits)
    
        testFunc(args)
    
        success = True
        for state in range(nStates): # Compare final with expected states
            resultState = getAmp(Qubits,state)
            expectState = Complex(*list(map(float,testFile.readline().split(','))))
            if abs(resultState.real - expectState.real) > tol or abs(resultState.imag - expectState.imag) > tol:
                logFile.write('{:f} {:f}\n{:f} {:f}\ndiff {:f} {:f}\n'.format(resultState.real, resultState.imag,
                                                                              expectState.real, expectState.imag,
                                                                              abs(expectState.real-resultState.real),
                                                                              abs(expectState.imag-resultState.imag))
                             )
                      
                success = False
                break
            
        if success:
            print('.',end='')
        else:
            print('F',end='')
            
        destroyQureg(Qubits)
        
    del testFile


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
    
Env = createQuESTEnv()

# Load the test sets
from testset import *

logFile = open('QuESTTest.log','w')
tol = 1e-6
qubitTypes = {"Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps}

unitPath = "unitPy/"

if (len(argv) == 1):
    testsToRun = testSets['all']
else:
    testsToRun = []
    for test in argv[1:]:
        testsToRun += testSets.get(test,[test])
        
for test in testsToRun:
    if test in testSets:
        try:
            run_std_test(tests[test],test)
        except KeyError:
            print("Error: function '{}' does not exist, are you sure you wrote it correctly?".format(test))
    else:
        run_cust_test(test)
        
logFile.close()

