from ctypes import *
from QuESTTypes import *
from QuESTFunc import *
from sys import argv
import os.path

class QuESTTestFile:
    """ Class containing test file information """
    
    def __init__(self,filename):
        self.File = open(filename,'r')
        
        self.name = filename[filename.rfind('/')+1:] # Remove path
        self.nLine= 0
        try:
            temp = self.readline()
            self.nTests = int(temp)
        except:
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


def run_test(testFuncsList,name=''):

    print('Running tests '+name+":", end='')
    
    for testFunc in testFuncsList:

        logFile.write('\nRunning test {}\n'.format(testFunc.funcname))
        
        testPath = unitPath+testFunc.funcname+'.test'
        
        if os.path.isfile(testPath) :
            testFile = QuESTTestFile(testPath)
        else:
            logFile.write(fnfWarning.format(testPath))
            print('F',end='')
            continue
        
        for test in range(testFile.nTests):
            line = testFile.readline()
        
            qubitType,nBits,*args = line.split()
        
            nBits = int(nBits)
                
            #Upcase qubitType
            qubitType = qubitType.upper()
            # Initialise Qubits
            Qubits = createQureg([nBits, Env])
        
            args.insert(0,Qubits)
        
            nStates = getNumAmps([Qubits])
            
            if qubitType not in qubitTypes:
                raise IOError(file.format(message = 'Unrecognised qubit initialisation state "'+qubitType+'"',
                                          file = testFile.name, line=testFile.nLine))
            
            elif qubitType == "C": # Handle custom initialisation
                qAmps, args = args[0:nBits-1], args[nBits:] # Slice custom amplitudes off the front of args
                qAmpsReal, qAmpsImag = qAmps[0::2], qAmps[1::2] # Split into real and imag
            
                setAmps([Qubits, 0, byref(qAmpsReal), byref(qAmpsImag), nBits])
                
                del qAmps, qAmpsReal, qAmpsImag
            
            else:
                qubitTypes[qubitType]([Qubits])
        
            testFunc(args)
        
            success = True
            for state in range(nStates): # Compare final with expected states
                resultState = Complex(getRealAmp([Qubits,state]),getImagAmp([Qubits,state]))
                expectState = Complex(*list(map(float,testFile.readline().split())))
                logFile.write('{:f} {:f}\n {:f} {:f}\n diff {:f} {:f}\n'.format(resultState.real, resultState.imag,
                                                                              expectState.real, expectState.imag,
                                                                              abs(expectState.real-resultState.real),
                                                                              abs(expectState.imag-resultState.imag))
                              )
                if abs(resultState.real - expectState.real) > tol or abs(resultState.imag - expectState.imag) > tol:
                    success = False
                    break
                
            if success:
                print('.',end='')
            else:
                print('F',end='')
                
            destroyQureg([Qubits])
            
        del testFile
    print()
    
Env = createQuESTEnv()

# Load the test sets
import testset

logFile = open('QuESTTest.log','w')
tol = 1e-6
qubitTypes = {"Z":initZeroState,"P":initPlusState,"D":initDebugState,"C":setAmps}

unitPath = "unitPy/"

testFuncsList = [hadamard]

if (len(argv) == 1):
    testsToRun = ['all']
else:
    testsToRun = argv[1:]

for test in testsToRun:
    run_test(testset.tests[test],test)
    
logFile.close()
