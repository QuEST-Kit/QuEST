from ctypes import *
import os.path

QuESTLib = None
QuESTVer = ""
def init_QuESTLib(QuESTPath = ""):
    global QuESTLib
    QuESTPath = QuESTPath.rstrip('/ ') + "/libQuEST.so"
    if not os.path.isfile(QuESTPath):
        raise OSError(fnfWarning.format(QuESTPath))
    QuESTLib = CDLL(QuESTPath)

# Declare several warnings which may occur
argWarning  = 'Bad argument list in {0:s} expected {1:d}, recieved {2:d} \n'
argWarningGen  = 'Bad argument list in {} expected {}, recieved {} \n'
fileWarning = '{message} in {file} at line {line} \n'
fnfWarning  = 'File {} not found \n'
funWarning  = 'Function {} does not exist \n'
typeWarning = 'Unrecognised type {} requested in function {} \n'
nQubitsWarning = 'Warning: Printing {nStates:,} states may take a long time to run and will result in large output.\n Estimated Size: {sizeEst} MB'
