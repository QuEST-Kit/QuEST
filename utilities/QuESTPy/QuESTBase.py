"""
Module to provide various warnings and set up import of QuESTPy
"""

from ctypes import CDLL
import os.path
import importlib

QuESTLib = None
QuESTVer = ""
def init_QuESTLib(QuESTPath="", ext=".so"):
    """
    Import the QuEST library
    """
    global QuESTLib
    if not os.path.isdir(QuESTPath): # If we don't have a direct path to the file
        QuESTPath = QuESTPath.rstrip('/ ') + "/libQuEST" + ext
    if not os.path.isfile(QuESTPath):
        raise FileNotFoundError(fnfWarning.format(QuESTPath))
    QuESTLib = CDLL(QuESTPath)

#Autoinitialise if QuESTLibDir
if importlib.util.find_spec('.QuESTLibDir', package='QuESTPy') is not None:
    _defLib = importlib.import_module('.QuESTLibDir', package='QuESTPy')
    if _defLib.defaultQuESTDir and \
       os.path.isfile(_defLib.defaultQuESTDir.rstrip('/ ') + "/libQuEST.so"):
        init_QuESTLib(_defLib.defaultQuESTDir)

# Declare several warnings which may occur
argWarning = 'Bad argument list in {0:s} expected {1:d}, recieved {2:d} \n'
argWarningGen = 'Bad argument list in {} expected {}, recieved {} \n'
fileWarning = '{message} in {file} at line {line} \n'
fnfWarning = 'File {} not found \n'
funWarning = 'Function {} does not exist \n'
typeWarning = 'Unrecognised type {} requested in function {} \n'
nQubitsWarning = ("Warning: Printing {nStates:,} states may take a long time to run"
                  "and will result in large output."
                  "\n Estimated Size: {sizeEst} {unit}")
compArrKeyWarning = "Key of type {} cannot be selected from ComplexArray"
