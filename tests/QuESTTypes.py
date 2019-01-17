from ctypes import *

# Declare several warnings which may occur
argWarning  = 'Bad argument list in {0:s} expected {1:d}, recieved {2:d} \n'
fileWarning = '{message} in {file} at line {line} \n'
fnfWarning  = 'File {} not found \n'
funWarning  = 'Function {} does not exist \n'
typeWarning = 'Unrecognised type {} requested in function {} \n'

QuESTLib = CDLL('./QuEST.so')

qreal = c_double

class QASMLogger(Structure):
    _fields_ = [("buffer",c_char_p),
               ("bufferSize",c_int),
               ("bufferFill",c_int),
               ("isLogging",c_int)]

class ComplexArray(Structure):
    _fields_ = [("real", POINTER(qreal)),
               ("imag", POINTER(qreal))]

class Complex(Structure):
    __str__ = lambda self:"({},{})".format(self.real,self.imag)
    _fields_ = [("real",qreal),
                ("imag",qreal)]

class ComplexMatrix2(Structure):
    __str__ = lambda self:"["+"({},{}),"*3+"({},{})]".format(
        self.r0c0.real,self.r0c0.imag,
        self.r1c1.real,self.r1c1.imag,
        self.r1c0.real,self.r1c0.imag,
        self.r1c1.real,self.r1c1.imag)
    _fields_ = [("r0c0",Complex),("r0c1",Complex),
                ("r1c0",Complex),("r1c1",Complex)]

class Vector(Structure):
    __str__ = lambda self:"[{},{},{}]".format(self.x,self.y,self.z)
    _fields_ = [("x",qreal),("y",qreal),("z",qreal)]

class Qureg(Structure):
    _fields_ = [("isDensityMatrix", c_int),
                ("numQubitsRepresented", c_int),
                ("numQubitsInStateVec", c_int),
                ("numAmpsPerChunk",c_longlong),
                ("numAmpsTotal",   c_longlong),
                ("chunkId", c_int),
                ("numChunks", c_int),
                ("stateVec", ComplexArray),
                ("pairStateVec", ComplexArray),
                ("deviceStateVec", ComplexArray),
                ("firstLevelReduction",POINTER(qreal)),("secondLevelReduction",POINTER(qreal)),
                ("qasmLog",POINTER(QASMLogger))]

class QuESTEnv(Structure):
    _fields_ = [("rank",c_int),("numRanks",c_int)]

def stringToList(a):
    a = a.split(',')
    try :
        return list(map(float, a))
    except ValueError:
        raise IOError('Bad array in input file')

def stringToComplex(a):
    a=a.lstrip('(').rstrip(')')
    return list(map(float,a.split(',')))

def argVector(arg):
    Vector(*stringToList(arg))

def argComplexMatrix2(arg):
    vals = stringToList(arg)
    elements = []
    for i in range(0,len(vals),2):
        elements.append(Complex(vals[i],vals[i+1]))
    return ComplexMatrix2(*elements)

def argComplex(arg):
    Complex(*stringToList(arg))

def argComplexArray(arg):
    vals = stringToList(arg)
    real = vals[0::2]
    imag = vals[1::2]
    return ComplexArray(byref(real),byref(imag))

class QuESTTestee:
    basicTypeConv = {"c_int":int, "c_long":int, "c_longlong":int, "qreal":float, "Vector":argVector, "ComplexMatrix2":argComplexMatrix2, "ComplexArray":argComplexArray, "Complex":argComplex }

    funcsList = []
    
    def __init__(self, funcname=None, retType=None, argType=[], defArg=[]):
        self.funcname = funcname
        if not QuESTLib[funcname]:
            raise IOError(funcname+' not found in QuEST API')
        self.thisFunc = QuESTLib[funcname]

        if self.funcname not in list_funcnames():
            QuESTTestee.funcsList.append(self)
        else:
            raise IOError(funcname+' already defined')
        
        self.thisFunc.restype = retType
        self.thisFunc.argtypes = argType
        self.nArgs = len(argType) or 0
        self.defArg = defArg
        
        if self.defArg is not None and len(self.defArg) != self.nArgs:
            raise IOError(argWarning.format(self.funcname, self.nArgs, len(self.defArg)))
        
    def __call__(self,*argsList):
        # If packed as list, otherwise receive as variables
        if len(argsList) == 1 and isinstance(argsList[0],list):
            specArg = argsList[0]
        else:
            specArg = list(argsList)
        self.fix_types(specArg)
        if (len(specArg) == 0 and self.nArgs != 0) or (self.nArgs == 0):
            return self.thisFunc(*self.defArg)
        elif isinstance(specArg,list) and len(specArg) == self.nArgs:
            return self.thisFunc(*specArg)
        else:
            raise IOError(argWarning.format(self.funcname, self.nArgs, len(specArg)))

    def fix_types(self,args):
        for i in range(self.nArgs):
            reqType = self.thisFunc.argtypes[i]
            reqTypeName = self.thisFunc.argtypes[i].__name__
            if isinstance(args[i],reqType):
                pass
            elif reqTypeName in QuESTTestee.basicTypeConv: 
                args[i] = QuESTTestee.basicTypeConv[reqTypeName](args[i])
            else:
                raise IOError(typeWarning.format(reqTypeName, self.funcname))

def list_funcs():
    return QuESTTestee.funcsList

def list_funcnames():
    return list(map(lambda x: x.funcname, QuESTTestee.funcsList))
    
# Define some simple basic constants
complex0 = Complex(0.,0.)
complex1 = Complex(1.,0.)
complexi = Complex(0.,1.)
unitMatrix = ComplexMatrix2(complex1,complex0,complex0,complex1)
xDir = Vector(1.,0.,0.)
yDir = Vector(0.,1.,0.)
zDir = Vector(0.,0.,1.)
