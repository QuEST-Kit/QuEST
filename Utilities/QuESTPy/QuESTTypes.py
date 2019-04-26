from .QuESTBase import *

if QuESTLib is None:
    raise IOError('QuESTLib not initialised')

import random
import math

# Set QuEST default floating point type
QuESTPrecFunc = QuESTLib['getQuEST_PREC']
QuESTPrecFunc.restype = c_int
QuESTPrec = QuESTPrecFunc()

if   QuESTPrec == 1: qreal = c_float
elif QuESTPrec == 2: qreal = c_double
elif QuESTPrec == 4: qreal = c_longdouble
else: raise TypeError('Unable to determine precision of qreal')

class QASMLogger(Structure):
    _fields_ = [("buffer",c_char_p),
               ("bufferSize",c_int),
               ("bufferFill",c_int),
               ("isLogging",c_int)]
    
class ComplexArray(Structure):
    _fields_ = [("real", POINTER(qreal)),
               ("imag", POINTER(qreal))]

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start or 0
            stop  = key.stop
            step  = key.step  or 1
            return [Complex(self.real[i],self.imag[i]) for i in range(start, stop, step)]
        if isinstance(key, int):
            return Complex(self.real[key], self.imag[key])

class Complex(Structure):
    if qreal is c_float:          __repr__ = lambda self:"({:9.7f},{:9.7f})".format(self.real,self.imag)
    elif qreal is c_double:       __repr__ = lambda self:"({:15.13f},{:15.13f})".format(self.real,self.imag)
    elif qreal is c_longdouble:   __repr__ = lambda self:"({:20.18f},{:20.18f})".format(self.real,self.imag)
    __add__ = lambda self, b: Complex(self.real+b.real, self.imag+b.imag)
    __sub__ = lambda self, b: Complex(self.real-b.real, self.imag-b.imag)
    __mul__ = lambda self, b: Complex(self.real*b.real - self.imag*b.imag, self.real*b.imag + self.imag*b.real)
    __truediv__ = lambda self, b: Complex(self.real*b.real + self.imag*b.imag / (b.real*b.real + b.imag*b.imag),
                                      self.imag*b.real - self.real*b.imag / (b.real*b.real + b.imag*b.imag))
    conj = lambda self: Complex(self.real, -self.imag)
    __abs__ = lambda self: math.sqrt( (self*self.conj()).real )
    _fields_ = [("real",qreal),
                ("imag",qreal)]

class ComplexMatrix2(Structure):
    __repr__ = lambda self:"[{},{},{},{})]".format(self.r0c0,self.r0c1,self.r1c0,self.r1c1)
    __abs__ = lambda self: abs(self.r0c0*self.r1c1 - self.r1c0*self.r0c1)
    _fields_ = [("r0c0",Complex),("r0c1",Complex),
                ("r1c0",Complex),("r1c1",Complex)]

class Vector(Structure):
    __str__ = lambda self:"[{},{},{}]".format(self.x, self.y, self.z)
    __add__ = lambda self, b: Vector(self.x+b.x, self.y+b.y, self.z+b.z)
    __sub__ = lambda self, b: Vector(self.x-b.x, self.y-b.y, self.z-b.z)
    _fields_ = [("x",qreal),("y",qreal),("z",qreal)]

class Qureg(Structure):
    def __str__(self):
        if not self._size_warn(): return ""
        retString = ""
        for state in range(self.numAmpsPerChunk):
            retString += Complex(self.stateVec.real[state], self.stateVec.imag[state]).__str__()+"\n"
        return retString.rstrip()

    def _force_string(self):
        retString = ""
        for state in range(self.numAmpsPerChunk):
            retString += Complex(self.stateVec.real[state], self.stateVec.imag[state]).__str__()+"\n"
        return retString.rstrip()

    def _state_vec(self):
        for state in range(self.numAmpsPerChunk):
            yield Complex(self.stateVec.real[state], self.stateVec.imag[state]).__str__()+"\n"
    
    def _size_warn(self, numElem = None, maxElem = 2**5):
        if numElem is None:
            numElem = self.numAmpsTotal
        if numElem > maxElem:
            size = self._size_est(numElem)
            print(nQubitsWarning.format(nStates = numElem, sizeEst = size[0], unit = size[1]), flush=True)
            while True:
                user = input('Do you wish to continue? [y/N]')
                if user is None or user in "Nn":
                    return False
                elif user in "Yy":
                    return True
                else:
                    continue
        return True

    def __setitem__(self, a, b):
        raise TypeError('Manipulating the state vector is not permitted')

    def __set__(self, a, b):
        raise TypeError('Manipulating the state vector is not permitted')

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start or 0
            stop  = key.stop
            step  = key.step  or 1
            if not self._size_warn((stop-start)//step): return None
                    
        return self.stateVec[key]
    
    def _size_est(self, numElem, unit = "MB"):
        charsInStateVecLine = len(str(Complex(0.,0.)))
        sizes = ["YB","ZB","EB","PB","TB","GB","MB","kB","B"]
        i = 0
        size = [ (charsInStateVecLine*numElem), sizes.pop() ]

        while (size[0] > 1024):
            size = [ size[0] >> 10, sizes.pop() ]
        return size

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
    """ Turn a comma-separated string into a list of floats """
    if not isinstance(a, str): raise TypeError(argWarningGen.format('stringToList',str.__name__,type(a).__name__))
    a = a.strip(',').split(',')
    try :
        return list(map(float, a))
    except ValueError:
        raise IOError('Bad array in input file')

def stringToListInt(a):
    """ Turn a comma-separated string into a list of ints """
    
    if not isinstance(a, str): raise TypeError(argWarningGen.format('stringToList',str.__name__,type(a).__name__))
    a = a.split(',')
    try :
        return list(map(int, a))
    except ValueError:
        raise IOError('Bad array in input file')

def stringToComplex(a):
    """ Turn a comma-separated, round-bracketed pair of values into Complex """
    if not isinstance(a, str): raise TypeError(argWarningGen.format('stringToList',str.__name__,type(a).__name__))
    a=a.lstrip('(').rstrip(')')
    return Complex(* list(map(float,a.split(','))))

def argVector(arg):
    if   isinstance(arg, Vector): return arg
    elif isinstance(arg, list): return Vector(*arg)
    elif isinstance(arg, tuple): return Vector(*arg)
    elif isinstance(arg, str):  return Vector(*stringToList(arg))
    else : raise TypeError(argWarningGen.format('argVector','str, tuple or list',type(a).__name__))
    
def argComplexMatrix2(arg):
    if   isinstance(arg, ComplexMatrix2): return arg
    elif isinstance(arg, list): vals = arg
    elif isinstance(arg, tuple): vals = arg
    elif isinstance(arg, str):  vals = stringToList(arg)
    else : raise TypeError(argWarningGen.format('argComplexMatrix2','str, tuple or list',type(a).__name__))
    elements = []
    for i in range(0,len(vals),2):
        elements.append(Complex(vals[i],vals[i+1]))
    if len(elements) != 4: raise TypeError(argWarningGen.format(
            'argComplexMatrix2','4 arguments',len(elements)+" arguments"))
    return ComplexMatrix2(*elements)

def argComplex(arg):
    if   isinstance(arg, Complex): return arg
    elif isinstance(arg, list):  return Complex(*arg)
    elif isinstance(arg, tuple): return Complex(*arg)
    elif isinstance(arg, str):   return stringToComplex(arg)
    else : raise TypeError(argWarningGen.format('argComplex','str, tuple or list',type(a).__name__))

def argComplexArray(arg):
    if   isinstance(arg, ComplexArray): return arg
    elif isinstance(arg, list) or isinstance(arg, tuple):
        if all(map(lambda x: isinstance(x,Complex))):
            vals = stringToList(",".join(arg))
        elif all(map(lambda x: isinstance(x,float))):
            vals = arg
        else: raise TypeError(argWarningGen.format('argComplexArray','list of float/Complex',type(a).__name__))
    elif isinstance(arg, str):   vals = stringToList(arg)
    else : raise TypeError(argWarningGen.format('argComplexArray','str, tuple or list',type(a).__name__))

    real = vals[0::2]
    imag = vals[1::2]
    return ComplexArray(byref(real),byref(imag))

def argPointerQreal(arg):
    if isinstance(arg,str):
        arg = stringToList(arg)
    if isinstance(arg, float):
        arg = qreal(arg)
        
    if isinstance(arg,list) or isinstance(arg, tuple):
        newArg = (qreal*len(arg))()
        for i in range(len(arg)):
            newArg[i] = arg[i]
        return newArg
    elif isinstance(arg, qreal):
        return arg
    else : raise TypeError(argWarningGen.format('argVector','str, float, tuple or list',type(a).__name__))

def argPointerInt(arg):
    if isinstance(arg,str):
        arg = stringToListInt(arg)

    if isinstance(arg,list):
        return (c_int*len(arg))(*arg)
    elif isinstance(arg, c_int):
        return arg

def argPointerLongInt(arg):
    if isinstance(arg,str):
        arg = stringToListInt(arg)

    if isinstance(arg,list):
        return (c_long*len(arg))(*arg)
    elif isinstance(arg, c_int):
        return arg


class QuESTTestee:
    """ Extract function from QuEST C API and set it up as a Python callable object with type checking """
    _basicTypeConv = {"c_int":int, "c_long":int, "c_ulong":int, "c_longlong":int,
                     "c_float":float, "c_double":float, "c_longdouble":float,
                     "Vector":argVector, "ComplexMatrix2":argComplexMatrix2, "ComplexArray":argComplexArray,
                     "Complex":argComplex, "LP_c_double":argPointerQreal, "LP_c_int":argPointerInt,
                     "LP_c_long":argPointerLongInt }

    _funcsList = []
    _funcsDict = {}

    def __init__(self, funcname=None, retType=None, argType=[], defArg=[], denMat=None):
        self.funcname = funcname
        if not QuESTLib[funcname]:
            raise IOError(funcname+' not found in QuEST API')
        self.thisFunc = QuESTLib[funcname]

        if self.funcname not in list_funcnames():
            QuESTTestee._funcsList.append(self)
            QuESTTestee._funcsDict[self.funcname] = self
        else:
            return

        self.target = None
        self.targetType = None
        self.control = None
        self.controlType = None
        self.nControl = None
        
        # Handle shorthands
        for i in range(len(argType)):
            arg = argType[i]

            if isinstance(arg, tuple):
                arg, status = arg

                if status == "targetQubit":
                    self.target= i
                    self.targetType = "Qubit"
                elif status == "targetIndex":
                    self.target= i
                    self.targetType = "Index"
                elif status == "controlQubit":
                    self.control= i
                    self.controlType = "Single"
                elif status == "controlQubits":
                    self.control= i
                    self.controlType = "Multi"
                elif status == "numControlQubits":
                    self.nControl = i
            argType[i] = arg

        self.thisFunc.restype = retType
        self.thisFunc.argtypes = argType
                
        self.nArgs = len(argType) or 0
        self.defArg = defArg
        self.denMat = denMat

        if self.defArg is not None and len(self.defArg) != self.nArgs:
            raise IOError(argWarning.format(self.funcname, self.nArgs, len(self.defArg)))

    def __call__(self,*argsList):
        # If packed as list, otherwise receive as variables
        if len(argsList) == 1 and isinstance(argsList[0],list):
            specArg = argsList[0]
        else:
            specArg = list(argsList)

        if (len(specArg) == 0 and self.nArgs != 0) or (self.nArgs == 0):
            self._fix_types(self.defArg)
            return self.thisFunc(*self.defArg)
        elif isinstance(specArg,list) and len(specArg) == self.nArgs:
            self._fix_types(specArg)
            return self.thisFunc(*specArg)
        else:
            raise IOError(argWarning.format(self.funcname, self.nArgs, len(specArg)))

    def _fix_types(self,args):
        for i in range(self.nArgs):
            reqType = self.thisFunc.argtypes[i]
            reqTypeName = self.thisFunc.argtypes[i].__name__
            if isinstance(args[i],reqType):
                pass
            elif reqTypeName in QuESTTestee._basicTypeConv:
                args[i] = QuESTTestee._basicTypeConv[reqTypeName](args[i])
            else:
                print(args[i], reqTypeName)
                raise IOError(typeWarning.format(reqTypeName, self.funcname))

    def get_func(name):
            return QuESTTestee._funcsDict[name]
        
def dict_funcs():
    return QuESTTestee._funcsDict

def list_funcs():
    return QuESTTestee._funcsList

def list_funcnames():
   return list(map(lambda x: x.funcname, QuESTTestee._funcsList))

# Define some simple basic constants
complex0 = Complex(0.,0.)
complex1 = Complex(1.,0.)
complexi = Complex(0.,1.)
complexHalf = Complex(0.5,0.)
complexSqr2 = Complex(math.sqrt(2),0.0)
complexRSqr2 = Complex(1./math.sqrt(2),0.0)
idenMatrix = ComplexMatrix2(complex1,complex0,complex0,complex1)
xDir = Vector(1.,0.,0.)
yDir = Vector(0.,1.,0.)
zDir = Vector(0.,0.,1.)

def rand_norm_comp():
    newComplex = Complex(random.random(), random.random())
    norm = abs(newComplex)
    newComplex.imag /= norm
    newComplex.real /= norm
    return newComplex

def rand_norm_comp_pair():
    return rand_norm_comp()*complexRSqr2, rand_norm_comp()*complexRSqr2

def rand_unit_mat():
    elems = list(rand_norm_comp_pair())
    elems += [Complex(0,0)-elems[1].conj()]
    elems += [elems[0].conj()]
    newMat = ComplexMatrix2(*elems)

    return newMat
