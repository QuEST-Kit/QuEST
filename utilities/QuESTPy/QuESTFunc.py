from .QuESTTypes import *
import copy

# Aliases for better legibility
_targetQubit = (c_int, "targetQubit")
_controlQubit = (c_int, "controlQubit")
_controlQubits = (POINTER(c_int), "controlQubits")
_numControlQubits = (c_int, "numControlQubits")
_stateIndex = (c_longlong, "targetIndex")


# Data Operations
createQuESTEnv     = QuESTTestee ('createQuESTEnv',QuESTEnv)
destroyQuESTEnv    = QuESTTestee ('destroyQuESTEnv',None,[QuESTEnv],[None])
createQureg        = QuESTTestee ('createQureg',Qureg,[c_int,QuESTEnv],[1,None])
createDensityQureg = QuESTTestee ('createDensityQureg',Qureg,[c_int,QuESTEnv],[1,None])
destroyQureg       = QuESTTestee ('destroyQureg',None,[Qureg,QuESTEnv],[None,None])

# Utility Operations
cloneQureg             = QuESTTestee ("cloneQureg", retType=None, argType=[Qureg, Qureg], defArg=[None, None])
# getEnvironmentString = QuESTTestee ("getEnvironmentString", retType=None, argType=[QuESTEnv, Qureg, c_char*200], defArg=[None])
seedQuEST              = QuESTTestee ("seedQuEST", retType=None, argType=[POINTER(c_long),c_int], defArg=[None,None])
seedQuESTDefault       = QuESTTestee ("seedQuESTDefault", retType=None)

# Reporting Operations
reportQuESTEnv      = QuESTTestee ("reportQuESTEnv", retType=None, argType=[QuESTEnv], defArg=[None])
reportQuregParams   = QuESTTestee ("reportQuregParams", retType=None, argType=[Qureg], defArg=[None])
reportState         = QuESTTestee ("reportState", retType=None, argType=[Qureg], defArg=[None])
reportStateToScreen = QuESTTestee ("reportStateToScreen", retType=None, argType=[Qureg,QuESTEnv,c_int], defArg=[None,None,0]) 

# QASM Operations
clearRecordedQASM       = QuESTTestee ("clearRecordedQASM", retType=None, argType=[Qureg], defArg=[None])
printRecordedQASM       = QuESTTestee ("printRecordedQASM", retType=None, argType=[Qureg], defArg=[None])
startRecordingQASM      = QuESTTestee ("startRecordingQASM", retType=None, argType=[Qureg], defArg=[None])
stopRecordingQASM       = QuESTTestee ("stopRecordingQASM", retType=None, argType=[Qureg], defArg=[None])
writeRecordedQASMToFile = QuESTTestee ("writeRecordedQASMToFile", retType=None, argType=[Qureg,c_char_p], defArg=[None,None]) 

# Parallel Operations
syncQuESTEnv     = QuESTTestee ("syncQuESTEnv", retType=None, argType=[QuESTEnv], defArg=[None])
syncQuESTSuccess = QuESTTestee ("syncQuESTSuccess", retType=c_int, argType=[c_int], defArg=[None]) 

# Initialisation Operations
initZeroState      = QuESTTestee ("initZeroState",      retType=None, argType=[Qureg], defArg=[None])
initPlusState      = QuESTTestee ("initPlusState",      retType=None, argType=[Qureg], defArg=[None])
initClassicalState = QuESTTestee ("initClassicalState", retType=None, argType=[Qureg,c_longlong], defArg=[None,None])
initPureState      = QuESTTestee ("initPureState",      retType=None, argType=[Qureg,Qureg], defArg=[None,None])
initStateFromAmps  = QuESTTestee ("initStateFromAmps",  retType=None, argType=[Qureg,POINTER(qreal),POINTER(qreal)], defArg=[None,None,None], denMat = False)
initStateDebug     = QuESTTestee ("initStateDebug",     retType=None, argType=[Qureg], defArg=[None])
initDebugState     = initStateDebug  # Alias
setAmps            = QuESTTestee ("setAmps",            retType=None, argType=[Qureg,_stateIndex,POINTER(qreal),POINTER(qreal),c_longlong], defArg=[None,None,None,None,None], denMat = False) 
setDensityAmps     = QuESTTestee ("setDensityAmps",     retType=None, argType=[Qureg,POINTER(qreal),POINTER(qreal)], defArg=[None,None,None], denMat = True) 

# Basic Operations
hadamard         = QuESTTestee ("hadamard",         retType=None, argType=[Qureg,_targetQubit], defArg=[None,0])
pauliX           = QuESTTestee ("pauliX",           retType=None, argType=[Qureg,_targetQubit], defArg=[None,0])
pauliY           = QuESTTestee ("pauliY",           retType=None, argType=[Qureg,_targetQubit], defArg=[None,0])
pauliZ           = QuESTTestee ("pauliZ",           retType=None, argType=[Qureg,_targetQubit], defArg=[None,0])
sGate            = QuESTTestee ("sGate",            retType=None, argType=[Qureg,_targetQubit], defArg=[None,0])
tGate            = QuESTTestee ("tGate",            retType=None, argType=[Qureg,_targetQubit], defArg=[None,0]) 
compactUnitary   = QuESTTestee ("compactUnitary",   retType=None, argType=[Qureg,_targetQubit,Complex,Complex], defArg = [None, 0, *rand_norm_comp_pair()])
phaseShift       = QuESTTestee ("phaseShift",       retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,random.uniform(0.,360.)])
rotateAroundAxis = QuESTTestee ("rotateAroundAxis", retType=None, argType=[Qureg,_targetQubit,qreal,Vector], defArg=[None,0,random.uniform(0.,360.),xDir])
rotateX          = QuESTTestee ("rotateX",          retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,random.uniform(0.,360.)])
rotateY          = QuESTTestee ("rotateY",          retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,random.uniform(0.,360.)])
rotateZ          = QuESTTestee ("rotateZ",          retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,random.uniform(0.,360.)])
unitary          = QuESTTestee ("unitary",          retType=None, argType=[Qureg,_targetQubit,ComplexMatrix2], defArg=[None,0,rand_unit_mat()]) 

# Controlled Operations
controlledCompactUnitary   = QuESTTestee ("controlledCompactUnitary",   retType=None, argType=[Qureg,_controlQubit,_targetQubit,Complex,Complex], defArg=[None,1,0,*rand_norm_comp_pair()])
controlledNot              = QuESTTestee ("controlledNot",              retType=None, argType=[Qureg,_controlQubit,_targetQubit], defArg=[None,1,0])
controlledPauliY           = QuESTTestee ("controlledPauliY",           retType=None, argType=[Qureg,_controlQubit,_targetQubit], defArg=[None,1,0])
controlledPhaseFlip        = QuESTTestee ("controlledPhaseFlip",        retType=None, argType=[Qureg,_controlQubit,_targetQubit], defArg=[None,1,0])
controlledPhaseShift       = QuESTTestee ("controlledPhaseShift",       retType=None, argType=[Qureg,_controlQubit,_targetQubit,qreal], defArg=[None,1,0,90.])
controlledRotateAroundAxis = QuESTTestee ("controlledRotateAroundAxis", retType=None, argType=[Qureg,_controlQubit,_targetQubit,qreal,Vector], defArg=[None,1,0,90.,xDir])
controlledRotateX          = QuESTTestee ("controlledRotateX",          retType=None, argType=[Qureg,_controlQubit,_targetQubit,qreal], defArg=[None,1,0,90.])
controlledRotateY          = QuESTTestee ("controlledRotateY",          retType=None, argType=[Qureg,_controlQubit,_targetQubit,qreal], defArg=[None,1,0,90.])
controlledRotateZ          = QuESTTestee ("controlledRotateZ",          retType=None, argType=[Qureg,_controlQubit,_targetQubit,qreal], defArg=[None,1,0,90.])
controlledUnitary          = QuESTTestee ("controlledUnitary",          retType=None, argType=[Qureg,_controlQubit,_targetQubit,ComplexMatrix2], defArg=[None,1,0,rand_unit_mat()]) 

# Multi-controlled Operations
multiControlledPhaseFlip  = QuESTTestee ("multiControlledPhaseFlip",  retType=None, argType=[Qureg,_controlQubits,_numControlQubits], defArg=[None,[0],1])
multiControlledPhaseShift = QuESTTestee ("multiControlledPhaseShift", retType=None, argType=[Qureg,_controlQubits,_numControlQubits,qreal], defArg=[None,[0],1,random.uniform(0.,360.)])
multiControlledUnitary    = QuESTTestee ("multiControlledUnitary",    retType=None, argType=[Qureg,_controlQubits,_numControlQubits,_targetQubit,ComplexMatrix2], defArg=[None,[1],1,0,rand_unit_mat()]) 

# Density Matrix Operations
addDensityMatrix             = QuESTTestee ("addDensityMatrix",             retType=None, argType=[Qureg,qreal,Qureg], defArg=[None,50.,None], denMat=True)
applyOneQubitDephaseError    = QuESTTestee ("applyOneQubitDephaseError",    retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,0.25], denMat=True)
applyOneQubitDepolariseError = QuESTTestee ("applyOneQubitDepolariseError", retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,0.25], denMat=True)
applyOneQubitDampingError    = QuESTTestee ("applyOneQubitDampingError",    retType=None, argType=[Qureg,_targetQubit,qreal], defArg=[None,0,0.25], denMat=True)
applyTwoQubitDephaseError    = QuESTTestee ("applyTwoQubitDephaseError",    retType=None, argType=[Qureg,_targetQubit,_controlQubit,qreal], defArg=[None,0,1,0.25], denMat=True)
applyTwoQubitDepolariseError = QuESTTestee ("applyTwoQubitDepolariseError", retType=None, argType=[Qureg,_targetQubit,_controlQubit,qreal], defArg=[None,0,1,0.25], denMat=True) 

# Examination and Mathematical Operations
calcFidelity      = QuESTTestee ("calcFidelity",      retType=qreal, argType=[Qureg,Qureg], defArg=[None,None])
calcInnerProduct  = QuESTTestee ("calcInnerProduct",  retType=Complex, argType=[Qureg,Qureg], defArg=[None,None])
calcProbOfOutcome = QuESTTestee ("calcProbOfOutcome", retType=qreal, argType=[Qureg,_targetQubit,c_int], defArg=[None,0,1])
calcTotalProb     = QuESTTestee ("calcTotalProb",     retType=qreal, argType=[Qureg], defArg=[None])
calcPurity        = QuESTTestee ("calcPurity",        retType=qreal, argType=[Qureg], defArg=[None],denMat=True)
collapseToOutcome = QuESTTestee ("collapseToOutcome", retType=None, argType=[Qureg,_targetQubit,c_int], defArg=[None,0,0]) 
getAmp            = QuESTTestee ("getAmp",            retType=Complex, argType=[Qureg,_stateIndex], defArg=[None,0], denMat = False)
getDensityAmp     = QuESTTestee ("getDensityAmp",     retType=Complex, argType=[Qureg,_stateIndex,_stateIndex], defArg=[None,0,0],denMat=True)
getNumAmps        = QuESTTestee ("getNumAmps",        retType=c_int, argType=[Qureg], defArg=[None], denMat=False)
getImagAmp        = QuESTTestee ("getImagAmp",        retType=qreal, argType=[Qureg,_stateIndex], defArg=[None,0], denMat = False)
getProbAmp        = QuESTTestee ("getProbAmp",        retType=qreal, argType=[Qureg,_stateIndex], defArg=[None,0], denMat = False)
getRealAmp        = QuESTTestee ("getRealAmp",        retType=qreal, argType=[Qureg,_stateIndex], defArg=[None,0], denMat = False)
getNumQubits      = QuESTTestee ("getNumQubits",      retType=c_int, argType=[Qureg], defArg=[None])
measure           = QuESTTestee ("measure",           retType=c_int, argType=[Qureg,_targetQubit], defArg=[None,0])
measureWithStats  = QuESTTestee ("measureWithStats",  retType=c_int, argType=[Qureg,_targetQubit,POINTER(qreal)], defArg=[None,0,None])
