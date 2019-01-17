from QuESTFunc import *

tests = {}
testSets = {}
# Underscored names are internal lists of functions, plain names are the calling name from CLI
tests["init_operations"] = [ initZeroState, initPlusState, initClassicalState, initPureState, initStateFromAmps, initDebugState, setAmps ]

tests["stnd_operations"] = [ compactUnitary, hadamard, pauliX, pauliY, pauliZ, phaseShift, rotateAroundAxis, rotateX, rotateY, rotateZ, sGate, tGate, unitary ]
    
tests["cont_operations"] = [ controlledCompactUnitary, controlledNot, controlledPauliY, controlledPhaseFlip, controlledPhaseShift, controlledRotateAroundAxis, controlledRotateX, controlledRotateY, controlledRotateZ, controlledUnitary ]
        
tests["mcon_operations"] = [ multiControlledPhaseFlip, multiControlledPhaseShift, multiControlledUnitary ]
        
tests["denm_operations"] = [ addDensityMatrix, applyOneQubitDephaseError, applyOneQubitDepolariseError, applyTwoQubitDephaseError, applyTwoQubitDepolariseError ]
        
tests["math_operations"] = [ calcFidelity, calcInnerProduct, calcProbOfOutcome, calcTotalProb, calcPurity, collapseToOutcome, getAmp, getDensityAmp, getNumAmps, getImagAmp, getProbAmp, getRealAmp, getNumQubits, measure, measureWithStats ]

# Build individual function tests list
for testFunc in list_funcs():
    name = testFunc.funcname
    testSets[name] = [name]
    tests[name] = [testFunc]
del name


testSets["all"] = ["stnd_operations", "cont_operations", "mcon_operations", "denm_operations", "math_operations"]
testSets["qubit_operations"] = ["stnd_operations","cont_operations", "mcon_operations"]
