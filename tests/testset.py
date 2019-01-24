from QuESTFunc import *

tests = {}
testSets = {}

tests["essential"] = [createQuESTEnv, destroyQuESTEnv, createQureg, createDensityQureg, destroyQureg, seedQuEST, initZeroState, initPlusState, initClassicalState, initDebugState, setAmps]


tests["init_operations"] = [ initZeroState, initPlusState, initClassicalState, initPureState, initStateFromAmps, initDebugState, setAmps ]

tests["stnd_operations"] = [ collapseToOutcome, compactUnitary, hadamard, pauliX, pauliY, pauliZ, phaseShift, rotateAroundAxis, rotateX, rotateY, rotateZ, sGate, tGate, unitary ]
    
tests["cont_operations"] = [ controlledCompactUnitary, controlledNot, controlledPauliY, controlledPhaseFlip, controlledPhaseShift, controlledRotateAroundAxis, controlledRotateX, controlledRotateY, controlledRotateZ, controlledUnitary ]
        
tests["mcon_operations"] = [ multiControlledPhaseFlip, multiControlledPhaseShift, multiControlledUnitary ]
        
tests["denm_operations"] = [ addDensityMatrix, getDensityAmp, calcPurity, applyOneQubitDephaseError, applyOneQubitDepolariseError, applyTwoQubitDephaseError, applyTwoQubitDepolariseError ]
        
tests["math_operations"] = [ calcFidelity, calcInnerProduct, calcProbOfOutcome, calcTotalProb, getAmp, getNumAmps, getImagAmp, getProbAmp, getRealAmp, getNumQubits, measure, measureWithStats]

tests["don't_generate"] = tests["essential"] + [ calcFidelity, calcInnerProduct, measure, measureWithStats, addDensityMatrix ]

testSets["essential"] = ["essential"]
testSets["qubit_operations"] = ["stnd_operations","cont_operations", "mcon_operations"]
testSets["all"] = ["stnd_operations", "cont_operations", "mcon_operations", "denm_operations", "math_operations"]

# Build individual function tests list
for testFunc in list_funcs():
    name = testFunc.funcname
    testSets[name] = [name]
    tests[name] = [testFunc]
del name


printSets = testSets["essential"] + testSets["all"]
