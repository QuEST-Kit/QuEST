from QuESTFunc import *

tests = {}
testSets = {}

tests["essential"] = [createQureg, createDensityQureg, destroyQureg, seedQuEST, initZeroState, initPlusState, initClassicalState, initDebugState, setAmps]


tests["init_operations"] = [ initZeroState, initPlusState, initClassicalState, initPureState, initStateFromAmps, initDebugState, setAmps ]

tests["stnd_operations"] = [ collapseToOutcome, compactUnitary, hadamard, pauliX, pauliY, pauliZ, phaseShift, rotateAroundAxis, rotateX, rotateY, rotateZ, sGate, tGate, unitary ]
    
tests["cont_operations"] = [ controlledCompactUnitary, controlledNot, controlledPauliY, controlledPhaseFlip, controlledPhaseShift, controlledRotateAroundAxis, controlledRotateX, controlledRotateY, controlledRotateZ, controlledUnitary ]
        
tests["mcon_operations"] = [ multiControlledPhaseFlip, multiControlledPhaseShift, multiControlledUnitary ]
        
tests["denm_operations"] = [ addDensityMatrix, getDensityAmp, calcPurity, applyOneQubitDephaseError, applyOneQubitDepolariseError, applyTwoQubitDephaseError, applyTwoQubitDepolariseError ]
        
tests["math_operations"] = [ calcFidelity, calcInnerProduct, calcProbOfOutcome, calcTotalProb, getAmp, getNumAmps, getImagAmp, getProbAmp, getRealAmp, getNumQubits, measure, measureWithStats]

tests["pyth_generate"] = [ addDensityMatrix ]
tests["don't_generate"] = tests["essential"] + [ destroyQuESTEnv, createQuESTEnv, calcFidelity, calcInnerProduct, measure, measureWithStats ]

testSets["essential"] = ["essential"]
testSets["qubit_operations"] = ["stnd_operations","cont_operations", "mcon_operations"]
testSets["stnd_operations"] = ["stnd_operations"]
testSets["cont_operations"] = ["cont_operations"]
testSets["mcon_operations"] = ["mcon_operations"]
testSets["denm_operations"] = ["denm_operations"]
testSets["math_operations"] = ["math_operations"]
testSets["stnd_algorithms"] = ["QFT"]
testSets["extended_tests"]  = ["rotate_test"]
testSets["benchmarks"]      = ["rotate_benchmark"]
testSets["all"] = ["stnd_operations", "cont_operations", "mcon_operations", "denm_operations", "math_operations"]
printSets = list(testSets.keys())[:]

# Build individual function tests list
for testFunc in list_funcs():
    name = testFunc.funcname
    testSets[name] = [name]
    tests[name] = [testFunc]
del name


