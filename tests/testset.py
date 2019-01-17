from QuESTFunc import *

tests = {}
tests["init_operations"] = [ initZeroState, initPlusState, initClassicalState, initPureState, initStateFromAmps, initDebugState, setAmps ],

tests["stnd_operations"] = [ compactUnitary, hadamard, pauliX, pauliY, pauliZ, phaseShift, rotateAroundAxis, rotateX, rotateY, rotateZ, sGate, tGate, unitary ]
    
tests["cont_operations"] = [ controlledCompactUnitary, controlledNot, controlledPauliY, controlledPhaseFlip, controlledPhaseShift, controlledRotateAroundAxis, controlledRotateX, controlledRotateY, controlledRotateZ, controlledUnitary ]
        
tests["mcon_operations"] = [ multiControlledPhaseFlip, multiControlledPhaseShift, multiControlledUnitary ]
        
tests["denm_operations"] = [ addDensityMatrix, applyOneQubitDephaseError, applyOneQubitDepolariseError, applyTwoQubitDephaseError, applyTwoQubitDepolariseError ]
        
tests["math_operations"] = [ calcFidelity, calcInnerProduct, calcProbOfOutcome, calcTotalProb, calcPurity, collapseToOutcome, getAmp, getDensityAmp, getNumAmps, getImagAmp, getProbAmp, getRealAmp, getNumQubits, measure, measureWithStats ]
        
tests["all"] = tests["stnd_operations"] + tests["cont_operations"] + tests["mcon_operations"] + tests["denm_operations"] + tests["math_operations"]

tests["qubit_operations"] = tests["stnd_operations"] + tests["cont_operations"] + tests["mcon_operations"]
    

