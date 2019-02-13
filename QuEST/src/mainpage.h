/** @file containing the mainpage doc */

/** 
\mainpage The QuEST API

This page is currently a mere grouping of the full API doc available at QuEST.h

\section sec_types Datatypes

- \ref QuESTEnv
- \ref Qureg
- \ref qreal
- \ref qcomp
- \ref Complex
- \ref ComplexMatrix2
<br><br>
- \ref toComplex
- \ref fromComplex

\section sec_alloc Allocation

- \ref createDensityQureg
- \ref createQuESTEnv
- \ref createQureg
- \ref destroyQuESTEnv
- \ref destroyQureg
- \ref seedQuEST
- \ref seedQuESTDefault

\section sec_init Initialisation

- \ref initClassicalState
- \ref initPlusState
- \ref initZeroState
- \ref initPureState
- \ref initStateFromAmps

\section sec_gates Gates

- \ref hadamard
- \ref pauliX
- \ref pauliY
- \ref pauliZ
- \ref unitary
- \ref compactUnitary
- \ref phaseShift
- \ref rotateX
- \ref rotateY
- \ref rotateZ
- \ref rotateAroundAxis
- \ref sGate
- \ref tGate
<br><br>
- \ref controlledCompactUnitary
- \ref controlledNot
- \ref controlledPauliY
- \ref controlledPhaseFlip
- \ref controlledPhaseShift
- \ref controlledRotateAroundAxis
- \ref controlledRotateX
- \ref controlledRotateY
- \ref controlledRotateZ
- \ref controlledUnitary
<br><br>
- \ref multiControlledPhaseFlip
- \ref multiControlledPhaseShift
- \ref multiControlledUnitary
<br><br>
- \ref measure
- \ref measureWithStats

\section sec_noise Noise

- \ref applyOneQubitDephaseError
- \ref applyTwoQubitDephaseError
- \ref applyOneQubitDepolariseError
- \ref applyTwoQubitDepolariseError

\section sec_othermods Modifications

- \ref addDensityMatrix
- \ref cloneQureg
- \ref collapseToOutcome
- \ref setAmps

\section sec_properties Properties

- \ref getAmp
- \ref getDensityAmp
- \ref getNumAmps
- \ref getNumQubits
- \ref getProbAmp
- \ref getImagAmp
- \ref getRealAmp

\section sec_calculations Calculations

- \ref calcPurity
- \ref calcTotalProb
- \ref calcFidelity
- \ref calcInnerProduct
- \ref calcProbOfOutcome

\section sec_qasm QASM

- \ref startRecordingQASM
- \ref stopRecordingQASM
- \ref clearRecordedQASM
- \ref printRecordedQASM
- \ref writeRecordedQASMToFile

\section sec_debug Debugging

- \ref reportQuESTEnv
- \ref reportQuregParams
- \ref reportState
- \ref reportStateToScreen

\section sec_other Other

- \ref syncQuESTEnv
- \ref syncQuESTSuccess
- \ref getEnvironmentString

*/