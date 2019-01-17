from QuESTFunc import *
import math

def run_tests():
    testRot()

def compareStates(a, b, tol):
    if (a.numQubitsRepresented != b.numQubitsRepresented):
        raise IndexError('A and B registers are not the same size')
    
    for state in range(getNumAmps(a)): # Compare final with expected states
        aState = getAmp(a,state)
        bState = getAmp(b,state)
        if abs(aState.real - bState.real) > tol or abs(aState.imag - bState.imag) > tol:
            return False
    return True

def compareReals(a, b, tol):
    if abs(a - b) > tol: return False
    return True

def testRot():
    COMPARE_PRECISION = 1e-6
    env = createQuESTEnv()

    passed=1

    numQubits=10
    rotQubit = 0

    angs = [1.2,-2.4,0.3]

    alpha = Complex(math.cos(angs[0]) * math.cos(angs[1]),
                    math.cos(angs[0]) * math.sin(angs[1]) )
    beta  = Complex(math.sin(angs[0]) * math.cos(angs[2]),
                    math.sin(angs[0]) * math.sin(angs[2]) )

    mq = createQureg(numQubits, env)
    mqVerif = createQureg(numQubits, env)

    initStateDebug(mq)
    initStateDebug(mqVerif)
    for i in range(numQubits):
        rotQubit=i
        compactUnitary(mq, rotQubit, alpha, beta)
    
    # note -- this is only checking if the state changed at all due to rotation,
    # not that it changed correctly
    if passed: passed = not compareStates(mq, mqVerif, COMPARE_PRECISION)

    # Rotate back the other way and check we arrive back at the initial state
    # (conjugate transpose of the unitary)
    alpha.imag *= -1
    beta.real  *= -1
    beta.imag  *= -1

    # (order of qubits operated upon doesn't matter)
    for i in range(numQubits):
        rotQubit=i
        compactUnitary(mq, rotQubit, alpha, beta)

    # unitaries are relatively imprecise (10* needed for signle precision)
    if passed: passed = compareStates(mq, mqVerif, 10*COMPARE_PRECISION)

    destroyQureg(mq, env)
    destroyQureg(mqVerif, env)

    # check for normalisation
    numQubits=25
    mq = createQureg(numQubits, env)
    
    initPlusState(mq)
    for i in range(numQubits):
        rotQubit=i
        compactUnitary(mq, rotQubit, alpha, beta)
    
    outcome = calcTotalProb(mq);    
    if passed: passed = compareReals(1.0, outcome, COMPARE_PRECISION)
    destroyQureg(mq, env)

    if passed:
        print('.',end='')
    else:
        print('F',end='')


    
