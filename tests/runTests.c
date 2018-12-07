# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_debug.h"

# define NUM_TESTS 38
# define PATH_TO_TESTS "unit/"
# define VERBOSE 0

// quad precision unit testing is no more stringent than double
# if QuEST_PREC==1
# define COMPARE_PRECISION 10e-5
# else
# define COMPARE_PRECISION 10e-13
# endif


QuESTEnv env;

void reportTest(Qureg qureg, char testName[200]){
    printf("\nTest: %s\n", testName);
    reportStateToScreen(qureg, env, 0);
}

/** returns 1 if equal to within precision */
int compareReals(qreal a, qreal b, qreal precision){
    qreal diff = a-b;
    if (diff<0) diff *= -1;
    if (diff>precision) return 0;
    else return 1;
}

int test_initZeroState(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    initZeroState(mq);

    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);

    passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_initPlusState(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    
    /*
     * state vectors
     */
    Qureg mq, mqVerif; 
    
    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    initPlusState(mq);

    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);
    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    
    // |+> state has every qubit 50% prob in 0
    initPlusState(mq);
    for (int q=0; q < numQubits; q++) {
        qreal prob = calcProbOfOutcome(mq, q, 0);
        if (passed) passed = compareReals(prob, 0.5, COMPARE_PRECISION);
    }
        
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);
    
   /*
    * density matrices
    */
    Qureg dens;
    dens = createDensityQureg(numQubits, env);
    
    initPlusState(dens);
    
    for (int q=0; q < numQubits; q++) {
        qreal prob = calcProbOfOutcome(dens, q, 0);            
        if (passed) passed = compareReals(prob, 0.5, COMPARE_PRECISION);
    }
    
    destroyQureg(dens, env);

    return passed;
}

int test_initClassicalState(char testName[200]){
    int passed=1;
    int numQubits=3;
    int numAmps=1 << numQubits;
    Qureg mq;
    
    /*
     * state-vectors
     */
    mq = createQureg(numQubits, env);

    // test every classical state
    for (long long int stateInd=0LL; stateInd < numAmps; stateInd++) {
        initClassicalState(mq, stateInd);
        
        // check that every other state has prob 0
        for (long long int i=0LL; i < numAmps; i++) {
            if (i == stateInd)
                passed = passed && (getProbAmp(mq,i) == 1.0);
            else
                passed = passed && (getProbAmp(mq,i) == 0.0);
        }
    }

    destroyQureg(mq, env);
    
    /*
     * density matrices
     */
    
    mq = createDensityQureg(numQubits, env);
    
    // test every classical state
    for (long long int stateInd=0LL; stateInd < numAmps; stateInd++) {
        initClassicalState(mq, stateInd);
        
        // check that every qubit has correct probabilities
        for (long long int q=0; q < numQubits; q++) {
            int bit = (stateInd & ( 1LL << q )) >> q;
            qreal probOf1 = calcProbOfOutcome(mq, q, 1);
            if (passed) passed = compareReals(probOf1, bit, COMPARE_PRECISION);
        }
    }
    
    destroyQureg(mq, env);
    
    return passed;
}

int test_initPureState(char testName[200]) {    
    int passed=1;
    int numQubits=3;
    
    /*
     * for statevectors / purestates
     */
    Qureg mq, mqVerif;
    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);
    
    initZeroState(mq);
    initStateDebug(mqVerif);
    
    initPureState(mq, mqVerif);
    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);
        
    
    /*
     * for density matrices / mixed states
     */
    Qureg dens, vec;
    vec = createQureg(numQubits, env);
    dens = createDensityQureg(numQubits, env);

    initPlusState(vec);
    initZeroState(dens);
    
    // @TODO this is causing error on 1-process MPI, due to call to copyVecIntoMatrixPairState
    // set dens = |+><+|
    initPureState(dens, vec);
        
    qreal prob;
    for (int q=0; q < numQubits; q++) {
        prob = calcProbOfOutcome(dens, q, 0);
        if (passed) passed = compareReals(prob, 0.5, COMPARE_PRECISION);
    }

    destroyQureg(vec, env);
    destroyQureg(dens, env);
    
    return passed;
}

int test_setAmps(char testName[200]) {
    
    int passed=1;
    int numQubits=3;
    
    Qureg qureg;
    qureg = createQureg(numQubits, env);
    
    // test writing total state vec
    qreal reals[8] = {1,2,3,4,5,6,7,8};
    qreal imags[8] = {8,7,6,5,4,3,2,1};
    setAmps(qureg, 0, reals, imags, qureg.numAmpsTotal);
    for (long long int i=0; i < 8; i++) {
        if (passed) passed = compareReals(getRealAmp(qureg,i), reals[i], 0);
        if (passed) passed = compareReals(getImagAmp(qureg,i), imags[i], 0);
    }
    
    // test writing only some of statevec
    initZeroState(qureg);
    setAmps(qureg, 2, reals+2, imags+2, 4); // write {3,4,5,6} to inds {2,3,4,5}
    for (long long int i=0; i < 8; i++) {
        
        // indices outside {2,3,4,5} are unchanged from |0> = {1,0,0,0}...
        if (i==0) {
            if (passed) passed = compareReals(getRealAmp(qureg,i), 1, 0);
            if (passed) passed = compareReals(getImagAmp(qureg,i), 0, 0);
        }
        else if (i<2 || i>=6) {
            if (passed) passed = compareReals(getRealAmp(qureg,i), 0, 0);
            if (passed) passed = compareReals(getImagAmp(qureg,i), 0, 0);
        }
        // otherwise, indices should be set to those in reals and imag
        else {
            if (passed) passed = compareReals(getRealAmp(qureg,i), reals[i], 0);
            if (passed) passed = compareReals(getImagAmp(qureg,i), imags[i], 0);
        }
    }
    
    destroyQureg(qureg, env);
    return passed;
}

int test_pauliX(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        pauliX(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_pauliY(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        pauliY(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_pauliZ(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        pauliZ(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_hadamard(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        hadamard(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_sGate(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        sGate(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_tGate(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        tGate(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_phaseShift(char testName[200]) {
    
    qreal pi = 3.1415926535897932384626;

    int passed=1;
    int numQubits=3;
    Qureg mq;

    // prepare |00>(|0> + |1>)/sqrt(2)
    mq = createQureg(numQubits, env);
    initZeroState(mq);
    hadamard(mq, 0);
    
    // enter state |00> (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2)
    // coeff of |0>:  1/sqrt(2)
    // coeff of |1>: - (1 + i)/2
    phaseShift(mq, 0, pi * 5/4.0 );
    
    if (passed) passed = compareReals(getRealAmp(mq, 0), 1/sqrt(2), COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmp(mq, 0),         0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getRealAmp(mq, 1),    -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmp(mq, 1),    -1/2.0, COMPARE_PRECISION);
    
    destroyQureg(mq, env);

    // also test MPI version
    numQubits = 4;
    mq = createQureg(numQubits, env);
    
    // prepare state (|0> + |1>)/sqrt(2) |111>
    initZeroState(mq);
    for (int i=0; i < 3; i++)
        pauliX(mq, i);
    hadamard(mq, 3);
    
    // enter state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |111>
    // coef of |1111> is - (1 + i)/2
    // index of |1111> is 2^4 - 1 = 15
    phaseShift(mq, 0, pi * 5/4.0 );
    if (passed) passed = compareReals(getRealAmp(mq, 15), -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmp(mq, 15), -1/2.0, COMPARE_PRECISION);
    
    
    destroyQureg(mq, env);
    
    return passed;
}






int test_controlledPhaseShift(char testName[200]) {
    int passed=1;
    
    qreal pi = 3.1415926535897932384626;
    
    Qureg mq;
    mq = createQureg(4, env);
    
    // prepare state (|0> + |1>)/sqrt(2) |010>
    initZeroState(mq);
    pauliX(mq, 1);
    hadamard(mq, 3);
    
    // confirm controlling first and third qubits does nothing (state |1010> = 2^1 + 2^3 = 10)
    controlledPhaseShift(mq, 0, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmp(mq, 10), 1/sqrt(2), COMPARE_PRECISION);
    controlledPhaseShift(mq, 2, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmp(mq, 10), 1/sqrt(2), COMPARE_PRECISION);
    
    // controlling 2nd qubit enters state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |010>
    controlledPhaseShift(mq, 1, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmp(mq, 10), -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmp(mq, 10), -1/2.0, COMPARE_PRECISION);
    
    // enter (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |011>
    pauliX(mq, 0);
    
    // enter (|0> + |1>)/sqrt(2) |011> where |0011> = 2^0 + 2^1 = 3
    controlledPhaseShift(mq, 0, 3, - pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmp(mq, 3), 1/sqrt(2), COMPARE_PRECISION);
    
    destroyQureg(mq, env);
    
    return passed;
}

int test_multiControlledPhaseShift(char testName[200]) {
    int passed=1;
    
    qreal pi = 3.1415926535897932384626;
    
    Qureg mq;
    mq = createQureg(4, env);
    
    // prepare state (|0> + |1>)/sqrt(2) |010>
    initZeroState(mq);
    pauliX(mq, 1);
    hadamard(mq, 3);
    
    // confirm controlling on 2nd,3rd,4th qubits does nothing (state |1010> = 2^1 + 2^3 = 10)
    int ctrls[] = {1,2,3};
    multiControlledPhaseShift(mq, ctrls , 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmp(mq, 10), 1/sqrt(2), COMPARE_PRECISION);
    
    // enter state (|0> + |1>)/sqrt(2) |110>
    pauliX(mq, 2);
    
    // controlling on 2nd,3rd,4th qubits enters state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |110>
    // index of state |1110> = 2^1 + 2^2 + 2^3 = 14
    multiControlledPhaseShift(mq, ctrls, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmp(mq, 14), -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmp(mq, 14), -1/2.0, COMPARE_PRECISION);
    
    destroyQureg(mq, env);
    
    return passed;
}

int test_controlledNot(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit, controlQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int j=0; j<3; j++){
        controlQubit=j;
        for (int i=0; i<3; i++){
            if (i==j) {count++; continue;}
            syncQuESTEnv(env);
            initStateDebug(mq);
            rotateQubit=i;
            controlledNot(mq, controlQubit, rotateQubit);

            sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
            initStateFromSingleFile(&mqVerif, filename, env);

            if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
        }
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_controlledPhaseFlip(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit, controlQubit;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int j=0; j<3; j++){
        controlQubit=j;
        for (int i=0; i<3; i++){
            if (i==j) {count++; continue;}
            initStateDebug(mq);
            rotateQubit=i;
            controlledPhaseFlip(mq, rotateQubit, controlQubit);

            sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
            initStateFromSingleFile(&mqVerif, filename, env);

            if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
        }
    }
    
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);
    
    return passed;
}

int test_multiControlledPhaseFlip(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=4;
    Qureg mq, mqVerif; 

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    int qubits[4]={0,1,2,3};
    initStateDebug(mq);
    multiControlledPhaseFlip(mq, qubits, 4);

    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);

    passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

    /* Tyson 21 July */
    if (passed) {
        int i;
        
        // prepare state |111>(|0> - |1>)/sqrt(2)
        initZeroState(mqVerif);
        for (i=0; i < 4; i++)
            pauliX(mqVerif, i);
        hadamard(mqVerif, 3);
    
        // prepare state |111>(|0> + |1>)/sqrt(2)
        initZeroState(mq);
        for (i=0; i < 3; i++)
            pauliX(mq, i);
        hadamard(mq, 3);
        
        // and transition to |111>(|0> - |1>)/sqrt(2)
        multiControlledPhaseFlip(mq, qubits, 4);
    
        // compare these states
        passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_compactUnitary(char testName[200]){
    int passed=1;

    int numQubits=10;
    int rotQubit;
    Qureg mq, mqVerif; 

    qreal angs[3];
    Complex alpha, beta;

    angs[0]=1.2; angs[1]=-2.4; angs[2]=0.3;
    alpha.real = cos(angs[0]) * cos(angs[1]);
    alpha.imag = cos(angs[0]) * sin(angs[1]);
    beta.real  = sin(angs[0]) * cos(angs[2]);
    beta.imag  = sin(angs[0]) * sin(angs[2]);
    
    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    initStateDebug(mq);
    initStateDebug(mqVerif);
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        compactUnitary(mq, rotQubit, alpha, beta);
    }
    // note -- this is only checking if the state changed at all due to rotation,
    // not that it changed correctly
    if (passed) passed = !compareStates(mq, mqVerif, COMPARE_PRECISION);

    // Rotate back the other way and check we arrive back at the initial state
    // (conjugate transpose of the unitary)
    alpha.imag *= -1;
    beta.real  *= -1;
    beta.imag  *= -1;

    // (order of qubits operated upon doesn't matter)
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        compactUnitary(mq, rotQubit, alpha, beta);
    }

    // unitaries are relatively imprecise (10* needed for signle precision)
    if (passed) passed = compareStates(mq, mqVerif, 10*COMPARE_PRECISION);

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    // check for normalisation
    numQubits=25;
    mq = createQureg(numQubits, env);
    initPlusState(mq);
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        compactUnitary(mq, rotQubit, alpha, beta);
    }
    qreal outcome = calcTotalProb(mq);    
    if (passed) passed = compareReals(1.0, outcome, COMPARE_PRECISION);
    destroyQureg(mq, env);


    return passed;
}

int test_unitary(char testName[200]){
    int passed=1;

    int numQubits=10;
    int rotQubit;
    Qureg mq, mqVerif; 

    qreal angs[3];
    Complex alpha, beta;
    ComplexMatrix2 u, uDagger;

    angs[0]=1.2; angs[1]=-2.4; angs[2]=0.3;
    alpha.real = cos(angs[0]) * cos(angs[1]);
    alpha.imag = cos(angs[0]) * sin(angs[1]);
    beta.real  = sin(angs[0]) * cos(angs[2]);
    beta.imag  = sin(angs[0]) * sin(angs[2]);

    u.r0c0 = (Complex) {.real=alpha.real, .imag=alpha.imag};
    u.r0c1 = (Complex) {.real=-beta.real, .imag=beta.imag}; 
    u.r1c0 = (Complex) {.real=beta.real, .imag=beta.imag};
    u.r1c1 = (Complex) {.real=alpha.real, .imag=-alpha.imag};

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    initStateDebug(mq);
    initStateDebug(mqVerif);
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        compactUnitary(mqVerif, rotQubit, alpha, beta);
        unitary(mq, rotQubit, u);
    }
    // assigning alpha/beta values to u such that compactUnitary should match unitary
    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

    // Rotate back the other way and check we arrive back at the initial state
    uDagger.r0c0.real=u.r0c0.real; uDagger.r0c0.imag=-u.r0c0.imag;
    uDagger.r0c1.real=u.r1c0.real; uDagger.r0c1.imag=-u.r1c0.imag;
    uDagger.r1c0.real=u.r0c1.real; uDagger.r1c0.imag=-u.r0c1.imag; 
    uDagger.r1c1.real=u.r1c1.real; uDagger.r1c1.imag=-u.r1c1.imag;

    for (int i=numQubits-1; i>=0; i--){
        rotQubit=i;
        unitary(mq, rotQubit, uDagger);
    }

    initStateDebug(mqVerif);
    
    // unitaries are relatively imprecise (10* needed for single precision)
    if (passed) passed = compareStates(mq, mqVerif, 10*COMPARE_PRECISION);

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    // check for normalisation
    numQubits = 25;
    mq = createQureg(numQubits, env);
    initPlusState(mq);
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        unitary(mq, rotQubit, uDagger);
    }
    qreal outcome = calcTotalProb(mq);    
    if (passed) passed = compareReals(1.0, outcome, COMPARE_PRECISION);
    destroyQureg(mq, env);

    return passed;
}

int test_controlledCompactUnitary(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotQubit, controlQubit;
    Qureg mq, mqVerif; 

    // assumes compactUnitary function is correct

    qreal ang1, ang2, ang3;
    ang1 = 1.2320;
    ang2 = 0.4230;
    ang3 = -0.65230;

    Complex alpha, beta;
    alpha.real = cos(ang1) * cos(ang2);
    alpha.imag = cos(ang1) * sin(ang2);
    beta.real  = sin(ang1) * cos(ang3);
    beta.imag  = sin(ang1) * sin(ang3);

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int j=0; j<3; j++){
        controlQubit=j;
        for (int i=0; i<3; i++){
            if (i==j){count++; continue;}
            initStateDebug(mq);
            rotQubit=i;
            controlledCompactUnitary(mq, controlQubit, rotQubit, alpha, beta);

            sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
            initStateFromSingleFile(&mqVerif, filename, env);

            if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
        }
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_controlledUnitary(char testName[200]){
    int passed=1;

    int numQubits=10;
    int rotQubit, controlQubit;

    ComplexMatrix2 u;
    Qureg mq, mqVerif; 

    // assumes controlledCompactUnitary function is correct

    qreal ang1, ang2, ang3;
    ang1 = 1.2320;
    ang2 = 0.4230;
    ang3 = -0.65230;

    Complex alpha, beta;
    alpha.real = cos(ang1) * cos(ang2);
    alpha.imag = cos(ang1) * sin(ang2);
    beta.real  = sin(ang1) * cos(ang3);
    beta.imag  = sin(ang1) * sin(ang3);

    u.r0c0 = (Complex) {.real=alpha.real, .imag=alpha.imag};
    u.r0c1 = (Complex) {.real=-beta.real, .imag=beta.imag}; 
    u.r1c0 = (Complex) {.real=beta.real, .imag=beta.imag};
    u.r1c1 = (Complex) {.real=alpha.real, .imag=-alpha.imag};

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    for (int j=0; j<numQubits; j++){
        controlQubit=j;
        for (int i=0; i<numQubits; i++){
            if (j==i) continue;
            initStateDebug(mq);
            initStateDebug(mqVerif);
            rotQubit=i;
            controlledCompactUnitary(mqVerif, controlQubit, rotQubit, alpha, beta);
            controlledUnitary(mq, controlQubit, rotQubit, u);

            if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
        }
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_multiControlledUnitary(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=10;
    int rotQubit, controlQubit;
    ComplexMatrix2 u;
    Qureg mq, mqVerif; 

    // assumes controlledCompactUnitary function is correct

    qreal ang1, ang2, ang3;
    ang1 = 1.2320;
    ang2 = 0.4230;
    ang3 = -0.65230;

    Complex alpha, beta;
    alpha.real = cos(ang1) * cos(ang2);
    alpha.imag = cos(ang1) * sin(ang2);
    beta.real  = sin(ang1) * cos(ang3);
    beta.imag  = sin(ang1) * sin(ang3);

    u.r0c0 = (Complex) {.real=alpha.real, .imag=alpha.imag};
    u.r0c1 = (Complex) {.real=-beta.real, .imag=beta.imag}; 
    u.r1c0 = (Complex) {.real=beta.real, .imag=beta.imag};
    u.r1c1 = (Complex) {.real=alpha.real, .imag=-alpha.imag};

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    // test mask contains one control qubit
    for (int j=0; j<numQubits; j++){
        controlQubit=j;
        for (int i=0; i<numQubits; i++){
            if (j==i) continue;
            initStateDebug(mq);
            initStateDebug(mqVerif);
            rotQubit=i;
            controlledCompactUnitary(mqVerif, controlQubit, rotQubit, alpha, beta);
            multiControlledUnitary(mq, &controlQubit, 1, rotQubit, u);

            if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
        }
    }

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    // randomly test a few different other multi control qubit masks 
    numQubits=4;
    int controlQubits[3];
    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    rotQubit=3;
    controlQubits[0]=0;
    controlQubits[1]=2;

    initStateDebug(mq);
    multiControlledUnitary(mq, controlQubits, 2, rotQubit, u);
    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);
    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

    rotQubit=1;
    controlQubits[0]=0;
    controlQubits[1]=2;
    controlQubits[2]=3;

    initStateDebug(mq);
    multiControlledUnitary(mq, controlQubits, 3, rotQubit, u);
    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);
    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

// @TODO add dens testing
int test_calcProbOfOutcome(char testName[200]){
    int passed=1;

    int numQubits=12;
    int qubit;
    qreal outcome;
    
    /*
     * state vector
     */

    Qureg mq; 
    mq = createQureg(numQubits, env);

    // test qubit = |0> 
    initZeroState(mq);
    for (qubit=0; qubit<numQubits; qubit++){
        outcome = calcProbOfOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);

        outcome = calcProbOfOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);
    }

    // test qubit = |1> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateOfSingleQubit(&mq, qubit, 1);
        outcome = calcProbOfOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);

        outcome = calcProbOfOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);
    }

    // test qubit = |+> 
    for (qubit=0; qubit<numQubits; qubit++){
        initPlusState(mq);
        outcome = calcProbOfOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);

        outcome = calcProbOfOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);
    }

    destroyQureg(mq, env);
    
   /*
    * density matrix
    */
    
    // @TODO


    return passed;
}

int test_collapseToOutcome(char testName[200]){
    int passed=1;

    int numQubits=3;
    int qubit;
    qreal prob;

    /*
     * state-vectors
     */
    
    Qureg mq, mqVerif;
    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    // test qubit = |0> 
    for (qubit=0; qubit<numQubits; qubit++){
        initZeroState(mq);
        initZeroState(mqVerif);
        prob = collapseToOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(1, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

        /* uncomment to test error is thrown
           initZeroState(&mq);
           initZeroState(&mqVerif);
           prob = collapseToOutcome(mq, qubit, 1);
           if (passed) passed = compareReals(0, prob, COMPARE_PRECISION);
           if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
           */
    }

    // test qubit = |1> 
    for (qubit=0; qubit<numQubits; qubit++){
        /* uncomment to test error is thrown
           initStateOfSingleQubit(&mq, qubit, 1);
           initStateOfSingleQubit(&mqVerif, qubit, 1);
           prob = collapseToOutcome(mq, qubit, 0);
           if (passed) passed = compareReals(0, prob, COMPARE_PRECISION);
           if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
           */

        initStateOfSingleQubit(&mq, qubit, 1);
        initStateOfSingleQubit(&mqVerif, qubit, 1);
        prob = collapseToOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(1, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    // test qubit = |+> 
    for (qubit=0; qubit<numQubits; qubit++){
        initPlusState(mq);
        initStateOfSingleQubit(&mqVerif, qubit, 0);
        prob = collapseToOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(0.5, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

        initPlusState(mq);
        initStateOfSingleQubit(&mqVerif, qubit, 1);
        prob = collapseToOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(0.5, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);
    
   /*
    * density matrices
    */
    
    Qureg reg, regVerif;
    reg = createDensityQureg(numQubits, env);
    regVerif = createDensityQureg(numQubits, env);
    
    // test pure state collapse correctly
    
    initZeroState(regVerif); // make |+...+>|0>
    for (qubit=1; qubit<numQubits; qubit++)
        hadamard(regVerif,qubit);
    
    initPlusState(reg); // make |+...+>|+>
    collapseToOutcome(reg, 0, 0); // collapse to |+...+>|0>
    if (passed) passed = compareStates(reg, regVerif, COMPARE_PRECISION);
    
    // test mixtures of orthogonal pure states collapse correctly (to a pure state)
    
    initZeroState(reg); // |0...0>
    initClassicalState(regVerif, (1<<numQubits)-1); // |1...1>
    addDensityMatrix(reg, .5, regVerif); // .5 |0><0| + .5 |1><1|
    collapseToOutcome(reg, 0, 1); // collapse to |1...1><1...1|
    if (passed) passed = compareStates(reg, regVerif, COMPARE_PRECISION);
    if (passed) passed = compareReals(calcProbOfOutcome(reg, 0, 1), 1, COMPARE_PRECISION);
    if (passed) passed = compareReals(calcProbOfOutcome(reg, 1, 1), 1, COMPARE_PRECISION); // non-measured qubit also collapsed
    if (passed) passed = compareReals(calcPurity(reg), 1, COMPARE_PRECISION);
    
    // test mixed states of non-orthogonal pure states collapse correctly
    
    // ... when measurement does nothing ...
    initZeroState(reg); // |00>|0>
    initZeroState(regVerif);
    hadamard(regVerif,1);
    hadamard(regVerif,2); // |++>|0>
    addDensityMatrix(reg, 1-0.3, regVerif); // 0.3|000><000| + 0.7|++0><++0|
    cloneQureg(regVerif, reg);
    collapseToOutcome(reg, 0, 0);
    if (passed) passed = compareStates(reg, regVerif, COMPARE_PRECISION);
    if (passed) passed = (calcPurity(reg) < 1.0); // we should still be mixed
    
    // ... when measurement kills a mixture and collapses a superposition (here producing a pure state)
    initZeroState(reg);      // |000><000|
    initPlusState(regVerif); // |+++><+++|
    addDensityMatrix(reg, 1-0.4, regVerif); // 0.4 |000><000| + 0.6 |+++><+++|
    collapseToOutcome(reg, 0, 1); // |++1><++1|
    initClassicalState(regVerif,1); hadamard(regVerif,1); hadamard(regVerif,2); // |++1><++1|
    if (passed) passed = compareStates(reg, regVerif, COMPARE_PRECISION);
    if (passed) passed = compareReals(calcPurity(reg), 1, COMPARE_PRECISION);
    
    // ... when measurement collapses superposition but preserves a mixture (pure states are orthogonal)

    initPlusState(reg); // |+++><+++|
    initClassicalState(regVerif, (1<<numQubits)-1); // |111><111|
    addDensityMatrix(reg, 1-0.1, regVerif); // 0.1 |+++><+++| + 0.9 |111><111|
    collapseToOutcome(reg, 0, 1); // 0.1 |++1><++1| + 0.9 |111><111|
    if (passed) passed = compareReals(calcProbOfOutcome(reg, 0, 1), 1, COMPARE_PRECISION);
    if (passed) passed = (calcPurity(reg) < 1.0); // we should still be mixed
    
    destroyQureg(reg, env);
    destroyQureg(regVerif, env);
    
    return passed;
}

int test_measure(char testName[200]){
    int passed=1;

    int numQubits=4;
    Qureg mq, mqVerif;
    int qubit;
    int outcome;

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    // test qubit = |0> 
    for (qubit=0; qubit<numQubits; qubit++){
        initZeroState(mq);
        initZeroState(mqVerif);
        outcome = measure(mq, qubit);
        if (passed) passed = (outcome==0);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    // test qubit = |1> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateOfSingleQubit(&mq, qubit, 1);
        initStateOfSingleQubit(&mqVerif, qubit, 1);
        outcome = measure(mq, qubit);
        if (passed) passed = (outcome==1);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    // visual check:
    // test qubit = |+> 
    int nTrials=10;
    unsigned long int seedArray[] = {18239, 12391};
    int numSeeds = 2;
    seedQuEST(seedArray, numSeeds);
    for (qubit=0; qubit<numQubits; qubit++){
        if (env.rank==0) printf("  %d trials: measure qubit %d when in state |+>:\n", nTrials, qubit);
        if (env.rank==0) printf("    value of qubit = [");
        for (int i=0; i<nTrials; i++){
            initPlusState(mq);
            outcome = measure(mq, qubit);
            if (env.rank==0) printf(" %d", outcome);
        }
        if (env.rank==0) printf("]\n");
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_measureWithStats(char testName[200]){
    int passed=1;

    int numQubits=4;
    Qureg mq, mqVerif;
    int qubit;
    int outcome;
    qreal prob;

    mq = createQureg(numQubits, env);
    mqVerif = createQureg(numQubits, env);

    // test qubit = |0> 
    for (qubit=0; qubit<numQubits; qubit++){
        initZeroState(mq);
        initZeroState(mqVerif);
        prob=0;
        outcome = measureWithStats(mq, qubit, &prob);
        if (passed) passed = (outcome==0);
        if (passed) passed = compareReals(prob, 1, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    // test qubit = |1> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateOfSingleQubit(&mq, qubit, 1);
        initStateOfSingleQubit(&mqVerif, qubit, 1);
        prob=0;
        outcome = measureWithStats(mq, qubit, &prob);
        if (passed) passed = compareReals(prob, 1, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    // test qubit = |+> 
    for (qubit=0; qubit<numQubits; qubit++){
        initPlusState(mq);
        prob=0;
        outcome = measureWithStats(mq, qubit, &prob);
        if (passed) passed = compareReals(prob, 0.5, COMPARE_PRECISION);
    }
    destroyQureg(mq, env);
    destroyQureg(mqVerif, env);

    return passed;
}

int test_getRealAmp(char testName[200]){
    int passed=1;

    int numQubits=5;
    qreal ampEl=0, ampElVerif=0;

    Qureg mq; 
    mq = createQureg(numQubits, env);
    initStateDebug(mq);

    for (int i=0; i<getNumAmps(mq); i++){
        ampElVerif = (i*2.0)/10.0;
        ampEl = getRealAmp(mq, i);
        if (passed) passed = (ampElVerif==ampEl);
    }
    destroyQureg(mq, env);

    return passed;
}

int test_getImagAmp(char testName[200]){
    int passed=1;

    int numQubits=5;
    qreal ampEl=0, ampElVerif=0;

    Qureg mq; 
    mq = createQureg(numQubits, env);

    initStateDebug(mq);
    for (int i=0; i<getNumAmps(mq); i++){
        ampElVerif = (i*2.0+1)/10.0;
        ampEl = getImagAmp(mq, i);
        if (passed) passed = (ampElVerif==ampEl);
    }
    destroyQureg(mq, env);

    return passed;
}

int test_getProbAmp(char testName[200]){
    int passed=1;

    int numQubits=5;
    qreal ampEl=0, ampElVerif=0;
    qreal realEl, imagEl;

    Qureg mq; 
    mq = createQureg(numQubits, env);

    initStateDebug(mq);
    for (int i=0; i<getNumAmps(mq); i++){
        realEl = (i*2.0)/10.0;
        imagEl = (i*2.0+1)/10.0;
        ampElVerif = realEl*realEl + imagEl*imagEl;
        ampEl = getProbAmp(mq, i);
        if (passed) passed = (ampElVerif==ampEl);
    }
    destroyQureg(mq, env);

    return passed;
}

int test_calcInnerProduct(char testName[200]) {
    
    qreal pi = 3.1415926535897932384626;

    int passed = 1;
    
    // create two 3-qubit pure states
    Qureg bra, ket;
    bra = createQureg(3, env);
    ket = createQureg(3, env);
    Complex prod;
    
    // when states are equal, <s|s> = 1
    initPlusState(bra);
    initPlusState(ket);
    prod = calcInnerProduct(bra, ket);
    if (passed) passed = compareReals(prod.real, 1, COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag, 0, COMPARE_PRECISION);
    
    // orthogonal states return <bra|ket> = 0
    initClassicalState(bra, 1);
    initClassicalState(ket, 2);
    prod = calcInnerProduct(bra, ket);
    if (passed) passed = compareReals(prod.real, 0, COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag, 0, COMPARE_PRECISION);
    
    // <000|+++> = 1/sqrt(2)^3
    initZeroState(bra);
    initPlusState(ket);
    prod = calcInnerProduct(bra, ket);
    if (passed) passed = compareReals(prod.real, pow(1/sqrt(2),3), COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag, 0, COMPARE_PRECISION);
    
    // test imag component is populated
    initClassicalState(bra, 1);         // <001|
    initZeroState(ket);                 // |000>
    hadamard(ket, 0);                   // |00+>
    phaseShift(ket, 0, pi * 5/4.0 );    // |a> = |00> (1/sqrt(2) |0> - 1/2 (1 + i) |1>)
    prod = calcInnerProduct(bra, ket);  // <001|a> = - 1/2 (1 + i)
    if (passed) passed = compareReals(prod.real, -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag, -1/2.0, COMPARE_PRECISION);
    
    // test bra has complex conjugated amps
    initClassicalState(ket, 1);         // |001>
    initZeroState(bra);                 // <000|
    hadamard(bra, 0);                   // <00+|
    phaseShift(bra, 0, pi * 5/4.0 );    // <a| = <00| (1/sqrt(2) <0| + 1/2 (i - 1) <1|)
    prod = calcInnerProduct(bra, ket);  // <a|001> = 1/2 (i - 1)
    if (passed) passed = compareReals(prod.real, -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag,  1/2.0, COMPARE_PRECISION);
    
    destroyQureg(bra, env);
    destroyQureg(ket, env);
    return passed;
}

int test_calcFidelity(char testName[200]) {
    int passed=1;
    int numQubits=5;
    qreal fid;
    
    Qureg pure;
    pure = createQureg(numQubits, env);
    
    /*
     * test pure fid = |<a|b>|^2 (trivially calcInnerProduct, so not rigorous)
     */
    Qureg otherPure;
    otherPure = createQureg(numQubits, env);
    
    initZeroState(pure);
    initPlusState(otherPure); // <0|+> = 1/sqrt(2^n)
    fid = calcFidelity(otherPure, pure); // |<0|+>|^2 = 1/2^n
    if (passed) passed = compareReals(fid, 1.0/pow(2.0,numQubits), COMPARE_PRECISION);
    
    destroyQureg(otherPure, env);
    
    /* 
     * test mixed fid = <a| b |a>
     */
    Qureg mixed;
    mixed = createDensityQureg(numQubits, env);
    
    // <0|0><0|0> = 1
    initZeroState(pure);
    initZeroState(mixed);
    fid = calcFidelity(mixed, pure); 
    if (passed) passed = compareReals(fid, 1.0, COMPARE_PRECISION);
    
    // <0|0...1><0...1|0> = 0
    initZeroState(pure);
    initClassicalState(mixed, 1);
    fid = calcFidelity(mixed, pure); 
    if (passed) passed = compareReals(fid, 0.0, COMPARE_PRECISION);
    
    // <111...|+><+|111...> = 1/2^n
    initClassicalState(pure, (1<<numQubits)-1);
    initPlusState(mixed);
    fid = calcFidelity(mixed, pure);
    if (passed) passed = compareReals(fid, 1/(qreal)(1<<numQubits), COMPARE_PRECISION);
    
    // <0| .2 |0><0| + .8 |..1><..1| |0> = .2
    Qureg otherMixed;
    otherMixed = createDensityQureg(numQubits, env);
    initZeroState(pure);
    initZeroState(mixed);
    initClassicalState(otherMixed, 1);
    addDensityMatrix(mixed, 0.8, otherMixed); // .2 |0><0| + .8 |..1><..1|
    fid = calcFidelity(mixed, pure); 
    if (passed) passed = compareReals(fid, 0.2, COMPARE_PRECISION);
        
    destroyQureg(otherMixed, env);
    destroyQureg(mixed, env);
    
    // finish
    destroyQureg(pure, env);
    return passed;
}

int test_addDensityMatrix(char testName[200]) {
    int passed=1;
    int numQubits=5;
    
    qreal prob;
    Qureg reg1, reg2;
    reg1 = createDensityQureg(numQubits, env);
    reg2 = createDensityQureg(numQubits, env);
    
    // prob_0( p1 |0...><0...| + (1-p1) |1...><1...| ) = p1
    qreal p1 = 0.3;
    initZeroState(reg1);
    initClassicalState(reg2, 1);
    addDensityMatrix(reg1, 1-p1, reg2);
    prob = calcProbOfOutcome(reg1, 0, 0);
    if (passed) passed = compareReals(prob, p1, COMPARE_PRECISION);
    
    // prob_0( p2 {p1 |0...><0...| + (1-p1) |1...><1...|} + (1-p2)|+><+| ) 
    // = p2 p1 + (1-p2) / sqrt(2)
    qreal p2 = 0.7;
    initPlusState(reg2);
    addDensityMatrix(reg1, 1-p2, reg2);
    prob = calcProbOfOutcome(reg1, 0, 0);
    qreal trueProb = p2*p1 + (1-p2)*0.5;
    if (passed) passed = compareReals(prob, trueProb, COMPARE_PRECISION);
    
    destroyQureg(reg1, env);
    destroyQureg(reg2, env);
    return passed;
}

int test_calcPurity(char testName[200]) {
    int passed=1;
    int numQubits=3;
    qreal purity;
    
    Qureg qureg;
    qureg = createDensityQureg(numQubits, env);
    
    // pure states are pure
    initPlusState(qureg);
    purity = calcPurity(qureg);
    if (passed) passed = compareReals(purity, 1, COMPARE_PRECISION);
    
    for (int state=0; state < 1<<numQubits; state++) {
        initClassicalState(qureg, state);
        if (passed) passed = compareReals(calcPurity(qureg), 1, COMPARE_PRECISION);
    }
    
    Qureg otherQureg;
    otherQureg = createDensityQureg(numQubits, env);
    
    // a|+><+| + b|+><+| = (a+b) |+><+| (pure)
    initPlusState(qureg);
    initPlusState(otherQureg);
    addDensityMatrix(qureg, 0.5, otherQureg);
    purity = calcPurity(qureg);
    if (passed) passed = compareReals(purity, 1, COMPARE_PRECISION);
    
    // mixture of orthogonal pure states (purity = p1^2 + (1-p1)^2)
    initClassicalState(qureg, 0);       // |0><0|
    initClassicalState(otherQureg, 1);  // |0...01><0...01|
    
    qreal p1 = 0.3;
    addDensityMatrix(qureg, 1-p1, otherQureg);    
    purity = calcPurity(qureg);
    if (passed) passed = compareReals(purity, p1*p1 + (1-p1)*(1-p1), COMPARE_PRECISION);
    
    // mixture of non-orthogonal pure states, where <a|b> = c.
    /* Let rho = p1 |a><a| + (1-p1) |b><b|
     * Then rho^2 =   p1^2|a><a|     + p1(1-p1)|a><a|b><b| 
     *              + (1-p1)^2|b><b| + p1(1-p1)|b><b|a><a|
     *            = p1^2 |a><a| + (1-p1)^2|b><b| + p1(1-p1)(c|a><b| + c* |b><a|)
     * By using that Tr(|b><a|) = <a|b> = c, and linearity, we have
     * Tr(rho^2) = p1^2 + (1-p1)^2 + p1(1-p1)( c c* + c* c)
     *           = p1^2 + (1-p1)^2 + 2p1(1-p1)|c|^2
     *           = d p1^2 - d p1 + 1,   where d = 2(1-|c|^2)
    */ 
    initZeroState(qureg);       // |0> |0...>
    initZeroState(otherQureg);  
    hadamard(otherQureg, 0);    // 1/sqrt(2)(|0> + |1>) |0...>
    
    // c = 1/sqrt(2), d = 2(1-1/2) = 1, Tr(rho^2) = p1^2 - p1 + 1
    addDensityMatrix(qureg, 1-p1, otherQureg);
    purity = calcPurity(qureg);
    if (passed) passed = compareReals(purity, p1*p1 - p1 + 1, COMPARE_PRECISION);
        
    destroyQureg(qureg, env);
    destroyQureg(otherQureg, env);
    return passed;
}

int test_calcTotalProb(char testName[200]) {
    int passed=1;
    int numQubits=3;
    qreal prob;
    
    /* 
     * state-vector
     */
    Qureg qureg;
    qureg = createQureg(numQubits, env);
    
    
    
    hadamard(qureg, 0);
    rotateY(qureg, 0, 0.1);
    rotateZ(qureg, 0, 0.4);
    controlledRotateY(qureg, 0, 1, 0.9);
    controlledRotateX(qureg, 1, 2, 1.45);
    int ctrls[] = {0,1,2};
    multiControlledPhaseFlip(qureg, ctrls, 3);
    controlledRotateAroundAxis(qureg, 1, 0, 0.3, (Vector) {.x=1,.y=2,.z=3});
    prob = calcTotalProb(qureg);
    if (passed) passed = compareReals(prob, 1, COMPARE_PRECISION);
    
    destroyQureg(qureg, env);
    
    /* 
     * density matrix
     */
    qureg = createDensityQureg(numQubits, env);

    hadamard(qureg, 0);
    rotateY(qureg, 0, 0.1);
    rotateZ(qureg, 0, 0.4);
    controlledRotateY(qureg, 0, 1, 0.9);
    controlledRotateX(qureg, 1, 2, 1.45); 
    multiControlledPhaseFlip(qureg, ctrls, 3);
    controlledRotateAroundAxis(qureg, 1, 0, 0.3, (Vector) {.x=1,.y=2,.z=3});
    prob = calcTotalProb(qureg);
        
    if (passed) passed = compareReals(prob, 1, COMPARE_PRECISION);
    
    destroyQureg(qureg, env);
    return passed;
}

int test_applyOneQubitDepolariseError(char testName[200]) {
    char filename[200];
    int passed=1;
    int numQubits=3;

    int targetQubit;
    qreal depolProb=0.375;
    
    Qureg qureg, quregVerif;
    qureg = createDensityQureg(numQubits, env);
    quregVerif = createDensityQureg(numQubits, env);

    for (int i=0; i<numQubits; i++){
        initStateDebug(qureg);
        targetQubit=i;
        applyOneQubitDepolariseError(qureg, targetQubit, depolProb);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, i);  
        initStateFromSingleFile(&quregVerif, filename, env);
        //reportStateToScreen(qureg, env, 0);
        //reportStateToScreen(quregVerif, env, 0);

        if (passed) passed = compareStates(qureg, quregVerif, COMPARE_PRECISION);
    }
    destroyQureg(qureg, env);
    destroyQureg(quregVerif, env);

    return passed;
}

int test_applyOneQubitDephaseError(char testName[200]) {
    char filename[200];
    int passed=1;
    int numQubits=3;

    int targetQubit;
    qreal dephaseProb=0.25;
    
    Qureg qureg, quregVerif;
    qureg = createDensityQureg(numQubits, env);
    quregVerif = createDensityQureg(numQubits, env);

    for (int i=0; i<numQubits; i++){
        initStateDebug(qureg);
        targetQubit=i;
        applyOneQubitDephaseError(qureg, targetQubit, dephaseProb);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, i);  
        initStateFromSingleFile(&quregVerif, filename, env);
        //reportStateToScreen(qureg, env, 0);
        //reportStateToScreen(quregVerif, env, 0);

        if (passed) passed = compareStates(qureg, quregVerif, COMPARE_PRECISION);
    }
    destroyQureg(qureg, env);
    destroyQureg(quregVerif, env);

    return passed;
}

int test_applyTwoQubitDepolariseError(char testName[200]) {
    if (env.rank==0) printf("NOTE: twoQubitDepolarise test currently assumes GPU version is correct and tests against that version's output\n");
    char filename[200];
    int passed=1;
    int numQubits=3;

    int qubit1, qubit2;
    qreal depolProb=0.46875;
    
    Qureg qureg, quregVerif;
    qureg = createDensityQureg(numQubits, env);
    quregVerif = createDensityQureg(numQubits, env);

    int count=0;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
        qubit1=i; qubit2=j;
        if (qubit1==qubit2) continue;

        initStateDebug(qureg);
        applyTwoQubitDepolariseError(qureg, qubit1, qubit2, depolProb);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&quregVerif, filename, env);
        //reportStateToScreen(qureg, env, 0);
        //reportStateToScreen(quregVerif, env, 0);

        if (passed) passed = compareStates(qureg, quregVerif, COMPARE_PRECISION);
        }
    }
    destroyQureg(qureg, env);
    destroyQureg(quregVerif, env);

    return passed;
}

int test_applyTwoQubitDephaseError(char testName[200]) {
    if (env.rank==0) printf("NOTE: twoQubitDephase test currently assumes GPU version is correct and tests against that version's output\n");
    char filename[200];
    int passed=1;
    int numQubits=3;

    int qubit1, qubit2;
    qreal dephaseProb=0.375;
    
    Qureg qureg, quregVerif;
    qureg = createDensityQureg(numQubits, env);
    quregVerif = createDensityQureg(numQubits, env);

    int count=0;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
        qubit1=i; qubit2=j;
        if (qubit1==qubit2) continue;

        initStateDebug(qureg);
        applyTwoQubitDephaseError(qureg, qubit1, qubit2, dephaseProb);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&quregVerif, filename, env);
        //reportStateToScreen(qureg, env, 0);
        //reportStateToScreen(quregVerif, env, 0);

        if (passed) passed = compareStates(qureg, quregVerif, COMPARE_PRECISION);
        }
    }
    destroyQureg(qureg, env);
    destroyQureg(quregVerif, env);

    return passed;
}



int main (int narg, char** varg) {
    env = createQuESTEnv();
    reportQuESTEnv(env);

    int (*tests[NUM_TESTS])(char[200]) = {
        test_controlledNot,
        test_initZeroState,
        test_initPlusState,
        test_initClassicalState,
        test_initPureState,
        test_setAmps,
        test_pauliX,
        test_pauliY,
        test_pauliZ,
        test_hadamard,
        test_sGate,
        test_tGate,
        test_phaseShift,
        test_controlledPhaseShift,
        test_multiControlledPhaseShift,
        test_controlledPhaseFlip,
        test_multiControlledPhaseFlip,
        test_compactUnitary,
        test_unitary,
        test_controlledCompactUnitary,
        test_controlledUnitary,
        test_multiControlledUnitary,
        test_calcProbOfOutcome,
        test_collapseToOutcome,
        test_measure,
        test_measureWithStats,
        test_getRealAmp,
        test_getImagAmp,
        test_getProbAmp,
        test_calcInnerProduct,
        test_calcFidelity,
        test_addDensityMatrix,
        test_calcPurity,
        test_calcTotalProb,
        test_applyOneQubitDephaseError,
        test_applyOneQubitDepolariseError,
        test_applyTwoQubitDephaseError,
        test_applyTwoQubitDepolariseError,
    };

    char testNames[NUM_TESTS][200] = {
        "controlledNot",
        "initZeroState",
        "initPlusState",
        "initClassicalState",
        "initPureState",
        "setAmps",
        "pauliX",
        "pauliY",
        "pauliZ",
        "hadamard",
        "sGate",
        "tGate",
        "phaseShift",
        "controlledPhaseShift",
        "multiControlledPhaseShift",
        "controlledPhaseFlip",
        "multiControlledPhaseFlip",
        "compactUnitary",
        "unitary",
        "controlledCompactUnitary",
        "controlledUnitary",
        "multiControlledUnitary",
        "calcProbOfOutcome",
        "collapseToOutcome",
        "measure",
        "measureWithStats",
        "getRealAmp",
        "getImagAmp",
        "getProbAmp",
        "calcInnerProduct",
        "calcFidelity",
        "addDensityMatrix",
        "calcPurity",
        "calcTotalProb",
        "applyOneQubitDephaseError",
        "applyOneQubitDepolariseError",
        "applyTwoQubitDephaseError",
        "applyTwoQubitDepolariseError",
    };
    int passed=0;
    if (env.rank==0) printf("\nRunning unit tests\n");
    for (int i=0; i<NUM_TESTS; i++){
        passed=(*tests[i])(testNames[i]);   
        passed=syncQuESTSuccess(passed);
        if (!passed){
            if (env.rank==0) printf("!!! FAILED in test %d -- %s\n", i, testNames[i]);
            destroyQuESTEnv(env);
            return 1;
        } else if (env.rank==0) printf("Passed test %d -- %s\n", i, testNames[i]);
    }
    destroyQuESTEnv(env);

    return 0;
}


