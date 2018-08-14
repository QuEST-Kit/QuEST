# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>

# include "QuEST.h"
# include "QuEST_precision.h"
# include "QuEST_debug.h"

# define NUM_TESTS 31
# define PATH_TO_TESTS "unit/"
# define VERBOSE 0

// quad precision unit testing no more stringent than double
# if QuEST_PREC==1
# define COMPARE_PRECISION 10e-5
# else
# define COMPARE_PRECISION 10e-13
# endif


QuESTEnv env;

void reportTest(QubitRegister qureg, char testName[200]){
    printf("\nTest: %s\n", testName);
    reportStateToScreen(qureg, env, 0);
}

/** returns 1 if equal to within precision */
int compareReals(REAL a, REAL b, REAL precision){
    REAL diff = a-b;
    if (diff<0) diff *= -1;
    if (diff>precision) return 0;
    else return 1;
}

int test_initStateZero(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    initStateZero(mq);

    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);

    passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_initStatePlus(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    initStatePlus(mq);

    sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
    initStateFromSingleFile(&mqVerif, filename, env);

    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_initClassicalState(char testName[200]){
    int passed=1;
    int numQubits=3;
    int numAmps=1 << numQubits;
    
    QubitRegister mq;
    createQubitRegister(&mq, numQubits, env);

    // test every classical state
    for (long long int stateInd=0LL; stateInd < numAmps; stateInd++) {
        initClassicalState(mq, stateInd);
        
        // check that every other state has prob 0
        for (long long int i=0LL; i < numAmps; i++) {
            if (i == stateInd)
                passed = passed && (getProbEl(mq,i) == 1.0);
            else
                passed = passed && (getProbEl(mq,i) == 0.0);
        }
    }

    destroyQubitRegister(mq, env);
    return passed;
}

int test_initPureState(char testName[200]) {
    
    // remember, this only tests initPureState for statevecs (not density matrices)
    
    int passed=1;
    int numQubits=3;
    
    QubitRegister mq, mqVerif;
    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);
    
    initStateZero(mq);
    initStateDebug(mqVerif);
    
    initPureState(mq, mqVerif);
    if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    
    return passed;
}

int test_initStateFromAmps(char testName[200]) {
    
    int passed=1;
    int numQubits=3;
    
    QubitRegister qureg;
    createQubitRegister(&qureg, numQubits, env);
    
    // test writing total state vec
    REAL reals[8] = {1,2,3,4,5,6,7,8};
    REAL imags[8] = {8,7,6,5,4,3,2,1};
    initStateFromAmps(qureg, 0, reals, imags, qureg.numAmpsTotal);
    for (long long int i=0; i < 8; i++) {
        if (passed) passed = compareReals(getRealAmpEl(qureg,i), reals[i], 0);
        if (passed) passed = compareReals(getImagAmpEl(qureg,i), imags[i], 0);
    }
    
    // test writing only some of statevec
    initStateZero(qureg);
    initStateFromAmps(qureg, 2, reals+2, imags+2, 4); // write {3,4,5,6} to inds {2,3,4,5}
    for (long long int i=0; i < 8; i++) {
        
        // indices outside {2,3,4,5} are unchanged from |0> = {1,0,0,0}...
        if (i==0) {
            if (passed) passed = compareReals(getRealAmpEl(qureg,i), 1, 0);
            if (passed) passed = compareReals(getImagAmpEl(qureg,i), 0, 0);
        }
        else if (i<2 || i>=6) {
            if (passed) passed = compareReals(getRealAmpEl(qureg,i), 0, 0);
            if (passed) passed = compareReals(getImagAmpEl(qureg,i), 0, 0);
        }
        // otherwise, indices should be set to those in reals and imag
        else {
            if (passed) passed = compareReals(getRealAmpEl(qureg,i), reals[i], 0);
            if (passed) passed = compareReals(getImagAmpEl(qureg,i), imags[i], 0);
        }
    }
    
    destroyQubitRegister(qureg, env);
    return passed;
}

int test_sigmaX(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        sigmaX(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_sigmaY(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        sigmaY(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_sigmaZ(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        sigmaZ(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_hadamard(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        hadamard(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_sGate(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        sGate(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_tGate(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    for (int i=0; i<3; i++){
        initStateDebug(mq);
        rotateQubit=i;
        tGate(mq, rotateQubit);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++);  
        initStateFromSingleFile(&mqVerif, filename, env);

        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_phaseShift(char testName[200]) {
    int passed=1;

    // prepare (|0> + |1>)/sqrt(2)
    QubitRegister mq;
    createQubitRegister(&mq, 1, env);
    initStatePlus(mq);
    
    // enter state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2)
    // coeff of |0>:  1/sqrt(2)
    // coeff of |1>: - (1 + i)/2
    REAL pi = 3.1415926535897932384626;
    phaseShift(mq, 0, pi * 5/4.0 );
    
    if (passed) passed = compareReals(getRealAmpEl(mq, 0), 1/sqrt(2), COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmpEl(mq, 0),         0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getRealAmpEl(mq, 1),    -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmpEl(mq, 1),    -1/2.0, COMPARE_PRECISION);
    destroyQubitRegister(mq, env);
    
    // also test MPI version
    createQubitRegister(&mq, 4, env);
    
    // prepare state (|0> + |1>)/sqrt(2) |111>
    initStateZero(mq);
    for (int i=0; i < 3; i++)
        sigmaX(mq, i);
    hadamard(mq, 3);
    
    // enter state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |111>
    // coef of |1111> is - (1 + i)/2
    // index of |1111> is 2^4 - 1 = 15
    phaseShift(mq, 0, pi * 5/4.0 );
    if (passed) passed = compareReals(getRealAmpEl(mq, 15), -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmpEl(mq, 15), -1/2.0, COMPARE_PRECISION);
    
    return passed;
}






int test_controlledPhaseShift(char testName[200]) {
    int passed=1;
    
    REAL pi = 3.1415926535897932384626;
    
    QubitRegister mq;
    createQubitRegister(&mq, 4, env);
    
    // prepare state (|0> + |1>)/sqrt(2) |010>
    initStateZero(mq);
    sigmaX(mq, 1);
    hadamard(mq, 3);
    
    // confirm controlling first and third qubits does nothing (state |1010> = 2^1 + 2^3 = 10)
    controlledPhaseShift(mq, 0, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmpEl(mq, 10), 1/sqrt(2), COMPARE_PRECISION);
    controlledPhaseShift(mq, 2, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmpEl(mq, 10), 1/sqrt(2), COMPARE_PRECISION);
    
    // controlling 2nd qubit enters state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |010>
    controlledPhaseShift(mq, 1, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmpEl(mq, 10), -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmpEl(mq, 10), -1/2.0, COMPARE_PRECISION);
    
    // enter (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |011>
    sigmaX(mq, 0);
    
    // enter (|0> + |1>)/sqrt(2) |011> where |0011> = 2^0 + 2^1 = 3
    controlledPhaseShift(mq, 0, 3, - pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmpEl(mq, 3), 1/sqrt(2), COMPARE_PRECISION);
    
    destroyQubitRegister(mq, env);
    
    return passed;
}

int test_multiControlledPhaseShift(char testName[200]) {
    int passed=1;
    
    REAL pi = 3.1415926535897932384626;
    
    QubitRegister mq;
    createQubitRegister(&mq, 4, env);
    
    // prepare state (|0> + |1>)/sqrt(2) |010>
    initStateZero(mq);
    sigmaX(mq, 1);
    hadamard(mq, 3);
    
    // confirm controlling on 2nd,3rd,4th qubits does nothing (state |1010> = 2^1 + 2^3 = 10)
    multiControlledPhaseShift(mq, (int[]) {1,2,3}, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmpEl(mq, 10), 1/sqrt(2), COMPARE_PRECISION);
    
    // enter state (|0> + |1>)/sqrt(2) |110>
    sigmaX(mq, 2);
    
    // controlling on 2nd,3rd,4th qubits enters state (|0> - 1/sqrt(2) (1 + i) |1>)/sqrt(2) |110>
    // index of state |1110> = 2^1 + 2^2 + 2^3 = 14
    multiControlledPhaseShift(mq, (int[]) {1,2,3}, 3, pi * 5/4.0);
    if (passed) passed = compareReals(getRealAmpEl(mq, 14), -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(getImagAmpEl(mq, 14), -1/2.0, COMPARE_PRECISION);
    
    destroyQubitRegister(mq, env);
    
    return passed;
}

int test_controlledNot(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit, controlQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_controlledPhaseFlip(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotateQubit, controlQubit;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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
    
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);
    
    return passed;
}

int test_multiControlledPhaseFlip(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=4;
    QubitRegister mq, mqVerif; 

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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
        initStateZero(mqVerif);
        for (i=0; i < 4; i++)
            sigmaX(mqVerif, i);
        hadamard(mqVerif, 3);
    
        // prepare state |111>(|0> + |1>)/sqrt(2)
        initStateZero(mq);
        for (i=0; i < 3; i++)
            sigmaX(mq, i);
        hadamard(mq, 3);
        
        // and transition to |111>(|0> - |1>)/sqrt(2)
        multiControlledPhaseFlip(mq, qubits, 4);
    
        // compare these states
        passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_compactUnitary(char testName[200]){
    int passed=1;

    int numQubits=10;
    int rotQubit;
    QubitRegister mq, mqVerif; 

    REAL angs[3];
    Complex alpha, beta;

    angs[0]=1.2; angs[1]=-2.4; angs[2]=0.3;
    alpha.real = cos(angs[0]) * cos(angs[1]);
    alpha.imag = cos(angs[0]) * sin(angs[1]);
    beta.real  = sin(angs[0]) * cos(angs[2]);
    beta.imag  = sin(angs[0]) * sin(angs[2]);
    
    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    // check for normalisation
    numQubits=25;
    createQubitRegister(&mq, numQubits, env);
    initStatePlus(mq);
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        compactUnitary(mq, rotQubit, alpha, beta);
    }
    REAL outcome = calcTotalProbability(mq);    
    if (passed) passed = compareReals(1.0, outcome, COMPARE_PRECISION);
    destroyQubitRegister(mq, env);


    return passed;
}

int test_unitary(char testName[200]){
    int passed=1;

    int numQubits=10;
    int rotQubit;
    QubitRegister mq, mqVerif; 

    REAL angs[3];
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

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    // check for normalisation
    numQubits = 25;
    createQubitRegister(&mq, numQubits, env);
    initStatePlus(mq);
    for (int i=0; i<numQubits; i++){
        rotQubit=i;
        unitary(mq, rotQubit, uDagger);
    }
    REAL outcome = calcTotalProbability(mq);    
    if (passed) passed = compareReals(1.0, outcome, COMPARE_PRECISION);
    destroyQubitRegister(mq, env);

    return passed;
}

int test_controlledCompactUnitary(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=3;
    int rotQubit, controlQubit;
    QubitRegister mq, mqVerif; 

    // assumes compactUnitary function is correct

    REAL ang1, ang2, ang3;
    ang1 = 1.2320;
    ang2 = 0.4230;
    ang3 = -0.65230;

    Complex alpha, beta;
    alpha.real = cos(ang1) * cos(ang2);
    alpha.imag = cos(ang1) * sin(ang2);
    beta.real  = sin(ang1) * cos(ang3);
    beta.imag  = sin(ang1) * sin(ang3);

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_controlledUnitary(char testName[200]){
    int passed=1;

    int numQubits=10;
    int rotQubit, controlQubit;

    ComplexMatrix2 u;
    QubitRegister mq, mqVerif; 

    // assumes controlledCompactUnitary function is correct

    REAL ang1, ang2, ang3;
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

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_multiControlledUnitary(char testName[200]){
    char filename[200];
    int passed=1;
    int count=1;

    int numQubits=10;
    int rotQubit, controlQubit;
    ComplexMatrix2 u;
    QubitRegister mq, mqVerif; 

    // assumes controlledCompactUnitary function is correct

    REAL ang1, ang2, ang3;
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

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    // randomly test a few different other multi control qubit masks 
    numQubits=4;
    int controlQubits[3];
    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

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

    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_findProbabilityOfOutcome(char testName[200]){
    int passed=1;

    int numQubits=12;
    QubitRegister mq; 
    int qubit;
    REAL outcome;

    createQubitRegister(&mq, numQubits, env);

    // test qubit = |0> 
    initStateZero(mq);
    for (qubit=0; qubit<numQubits; qubit++){
        outcome = findProbabilityOfOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);

        outcome = findProbabilityOfOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);
    }

    // test qubit = |1> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateOfSingleQubit(&mq, qubit, 1);
        outcome = findProbabilityOfOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);

        outcome = findProbabilityOfOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);
    }

    // test qubit = |+> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStatePlus(mq);
        outcome = findProbabilityOfOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);

        outcome = findProbabilityOfOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);
    }

    destroyQubitRegister(mq, env);

    return passed;
}

int test_collapseToOutcome(char testName[200]){
    int passed=1;

    int numQubits=3;
    QubitRegister mq, mqVerif;
    int qubit;
    REAL prob;

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    // test qubit = |0> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateZero(mq);
        initStateZero(mqVerif);
        prob = collapseToOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(1, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

        /* uncomment to test error is thrown
           initStateZero(&mq);
           initStateZero(&mqVerif);
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
        initStatePlus(mq);
        initStateOfSingleQubit(&mqVerif, qubit, 0);
        prob = collapseToOutcome(mq, qubit, 0);
        if (passed) passed = compareReals(0.5, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

        initStatePlus(mq);
        initStateOfSingleQubit(&mqVerif, qubit, 1);
        prob = collapseToOutcome(mq, qubit, 1);
        if (passed) passed = compareReals(0.5, prob, COMPARE_PRECISION);
        if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
    }
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_measure(char testName[200]){
    int passed=1;

    int numQubits=4;
    QubitRegister mq, mqVerif;
    int qubit;
    int outcome;

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    // test qubit = |0> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateZero(mq);
        initStateZero(mqVerif);
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
            initStatePlus(mq);
            outcome = measure(mq, qubit);
            if (env.rank==0) printf(" %d", outcome);
        }
        if (env.rank==0) printf("]\n");
    }
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_measureWithStats(char testName[200]){
    int passed=1;

    int numQubits=4;
    QubitRegister mq, mqVerif;
    int qubit;
    int outcome;
    REAL prob;

    createQubitRegister(&mq, numQubits, env);
    createQubitRegister(&mqVerif, numQubits, env);

    // test qubit = |0> 
    for (qubit=0; qubit<numQubits; qubit++){
        initStateZero(mq);
        initStateZero(mqVerif);
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
        initStatePlus(mq);
        prob=0;
        outcome = measureWithStats(mq, qubit, &prob);
        if (passed) passed = compareReals(prob, 0.5, COMPARE_PRECISION);
    }
    destroyQubitRegister(mq, env);
    destroyQubitRegister(mqVerif, env);

    return passed;
}

int test_getRealAmpEl(char testName[200]){
    int passed=1;

    int numQubits=5;
    REAL ampEl=0, ampElVerif=0;

    QubitRegister mq; 
    createQubitRegister(&mq, numQubits, env);
    initStateDebug(mq);

    for (int i=0; i<getNumAmps(mq); i++){
        ampElVerif = (i*2.0)/10.0;
        ampEl = getRealAmpEl(mq, i);
        if (passed) passed = (ampElVerif==ampEl);
    }
    destroyQubitRegister(mq, env);

    return passed;
}

int test_getImagAmpEl(char testName[200]){
    int passed=1;

    int numQubits=5;
    REAL ampEl=0, ampElVerif=0;

    QubitRegister mq; 
    createQubitRegister(&mq, numQubits, env);

    initStateDebug(mq);
    for (int i=0; i<getNumAmps(mq); i++){
        ampElVerif = (i*2.0+1)/10.0;
        ampEl = getImagAmpEl(mq, i);
        if (passed) passed = (ampElVerif==ampEl);
    }
    destroyQubitRegister(mq, env);

    return passed;
}

int test_getProbEl(char testName[200]){
    int passed=1;

    int numQubits=5;
    REAL ampEl=0, ampElVerif=0;
    REAL realEl, imagEl;

    QubitRegister mq; 
    createQubitRegister(&mq, numQubits, env);

    initStateDebug(mq);
    for (int i=0; i<getNumAmps(mq); i++){
        realEl = (i*2.0)/10.0;
        imagEl = (i*2.0+1)/10.0;
        ampElVerif = realEl*realEl + imagEl*imagEl;
        ampEl = getProbEl(mq, i);
        if (passed) passed = (ampElVerif==ampEl);
    }
    destroyQubitRegister(mq, env);

    return passed;
}

int test_calcInnerProduct(char testName[200]) {
    
    REAL pi = 3.1415926535897932384626;

    int passed = 1;
    
    // create two 3-qubit pure states
    QubitRegister bra, ket;
    createQubitRegister(&bra, 3, env);
    createQubitRegister(&ket, 3, env);
    Complex prod;
    
    // when states are equal, <s|s> = 1
    initStatePlus(bra);
    initStatePlus(ket);
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
    initStateZero(bra);
    initStatePlus(ket);
    prod = calcInnerProduct(bra, ket);
    if (passed) passed = compareReals(prod.real, pow(1/sqrt(2),3), COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag, 0, COMPARE_PRECISION);
    
    // test imag component is populated
    initClassicalState(bra, 1);         // <001|
    initStateZero(ket);                 // |000>
    hadamard(ket, 0);                   // |00+>
    phaseShift(ket, 0, pi * 5/4.0 );    // |a> = |00> (1/sqrt(2) |0> - 1/2 (1 + i) |1>)
    prod = calcInnerProduct(bra, ket);  // <001|a> = - 1/2 (1 + i)
    if (passed) passed = compareReals(prod.real, -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag, -1/2.0, COMPARE_PRECISION);
    
    // test bra has complex conjugated amps
    initClassicalState(ket, 1);         // |001>
    initStateZero(bra);                 // <000|
    hadamard(bra, 0);                   // <00+|
    phaseShift(bra, 0, pi * 5/4.0 );    // <a| = <00| (1/sqrt(2) <0| + 1/2 (i - 1) <1|)
    prod = calcInnerProduct(bra, ket);  // <a|001> = 1/2 (i - 1)
    if (passed) passed = compareReals(prod.real, -1/2.0, COMPARE_PRECISION);
    if (passed) passed = compareReals(prod.imag,  1/2.0, COMPARE_PRECISION);
    
    destroyQubitRegister(bra, env);
    destroyQubitRegister(ket, env);
    return passed;
}




int test_calcFidelity(char testName[200]) {
    int passed=1;
    int numQubits=5;
    REAL fid;
    
    QubitRegister pure;
    createQubitRegister(&pure, numQubits, env);
    
    /*
     * test pure fid = |<a|b>|^2 (trivially calcInnerProduct, so not rigorous)
     */
    QubitRegister otherPure;
    createQubitRegister(&otherPure, numQubits, env);
    
    initStateZero(pure);
    initStatePlus(otherPure); // <0|+> = 1/sqrt(2^n)
    fid = calcFidelity(otherPure, pure); // |<0|+>|^2 = 1/2^n
    if (passed) passed = compareReals(fid, 1.0/pow(2.0,numQubits), COMPARE_PRECISION);
    
    destroyQubitRegister(otherPure, env);
    
    /* 
     * test mixed fid = <a| b |a>
     */
    QubitRegister mixed;
    createDensityQubitRegister(&mixed, numQubits, env);
    
    // <0|0><0|0> = 1
    initStateZero(pure);
    initStateZero(mixed);
    fid = calcFidelity(mixed, pure); 
    if (passed) passed = compareReals(fid, 1.0, COMPARE_PRECISION);
    
    // <0|0...1><0...1|0> = 0
    initStateZero(pure);
    initClassicalState(mixed, 1);
    fid = calcFidelity(mixed, pure); 
    if (passed) passed = compareReals(fid, 0.0, COMPARE_PRECISION);
    
    // <0| .2 |0><0| + .8 |..1><..1| |0> = .2
    QubitRegister otherMixed;
    createDensityQubitRegister(&otherMixed, numQubits, env);
    initStateZero(pure);
    initStateZero(mixed);
    initClassicalState(otherMixed, 1);
    combineDensityMatrices(0.2, mixed, 0.8, otherMixed); // .2 |0><0| + .8 |..1><..1|
    fid = calcFidelity(mixed, pure); 
    //if (passed) passed = compareReals(fid, 0.2, COMPARE_PRECISION);
    printf("@TODO: calcFidelity needs combineDensityMatrices!\n");
    
    destroyQubitRegister(otherMixed, env);
    destroyQubitRegister(mixed, env);
    
    // finish
    destroyQubitRegister(pure, env);
    return passed;
}


int main (int narg, char** varg) {
    initQuESTEnv(&env);
    reportQuESTEnv(env);

    int (*tests[NUM_TESTS])(char[200]) = {
        test_controlledNot,
        test_initStateZero,
        test_initStatePlus,
        test_initClassicalState,
        test_initPureState,
        test_initStateFromAmps,
        test_sigmaX,
        test_sigmaY,
        test_sigmaZ,
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
        test_findProbabilityOfOutcome,
        test_collapseToOutcome,
        test_measure,
        test_measureWithStats,
        test_getRealAmpEl,
        test_getImagAmpEl,
        test_getProbEl,
        test_calcInnerProduct,
        test_calcFidelity
    };

    char testNames[NUM_TESTS][200] = {
        "controlledNot",
        "initStateZero",
        "initStatePlus",
        "initClassicalState",
        "initPureState",
        "initStateFromAmps",
        "sigmaX",
        "sigmaY",
        "sigmaZ",
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
        "findProbabilityOfOutcome",
        "collapseToOutcome",
        "measure",
        "measureWithStats",
        "getRealAmpEl",
        "getImagAmpEl",
        "getProbEl",
        "calcInnerProduct",
        "calcFidelity"
    };
    int passed=0;
    if (env.rank==0) printf("\nRunning unit tests\n");
    for (int i=0; i<NUM_TESTS; i++){
        passed=(*tests[i])(testNames[i]);   
        passed=syncQuESTSuccess(passed);
        if (!passed){
            if (env.rank==0) printf("!!! FAILED in test %d -- %s\n", i, testNames[i]);
            closeQuESTEnv(env);
            return 1;
        } else if (env.rank==0) printf("Passed test %d -- %s\n", i, testNames[i]);
    }
    closeQuESTEnv(env);

    return 0;
}


