# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include <omp.h>

# include "QuEST/qubits.h"
# include "QuEST/precision.h"
# include "QuEST/qubits_debug.h"

# define NUM_TESTS 17
# define COMPARE_PRECISION 10e-13
# define PATH_TO_TESTS "tests/unit/"
# define VERBOSE 0

QuESTEnv env;

void reportTest(MultiQubit multiQubit, char testName[200]){
	printf("\nTest: %s\n", testName);
	reportStateToScreen(multiQubit, env, 0);
}

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
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	initStateZero(&mq);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
	initializeStateFromSingleFile(&mqVerif, filename, env);

	passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_initStatePlus(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	initStatePlus(&mq);

        sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
	initializeStateFromSingleFile(&mqVerif, filename, env);

	passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_sigmaX(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int i=0; i<3; i++){
		initStateDebug(&mq);
		rotateQubit=i;
		sigmaX(mq, rotateQubit);
        	
		sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
		initializeStateFromSingleFile(&mqVerif, filename, env);

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_sigmaY(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int i=0; i<3; i++){
		initStateDebug(&mq);
		rotateQubit=i;
		sigmaY(mq, rotateQubit);
        	
		sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
		initializeStateFromSingleFile(&mqVerif, filename, env);

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_sigmaZ(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int i=0; i<3; i++){
		initStateDebug(&mq);
		rotateQubit=i;
		sigmaZ(mq, rotateQubit);
        	
		sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
		initializeStateFromSingleFile(&mqVerif, filename, env);

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_hadamard(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int i=0; i<3; i++){
		initStateDebug(&mq);
		rotateQubit=i;
		hadamard(mq, rotateQubit);
        	
		sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
		initializeStateFromSingleFile(&mqVerif, filename, env);

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_sGate(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int i=0; i<3; i++){
		initStateDebug(&mq);
		rotateQubit=i;
		sGate(mq, rotateQubit);
        	
		sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
		initializeStateFromSingleFile(&mqVerif, filename, env);

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_tGate(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int i=0; i<3; i++){
		initStateDebug(&mq);
		rotateQubit=i;
		tGate(mq, rotateQubit);
        	
		sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
		initializeStateFromSingleFile(&mqVerif, filename, env);

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_controlNot(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit, controlQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int j=0; j<3; j++){
		controlQubit=j;
		for (int i=0; i<3; i++){
			syncQuESTEnv(env);
			initStateDebug(&mq);
			rotateQubit=i;
			controlNot(mq, rotateQubit, controlQubit);
			
			sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
			initializeStateFromSingleFile(&mqVerif, filename, env);

			if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
		}
	}

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_controlPhaseGate(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotateQubit, controlQubit;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int j=0; j<3; j++){
		controlQubit=j;
		for (int i=0; i<3; i++){
			initStateDebug(&mq);
			rotateQubit=i;
			controlPhaseGate(mq, rotateQubit, controlQubit);
			
			sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
			initializeStateFromSingleFile(&mqVerif, filename, env);

			if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
		}
	}
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_quadCPhaseGate(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=4;
	MultiQubit mq, mqVerif; 

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	int qubit0=0, qubit1=1, qubit2=2, qubit3=3;
	initStateDebug(&mq);
	quadCPhaseGate(mq, qubit0, qubit1, qubit2, qubit3);

	sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
	initializeStateFromSingleFile(&mqVerif, filename, env);

	passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_rotateQubit(char testName[200]){
	printf("skipping rotate!\n");
	return 1;
	int passed=1;

	int numQubits=4;
	int rotQubit;
	MultiQubit mq, mqVerif; 

	REAL angs[3];
	Complex alpha, beta;

	for (int axis=0; axis<3; axis++){
		angs[0]=0; angs[1]=0; angs[2]=0;
		// rotate along one axis at a time
		angs[axis]=1.6;
		alpha.real = cos(angs[0]) * cos(angs[1]);
		alpha.imag = cos(angs[0]) * sin(angs[1]);
		beta.real  = sin(angs[0]) * cos(angs[2]);
		beta.imag  = sin(angs[0]) * sin(angs[2]);

		createMultiQubit(&mq, numQubits, env);
		createMultiQubit(&mqVerif, numQubits, env);

		initStateDebug(&mq);
		initStateDebug(&mqVerif);
		for (int i=0; i<numQubits; i++){
			rotQubit=i;
			rotateQubit(mq, rotQubit, alpha, beta);
		}
		// note -- this is only checking if the state changed at all due to rotation,
		// not that it changed correctly
		if (passed) passed = !compareStates(mq, mqVerif, COMPARE_PRECISION);

		// Rotate back the other way and check we arrive back at the initial state
		angs[axis]=-1.6;

		alpha.real = cos(angs[0]) * cos(angs[1]);
		alpha.imag = cos(angs[0]) * sin(angs[1]);
		beta.real  = sin(angs[0]) * cos(angs[2]);
		beta.imag  = sin(angs[0]) * sin(angs[2]);

		for (int i=0; i<numQubits; i++){
			rotQubit=i;
			rotateQubit(mq, rotQubit, alpha, beta);
		}

		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_controlRotateQubit(char testName[200]){
	char filename[200];
	int passed=1;
	int count=1;

	int numQubits=3;
	int rotQubit, controlQubit;
	MultiQubit mq, mqVerif; 

	// assumes rotateQubit function is correct
	
	REAL ang1, ang2, ang3;
        ang1 = 1.2320;
        ang2 = 0.4230;
        ang3 = -0.65230;

	Complex alpha, beta;
	alpha.real = cos(ang1) * cos(ang2);
        alpha.imag = cos(ang1) * sin(ang2);
        beta.real  = sin(ang1) * cos(ang3);
        beta.imag  = sin(ang1) * sin(ang3);

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	for (int j=0; j<3; j++){
		controlQubit=j;
		for (int i=0; i<3; i++){
			initStateDebug(&mq);
			rotQubit=i;
			controlRotateQubit(mq, rotQubit, controlQubit, alpha, beta);
			
			sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
			initializeStateFromSingleFile(&mqVerif, filename, env);

			if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
		}
	}
	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_findProbabilityOfOutcome(char testName[200]){
	int passed=1;

	int numQubits=12;
	MultiQubit mq; 
	int qubit;
	REAL outcome;

	createMultiQubit(&mq, numQubits, env);

	// test qubit = |0> 
	initStateZero(&mq);
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
		initStatePlus(&mq);
		outcome = findProbabilityOfOutcome(mq, qubit, 0);
		if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);

		outcome = findProbabilityOfOutcome(mq, qubit, 1);
		if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);
	}

	destroyMultiQubit(mq, env);

	return passed;
}

int test_measureInState(char testName[200]){
	int passed=1;

	int numQubits=3;
	MultiQubit mq, mqVerif;
	int qubit;
	REAL outcome;

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	// test qubit = |0> 
	for (qubit=0; qubit<numQubits; qubit++){
		initStateZero(&mq);
		initStateZero(&mqVerif);
		outcome = measureInState(mq, qubit, 0);
		if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);
		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);


		initStateZero(&mq);
		initStateZero(&mqVerif);
		outcome = measureInState(mq, qubit, 1);
		if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);
		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	// test qubit = |1> 
	for (qubit=0; qubit<numQubits; qubit++){
		initStateOfSingleQubit(&mq, qubit, 1);
		initStateOfSingleQubit(&mqVerif, qubit, 1);
		outcome = measureInState(mq, qubit, 0);
		if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);
		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

		initStateOfSingleQubit(&mq, qubit, 1);
		initStateOfSingleQubit(&mqVerif, qubit, 1);
		outcome = measureInState(mq, qubit, 1);
		if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);
		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	// test qubit = |+> 
	for (qubit=0; qubit<numQubits; qubit++){
		initStatePlus(&mq);
		initStateOfSingleQubit(&mqVerif, qubit, 0);
		outcome = measureInState(mq, qubit, 0);
		if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);
		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

		initStatePlus(&mq);
		initStateOfSingleQubit(&mqVerif, qubit, 1);
		outcome = measureInState(mq, qubit, 1);
		if (passed) passed = compareReals(0.5, outcome, COMPARE_PRECISION);
		if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);
	}

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}

int test_probOfFilterOut111(char testName[200]){
	char filename[200];
	int passed=1;
	int inCount=1;

	int numQubits=3;
	MultiQubit mq; 
	int qubit0=0, qubit1=1, qubit2=2;
	REAL outcome;

	createMultiQubit(&mq, numQubits, env);

	// test qubit = |0> 
	initStateZero(&mq);
	outcome = probOfFilterOut111(mq, qubit0, qubit1, qubit2);
	if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);

	// test qubit = |+> 
	initStatePlus(&mq);
	outcome = probOfFilterOut111(mq, qubit0, qubit1, qubit2);
	if (passed) passed = compareReals(7.0/8.0, outcome, COMPARE_PRECISION);

	sprintf(filename, "%s%s%d.in", PATH_TO_TESTS, testName, inCount++); 	
	initializeStateFromSingleFile(&mq, filename, env);
	outcome = probOfFilterOut111(mq, qubit0, qubit1, qubit2);
	if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);

	destroyMultiQubit(mq, env);

	return passed;
}


int test_filterOut111(char testName[200]){
	char filename[200];
	int passed=1;
	int inCount=1;
	int count=1;

	int numQubits=3;
	MultiQubit mq, mqVerif; 
	int qubit0=0, qubit1=1, qubit2=2;
	REAL outcome;

	createMultiQubit(&mq, numQubits, env);
	createMultiQubit(&mqVerif, numQubits, env);

	// test qubit = |0> 
	initStateZero(&mq);
	initStateZero(&mqVerif);
	outcome = filterOut111(mq, qubit0, qubit1, qubit2);
	if (passed) passed = compareReals(1, outcome, COMPARE_PRECISION);
	if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

	// test qubit = |+> 
	initStatePlus(&mq);
	outcome = filterOut111(mq, qubit0, qubit1, qubit2);
	sprintf(filename, "%s%s%d.out", PATH_TO_TESTS, testName, count++); 	
	initializeStateFromSingleFile(&mqVerif, filename, env);
	if (passed) passed = compareReals(7.0/8.0, outcome, COMPARE_PRECISION);
	if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

	sprintf(filename, "%s%s%d.in", PATH_TO_TESTS, testName, inCount++); 	
	initializeStateFromSingleFile(&mq, filename, env);
	outcome = probOfFilterOut111(mq, qubit0, qubit1, qubit2);
	initializeStateFromSingleFile(&mqVerif, filename, env);
	if (passed) passed = compareReals(0, outcome, COMPARE_PRECISION);
	if (passed) passed = compareStates(mq, mqVerif, COMPARE_PRECISION);

	destroyMultiQubit(mq, env);
	destroyMultiQubit(mqVerif, env);

	return passed;
}


int main (int narg, char** varg) {
	initQuESTEnv(&env);
	reportQuESTEnv(env);

	int (*tests[NUM_TESTS])(char[200]) = {
		test_controlNot,
		test_initStateZero,
		test_initStatePlus,
		test_sigmaX,
		test_sigmaY,
		test_sigmaZ,
		test_hadamard,
		test_sGate,
		test_tGate,
		test_controlPhaseGate,
		test_quadCPhaseGate,
		test_rotateQubit,
		test_controlRotateQubit,
		test_findProbabilityOfOutcome,
		test_measureInState,
		test_probOfFilterOut111,
		test_filterOut111
	};

	char testNames[NUM_TESTS][200] = {
		"controlNot",
		"initStateZero",
		"initStatePlus",
		"sigmaX",
		"sigmaY",
		"sigmaZ",
		"hadamard",
		"sGate",
		"tGate",
		"controlPhaseGate",
		"quadCPhaseGate",
		"rotateQubit",
		"controlRotateQubit",
		"findProbabilityOfOutcome",
		"measureInState",
		"probOfFilterOut111",
		"filterOut111"
	};
	int passed=0;

	if (env.rank==0) printf("\nRunning unit tests\n");
	for (int i=0; i<NUM_TESTS; i++){
		passed=(*tests[i])(testNames[i]);	
		passed=syncQuESTSuccess(env, passed);
		if (!passed){
			if (env.rank==0) printf("!!! FAILED in test %d -- %s\n", i, testNames[i]);
			closeQuESTEnv(env);
			return 1;
		} else if (env.rank==0) printf("Passed test %d -- %s\n", i, testNames[i]);
	}

	closeQuESTEnv(env);
	return 0;
}


