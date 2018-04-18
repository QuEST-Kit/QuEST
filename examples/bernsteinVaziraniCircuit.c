/** @file 
 * Implements the Bernstien--Vazirani circuit
 */

# include <stdio.h>
# include <math.h>
# include <omp.h>

# include "qubits.h" 


int main (int narg, char** varg) {

	
	/* 	
	 * PREPARE QuEST
	 */
	
	// model parameters
	int numQubits = 9;
	int secretNum = pow(2,4) + 1;
	
	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	
	// create register; let zeroth qubit be ancilla
	MultiQubit multiQubit; 
	createMultiQubit(&multiQubit, numQubits, env);
	initStateZero(&multiQubit);


	/* 	
	 * APPLY ALGORITHM
	 */
	 
	// NOT the ancilla
	sigmaX(multiQubit, 0);
	
	// CNOT secretNum bits with ancilla
	int bits = secretNum;
	int bit;
	for (int qb=1; qb < numQubits; qb++) {
		bit = bits % 2;
		bits /= 2;
		if (bit)
			controlledNot(multiQubit, 0, qb);
	}
	
	
	/* 	
	 * VERIFY FINAL STATE
	 */
	 
	 // calculate prob of solution state
	double successProb = 1.0;
	bits = secretNum;
	for (int qb=1; qb < numQubits; qb++) {
		bit = bits % 2;
		bits /= 2;
		successProb *= findProbabilityOfOutcome(
			multiQubit, qb, bit);
	}

	printf("solution reached with probability ");
	printf("%f", successProb);
	printf("\n");
	

	/*
	 * FREE MEMORY
	 */
	
	destroyMultiQubit(multiQubit, env); 
	closeQuESTEnv(env);
	return 0;
}
