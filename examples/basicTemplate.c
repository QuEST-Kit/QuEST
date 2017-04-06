/** @file 
 * Basic template for using the QUEST library. In general, leave the initialisation
 * and cleanup sections as they are and edit the rotations, measurement and phase gate
 * sections.
 */

# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include <omp.h>

# include "QUEST/qubits.h"

//! Max number of angles used to define qubit rotation
# define MaxAngles      10
//! Max number of qubits in the system
# define maxNumQubits   40
//! 1: print end qubit state to file, 0: don't print
# define REPORT_STATE 1

const long double Pi = 3.14159265358979323846264338327950288419716939937510;


//--------------------------------------------------------------
//---------------------- START OF main()  ----------------------
//--------------------------------------------------------------
int main (int narg, char** varg) {

	//
	// ===== INITIALISATION
	//
	
	// INIT ENVIRONMENT: ALWAYS REQUIRED ONCE AT BEGINNING OF PROGRAM
	// These two lines will automatically set up the environment (multinode,
	// openMP only etc)  
	QUESTEnv env;
	initQUESTEnv(&env);

	// model vars
	int numQubits;
	long int index;
	long int numAmps;
	
	// get number of qubits from command line argument
	if (narg >= 2) {
		numQubits = atoi(varg[1]);
		if (numQubits < 1 || numQubits > maxNumQubits) {
			printf(" *** error: argument %d out of range (1 -- %d)\n", numQubits,maxNumQubits);
			exit (EXIT_FAILURE);
		}
	} else {
		printf(" *** error: too few arguments, number of qubits expected\n");
		exit (EXIT_FAILURE);
	}

	// Reporting
	numAmps = 1L << numQubits;
	if (env.rank==0){
		printf("Demo of single qubit rotations.\n");
		printf("Number of qubits is %d.\n", numQubits);
		printf("Number of amps is %ld.\n", numAmps);
	}

	// CREATE QUBIT OBJECT: REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// Before doing any operations on a set of qubits, create the MultiQubit object that will be used to 
	// represent the qubits
	MultiQubit multiQubit; 
	createMultiQubit(&multiQubit, numQubits, env);

	// initialise the state to |0000..0>
	initStateVec (&multiQubit);


	//
	// ===== ROTATIONS
	//

	// INITIALISE QUBIT ROTATION
	// Edit these lines to change rotation angle
	double ang1,ang2,ang3;
	Complex alpha, beta;

	// define rotation angles
	double angles[MaxAngles][3] = {
		{ 1.2320,  0.4230, -0.6523},
		{ 2.1213,  0.0000,  3.6520},
		{-3.1213,  5.0230,  0.1230},
		{ 5.2341, -3.1001, -1.2340},
		{-0.1234, -0.9876,  4.1234}
	};
	int numAngles=5,iAngle;

	// rotate
	ang1 = angles[0][0];
	ang2 = angles[0][1];
	ang3 = angles[0][2];

	alpha.real = cos(ang1) * cos(ang2);
	alpha.imag = cos(ang1) * sin(ang2);
	beta.real  = sin(ang1) * cos(ang3);
	beta.imag  = sin(ang1) * sin(ang3);

	int rotQubit;

	// DO QUBIT ROTATION
	// Edit these lines to perform rotations as required
	for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
		// do rotation of each qubit
		rotateQubit(multiQubit,rotQubit,alpha,beta);
		// make sure rotations have finished on the entire state vector, across all nodes
		// if running on only one node, this is a noop
		syncQUESTEnv(env);
	}
	// END QUBIT ROTATION

	// Verification: check vector size is unchanged
        double totalProbability;
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);

        // report state vector to file
	if (REPORT_STATE){
		reportState(multiQubit);
        }


	//
	// ===== perform a measurement
	//
	int measureQubit;
	double stateProb,randProb;
	/* // keep time */
	/* wtime_start = system_timer (); */

	// measure
/*	
	for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
	//for (measureQubit=0; measureQubit<1; measureQubit++) {
		syncQUESTEnv(env);
		wtime_start = system_timer ();
		stateProb = findProbabilityOfZero (env.rank, numAmpsPerRank, numQubits, measureQubit, stateVecReal,stateVecImag);
		syncQUESTEnv(env);
		wtime_stop = system_timer ();
		if (env.rank==0) printf("   probability of 0 for qubit %d = %.14f\n", measureQubit, stateProb);
		if (env.rank==0) printf(" measurement qubit %d: elapsed time = %f [s]\n", measureQubit, wtime_stop - wtime_start);
	}
*/
	/* // keep time */
	/* wtime_stop = system_timer (); */

	/* // ----- timing report */
	/* printf(" measurement: elapsed time = %f [s]\n", wtime_stop - wtime_start); */


	//
	// ===== two qubit phase gate
	//
	/* // keep time */
	/* wtime_start = system_timer (); */
/*
	// two qubit phase gate
	if (numQubits >= 7) {
		wtime_start = system_timer (); 
		controlPhaseGate (env.rank, numAmpsPerRank, numQubits, 0, 2, stateVecReal, stateVecImag); 
		wtime_stop = system_timer (); 
		printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
		wtime_start = system_timer (); controlPhaseGate (env.rank, numAmpsPerRank, numQubits, 1, 3, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
		wtime_start = system_timer (); controlPhaseGate (env.rank, numAmpsPerRank, numQubits, 2, 4, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
		wtime_start = system_timer (); controlPhaseGate (env.rank, numAmpsPerRank, numQubits, 3, 5, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
		wtime_start = system_timer (); controlPhaseGate (env.rank, numAmpsPerRank, numQubits, 4, 6, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
	}
*/
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);
/*
if (env.rank==0){
	printf("\n\nIn rank %d, the following is the final state after rotations.\n\n",env.rank);
	printf("codeOutput=[\n");
	for(index=0; index<=numAmpsPerRank-1; index++) printf("%.8f %.8f\n",stateVecReal[index],stateVecImag[index]);
	printf("];\n\n");
}
syncQUESTEnv(env);

if (env.rank==1){
	printf("\n\nIn rank %d, the following is the final state after rotations.\n\n",env.rank);
	printf("codeOutput=[\n");
	for(index=0; index<=numAmpsPerRank-1; index++) printf("%.8f %.8f\n",stateVecReal[index],stateVecImag[index]);
	printf("];\n\n");
}
syncQUESTEnv(env);
*/
	/* // keep time */
	/* wtime_stop = system_timer (); */

	/* // ----- timing report */
	/* printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start); */


	//
	// ======== CLEANUP
	//
	
	// free memory

	// REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// When all operations on a set of qubits are completed, destroy the object
	destroyMultiQubit(multiQubit, env);


	// ALWAYS REQUIRED ONCE AT END OF PROGRAM: 
	// These two lines will perform any necessary cleanup of the environment (multinode,
	// openMP only etc)  
	closeQUESTEnv(env);

	return EXIT_SUCCESS;
}

//--------------------------------------------------------------
//----------------------- END OF main()  -----------------------
//--------------------------------------------------------------
