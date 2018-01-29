/** @file 
 * Basic template for using the QuEST library. In general, leave the initialisation
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

# include "QuEST/qubits.h"
# include "QuEST/precision.h"

//! Max number of angles used to define qubit rotation
# define MaxAngles      10
//! Max number of qubits in the system
# define maxNumQubits   40
//! 1: print end qubit state to file, 0: don't print
# define REPORT_STATE 1

# define REPORT_RANK 0
# define MAX_ROT 5
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
	QuESTEnv env;
	initQuESTEnv(&env);

	// model vars
	int numQubits;
	int rotQubit;
	REAL totalProbability;

/*	
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
*/
	numQubits = 5;

	// CREATE QUBIT OBJECT: REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// Before doing any operations on a set of qubits, create the MultiQubit object that will be used to 
	// represent the qubits
	MultiQubit multiQubit; 
	createMultiQubit(&multiQubit, numQubits, env);
	
	// Reporting
	if (env.rank==0) {
		printf("Demo of single qubit rotations.\n");
	}
	reportMultiQubitParams(multiQubit);
	reportQuESTEnv(env);

	// initialise the state to |0000..0>
	if (env.rank==0) printf("|00000>\n");
	initStateZero(&multiQubit);
	reportStateToScreen(multiQubit, env, REPORT_RANK);

	// initialise the state to |0000..0>
	if (env.rank==0) printf("|+++++>\n");
	initStatePlus(&multiQubit);
	reportStateToScreen(multiQubit, env, REPORT_RANK);

	// initialise the state to debug mode (not physical)
	if (env.rank==0) printf("Init debug\n");	
	initStateDebug(&multiQubit);
	reportStateToScreen(multiQubit, env, REPORT_RANK);

	REAL el=0.0;
	el = getRealAmpEl(multiQubit, 0);
	if (env.rank==0) printf("el at index 0 = %.15f\n", el);
	el = getRealAmpEl(multiQubit, 1);
	if (env.rank==0) printf("el at index 1 = %.15f\n", el);
	if (multiQubit.numQubits>4){
		el = getRealAmpEl(multiQubit, 20);
		if (env.rank==0) printf("el at index 20 = %.15f\n", el);
	}


	//
	// ==== sigmaX, sigmaY, sigmaZ
	//
	if (env.rank==0) printf("\nPerforming sigmaX\n");
	initStateDebug(&multiQubit);
	for (rotQubit=0; rotQubit<MAX_ROT; rotQubit++) {
		// do rotation of each qubit
		sigmaX(multiQubit,rotQubit);
	}
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);

	if (env.rank==0) printf("\nPerforming sigmaY\n");
	initStateDebug(&multiQubit);
	for (rotQubit=0; rotQubit<MAX_ROT; rotQubit++) {
		// do rotation of each qubit
		sigmaY(multiQubit,rotQubit);
	}
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);

	if (env.rank==0) printf("\nPerforming sigmaZ\n");
	initStateDebug(&multiQubit);
	for (rotQubit=0; rotQubit<MAX_ROT; rotQubit++) {
		// do rotation of each qubit
		sigmaZ(multiQubit,rotQubit);
	}
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);
	
	//
	// ==== hadamard
	//
	
	if (env.rank==0) printf("\nPerforming hadamard\n");
	initStateDebug(&multiQubit);
	//for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
	for (rotQubit=0; rotQubit<MAX_ROT; rotQubit++) {
		// do rotation of each qubit
		hadamard(multiQubit,rotQubit);
	}
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);


	//
	// ==== s gate
	//
	
	if (env.rank==0) printf("\nPerforming s gate\n");
	initStateDebug(&multiQubit);
	//for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
	for (rotQubit=0; rotQubit<MAX_ROT; rotQubit++) {
		// do rotation of each qubit
		sGate(multiQubit,rotQubit);
	}
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);


	//
	// ==== t gate
	//
	
	if (env.rank==0) printf("\nPerforming t gate\n");
	initStateDebug(&multiQubit);
	//for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
	for (rotQubit=0; rotQubit<MAX_ROT; rotQubit++) {
		// do rotation of each qubit
		tGate(multiQubit,rotQubit);
	}
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);


	//
	// ==== control not
	//
	
	if (env.rank==0) printf("\nPerforming control not\n");
	initStateDebug(&multiQubit);
	if (env.rank==0) printf("target=0, control=0\n");
	rotQubit=0;
	controlledNot(multiQubit,rotQubit,0);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	initStateDebug(&multiQubit);
	if (env.rank==0) printf("target=2, control=0\n");
	rotQubit=2;
	controlledNot(multiQubit,rotQubit,0);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);

	//
	// ==== control rotate
	//
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
        ang1 = angles[0][0];
        ang2 = angles[0][1];
        ang3 = angles[0][2];

        alpha.real = cos(ang1) * cos(ang2);
        alpha.imag = cos(ang1) * sin(ang2);
        beta.real  = sin(ang1) * cos(ang3);
        beta.imag  = sin(ang1) * sin(ang3);
	
	if (env.rank==0) printf("\nPerforming ordinary rotate for comparison to control rotate\n");
	initStateDebug(&multiQubit);
	rotQubit=3;
	compactUnitary(multiQubit,rotQubit,alpha,beta);
	reportStateToScreen(multiQubit, env, REPORT_RANK);

	if (env.rank==0) printf("\nPerforming control rotate\n");
	initStateDebug(&multiQubit);
	if (env.rank==0) printf("target=3, control=3\n");
	rotQubit=3;
	controlledCompactUnitary(multiQubit,rotQubit,3,alpha,beta);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	initStateDebug(&multiQubit);
	if (env.rank==0) printf("target=3, control=0\n");
	controlledCompactUnitary(multiQubit,rotQubit,0,alpha,beta);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);


	//
	// ====== Measure
	//
	if (env.rank==0) printf("Performing measurement on state |00000>\n");
	initStateZero(&multiQubit);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
	int measureQubit;	
	REAL qProbability;
        for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
                qProbability = findProbabilityOfOutcome(multiQubit, measureQubit, 0);
                if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f\n", measureQubit, qProbability);
        }
        for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
                qProbability = findProbabilityOfOutcome(multiQubit, measureQubit, 1);
                if (env.rank==0) printf("Probability of 1 for qubit %d = %.14f\n", measureQubit, qProbability);
        }
	
	if (env.rank==0) printf("Performing measurement on state |+++++>\n");
	initStatePlus(&multiQubit);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
        for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
                qProbability = findProbabilityOfOutcome(multiQubit, measureQubit, 0);
                if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f\n", measureQubit, qProbability);
        }
        for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
                qProbability = findProbabilityOfOutcome(multiQubit, measureQubit, 1);
                if (env.rank==0) printf("Probability of 1 for qubit %d = %.14f\n", measureQubit, qProbability);
        }

        if (env.rank==0) printf("Measuring probability of qubit 0 to be in state 0 and then setting that qubit to 0\n");
	measureQubit=0;
        qProbability = collapseToOutcome(multiQubit, measureQubit, 0);
        if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f\n", measureQubit, qProbability);
	reportStateToScreen(multiQubit, env, REPORT_RANK);

        if (env.rank==0) printf("Performing single qubit measurement\n");
        // Do measurement on state with q0=0
        for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
                qProbability = findProbabilityOfOutcome(multiQubit, measureQubit, 0);
                if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f\n", measureQubit, qProbability);
        }
        totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);	


	if (env.rank==0) printf("|+++++>\n");
	initStatePlus(&multiQubit);
	reportStateToScreen(multiQubit, env, REPORT_RANK);
        if (env.rank==0) printf("Measuring probability of qubit 0 to be in state 1 and then setting that qubit to 1\n");
	measureQubit=0;
        qProbability = collapseToOutcome(multiQubit, measureQubit, 1);
        if (env.rank==0) printf("Probability of 1 for qubit %d = %.14f\n", measureQubit, qProbability);
	reportStateToScreen(multiQubit, env, REPORT_RANK);

        if (env.rank==0) printf("Performing single qubit measurement\n");
        // Do measurement on state with q0=1
        for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
                qProbability = findProbabilityOfOutcome(multiQubit, measureQubit, 1);
                if (env.rank==0) printf("Probability of 1 for qubit %d = %.14f\n", measureQubit, qProbability);
        }
        totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);	


        // report state vector to file
	if (REPORT_STATE){
		reportState(multiQubit);
        }

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
	closeQuESTEnv(env);

	return EXIT_SUCCESS;
}

//--------------------------------------------------------------
//----------------------- END OF main()  -----------------------
//--------------------------------------------------------------
