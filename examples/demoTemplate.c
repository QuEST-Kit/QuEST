// ==================================================================== //
//                                                                      //
//  demo.c -- qubit operarion demonstrator for QuEST                    //
//                                                                      //
// ==================================================================== //

# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>

# include "qubits.h"


# define MaxAngles      10
# define maxNumQubits   34
# define N_TRIALS 20
# define REPORT_STATE 0

const long double Pi = 3.14159265358979323846264338327950288419716939937510;


// ==================================================================== //
//                                                                      //
//     system_timer -- precision walltime function, computes            //
//                     walltime based on the gettimeofday() function    //
//                                                                      //
// ==================================================================== //

double system_timer (void) {

# include <stdlib.h>
# include <sys/time.h>

  struct timeval time;

  gettimeofday (&time, NULL);

  return time.tv_sec + time.tv_usec / 1000000.0;

}


//--------------------------------------------------------------
//---------------------- START OF main()  ----------------------
//--------------------------------------------------------------
int main (int narg, char** varg) {

  int numQubits;
  long long int index;
  double *stateVecReal;
  double *stateVecImag;

  double ang1,ang2,ang3;
  double alphaRe,alphaIm;
  double betaRe,betaIm;

  double stateProb,randProb;

  // timing variables
  double wtime_start,
         wtime_stop;
  double *timingVec;
  int trial;

  // define rotation angles
  double angles[MaxAngles][3] = {
    { 1.2320,  0.4230, -0.6523},
    { 2.1213,  0.0000,  3.6520},
    {-3.1213,  5.0230,  0.1230},
    { 5.2341, -3.1001, -1.2340},
    {-0.1234, -0.9876,  4.1234}
  };
  int numAngles=5,iAngle;

  int rotQubit,measureQubit;
  long long int numAmps;

  // number of qubits is command line argument
  if (narg == 2) {
    numQubits = atoi(varg[1]);
    if (numQubits < 1 || numQubits > maxNumQubits) {
      printf(" *** error: argument %d out of range (1 -- %d)\n", numQubits,maxNumQubits);
      exit (EXIT_FAILURE);
    }
  } else {
    printf(" *** error: too few arguments, number of qubits expected\n");
    exit (EXIT_FAILURE);
  }


  numAmps = 1L << numQubits;

  printf("Demo of single qubit rotations.\n");
  printf("Number of qubits is %d.\n", numQubits);
  printf("Number of amps is %lld.\n", numAmps);


  // allocate memory
  
  timingVec = (double*) malloc(N_TRIALS*numQubits*sizeof(double));
  stateVecReal = (double *) malloc(numAmps * sizeof(double));
  stateVecImag = (double *) malloc(numAmps * sizeof(double));
  if ((!stateVecReal || !stateVecImag) && numAmps) {
    printf("Could not allocate memory!");
    exit (EXIT_FAILURE);
  }
    
  // initialise the state to |0000..0>
  initStateVec (numQubits, stateVecReal,stateVecImag);


  //
  // ===== apply rotations
  //
  /* // keep time */
  /* wtime_start = system_timer (); */

  // rotate
  ang1 = angles[0][0];
  ang2 = angles[0][1];
  ang3 = angles[0][2];

  alphaRe = cos(ang1) * cos(ang2);
  alphaIm = cos(ang1) * sin(ang2);
  betaRe  = sin(ang1) * cos(ang3);
  betaIm  = sin(ang1) * sin(ang3);

  for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
  //for (rotQubit=0; rotQubit<10; rotQubit++) {
  for (trial=0; trial<N_TRIALS; trial++){
    wtime_start = system_timer ();
    rotateQubit(numQubits,rotQubit,alphaRe,alphaIm,betaRe,betaIm,stateVecReal,stateVecImag);
    wtime_stop = system_timer ();
    if (trial==0) printf(" rotation qubit %d: elapsed time = %f [s]\n", rotQubit, wtime_stop - wtime_start);
    timingVec[trial*numQubits + rotQubit]=wtime_stop-wtime_start;
  }
  }

  	char filename[100];
	FILE *timing;
        sprintf(filename, "timing.csv");
	timing = fopen(filename, "w");
	fprintf(timing, "qubit, time(s)\n");

	double totTime;
	for(index=0; index<numQubits; index++){
		totTime=0;
		for (trial=0; trial<N_TRIALS; trial++){
			totTime+=timingVec[trial*numQubits + index];
		}
		fprintf(timing, "%lld, %.8f\n", index, totTime/(double)N_TRIALS);
	}
	fclose(timing);

  /* // keep time */
  /* wtime_stop = system_timer (); */

  /* // ----- timing report */
  /* printf(" rotation: elapsed time = %f [s]\n", wtime_stop - wtime_start); */

  
  // ----- report vector
  if (numQubits < 7) {
    printf("\n\nThe following is the final state after rotations.\n\n");
    printf("codeOutput=[\n");
    for(index=0; index<=numAmps-1; index++) printf("%.8f %.8f\n",stateVecReal[index],stateVecImag[index]);
    printf("];\n\n");
  }

 if (REPORT_STATE){ 
	FILE *state;
	sprintf(filename, "state_orig.csv");
	state = fopen(filename, "w");
	fprintf(state, "real, imag\n");
	for(index=0; index<=numAmps-1; index++){
	       	fprintf(state, "%.12f, %.12f\n",stateVecReal[index],stateVecImag[index]);
	}
	fclose(state);
  }



  //
  // ===== perform a measurement
  //
  /* // keep time */
  /* wtime_start = system_timer (); */

  // measure
  /*
  for (measureQubit=0; measureQubit<numQubits; measureQubit++) {
    wtime_start = system_timer ();
    stateProb = findProbabilityOfZero (numQubits, measureQubit, stateVecReal,stateVecImag);
    wtime_stop = system_timer ();
    printf("   probability of 0 for qubit %d = %.14f\n", measureQubit, stateProb);
    printf(" measurement qubit %d: elapsed time = %f [s]\n", measureQubit, wtime_stop - wtime_start);
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

  // two qubit phase gate
  /*
  if (numQubits >= 7) {
    wtime_start = system_timer (); 
    controlPhaseGate (numQubits, 0, 2, stateVecReal, stateVecImag); 
    wtime_stop = system_timer (); 
    printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
    wtime_start = system_timer (); controlPhaseGate (numQubits, 1, 3, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start); 
    wtime_start = system_timer (); controlPhaseGate (numQubits, 2, 4, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
    wtime_start = system_timer (); controlPhaseGate (numQubits, 3, 5, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
    wtime_start = system_timer (); controlPhaseGate (numQubits, 4, 6, stateVecReal, stateVecImag); wtime_stop = system_timer (); printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start);
  }
  */
/*
    printf("\n\nThe following is the final state after rotations.\n\n");
    printf("codeOutput=[\n");
    for(index=0; index<=numAmps-1; index++) printf("%.8f %.8f\n",stateVecReal[index],stateVecImag[index]);
    printf("];\n\n");
*/

  /* // keep time */
  /* wtime_stop = system_timer (); */

  /* // ----- timing report */
  /* printf(" two qubit phase gate: elapsed time = %f [s]\n", wtime_stop - wtime_start); */


  // free memory
  free(stateVecReal);
  free(stateVecImag);
  free(timingVec);


  return EXIT_SUCCESS;
}

//--------------------------------------------------------------
//----------------------- END OF main()  -----------------------
//--------------------------------------------------------------
