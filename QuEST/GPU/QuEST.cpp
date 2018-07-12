// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file QuEST.cpp
 * The core of the QuEST Library.
 */

# include <math.h> 
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include "QuEST_precision.h"
# include "QuEST.h"
# include "QuEST_internal.h"

# include "mt19937ar.h" // MT random number generation for random measurement seeding
# include <sys/param.h>
# include <unistd.h>

# include <sys/types.h> // include getpid
# include <sys/time.h>

# define DEBUG 0

const char* errorCodes[] = {
    "Success",                                              // 0
    "Invalid target qubit. Note qubits are zero indexed.",  // 1
    "Invalid control qubit. Note qubits are zero indexed.", // 2 
    "Control qubit cannot equal target qubit.",             // 3
    "Invalid number of control qubits",                     // 4
    "Invalid unitary matrix.",                              // 5
    "Invalid rotation arguments.",                          // 6
    "Invalid system size. Cannot print output for systems greater than 5 qubits.", // 7
    "Can't collapse to state with zero probability.", // 8
    "Invalid number of qubits.", // 9
    "Invalid measurement outcome -- must be either 0 or 1." // 10
};

#ifdef __cplusplus
extern "C" {
#endif

/** Print the current state vector of probability amplitudes for a set of qubits to file.
 * File format:
 * @verbatim
real, imag
realComponent1, imagComponent1
realComponent2, imagComponent2
...
realComponentN, imagComponentN
@endverbatim
 *
 * File naming convention:
 *
 * For each node that the program runs on, a file 'state_rank_[node_rank].csv' is generated. If there is
 * more than one node, ranks after the first do not include the header
 * @verbatim
real, imag
@endverbatim
 * so that files are easier to combine.
 * @param[in,out] multiQubit object representing the set of qubits
 */
void reportState(MultiQubit multiQubit){
	FILE *state;
	char filename[100];
	long long int index;
	sprintf(filename, "state_rank_%d.csv", multiQubit.chunkId);
	state = fopen(filename, "w");
	if (multiQubit.chunkId==0) fprintf(state, "real, imag\n");

	for(index=0; index<multiQubit.numAmpsDividedByNumChunks; index++){
		fprintf(state, "%.12f, %.12f\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
	}
	fclose(state);
}

/** Report metainformation about a set of qubits: number of qubits, number of probability amplitudes.
 * @param[in,out] multiQubit object representing the set of qubits
 * @param[in] env object representing the execution environment (local, multinode etc)
 */
void reportMultiQubitParams(MultiQubit multiQubit){
	long long int numAmps = 1L << multiQubit.numQubits;
	long long int numAmpsPerRank = numAmps/multiQubit.numChunks;
	if (multiQubit.chunkId==0){
                printf("QUBITS:\n");
                printf("Number of qubits is %d.\n", multiQubit.numQubits);
                printf("Number of amps is %lld.\n", numAmps);
		printf("Number of amps per rank is %lld.\n", numAmpsPerRank);
    }
}

int getNumQubits(MultiQubit multiQubit){
    return multiQubit.numQubits;
}

int getNumAmps(MultiQubit multiQubit){
    return multiQubit.numAmpsDividedByNumChunks*multiQubit.numChunks;
}

REAL getProbEl(MultiQubit multiQubit, long long int index){
    REAL real;
    REAL imag;
    real = getRealAmpEl(multiQubit, index);
    imag = getImagAmpEl(multiQubit, index);
    return real*real + imag*imag;
}

void rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector axis){

    double mag = sqrt(pow(axis.x,2) + pow(axis.y,2) + pow(axis.z,2));
    Vector unitAxis = {axis.x/mag, axis.y/mag, axis.z/mag};

    Complex alpha, beta;
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;    
    beta.real = sin(angle/2.0)*unitAxis.y;
    beta.imag = -sin(angle/2.0)*unitAxis.x;
    compactUnitary(multiQubit, rotQubit, alpha, beta);
}

void rotateX(MultiQubit multiQubit, const int rotQubit, REAL angle){

    Vector unitAxis = {1, 0, 0};
    rotateAroundAxis(multiQubit, rotQubit, angle, unitAxis);
}

void rotateY(MultiQubit multiQubit, const int rotQubit, REAL angle){

    Vector unitAxis = {0, 1, 0};
    rotateAroundAxis(multiQubit, rotQubit, angle, unitAxis);
}

void rotateZ(MultiQubit multiQubit, const int rotQubit, REAL angle){

    Vector unitAxis = {0, 0, 1};
    rotateAroundAxis(multiQubit, rotQubit, angle, unitAxis);
}

void controlledRotateAroundAxis(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle, Vector axis){

    double mag = sqrt(pow(axis.x,2) + pow(axis.y,2) + pow(axis.z,2));
    Vector unitAxis = {axis.x/mag, axis.y/mag, axis.z/mag};

    Complex alpha, beta;
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;    
    beta.real = sin(angle/2.0)*unitAxis.y;
    beta.imag = -sin(angle/2.0)*unitAxis.x;
    controlledCompactUnitary(multiQubit, controlQubit, targetQubit, alpha, beta);
}

void controlledRotateX(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {1, 0, 0};
    controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, unitAxis);
}

void controlledRotateY(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {0, 1, 0};
    controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, unitAxis);
}

void controlledRotateZ(MultiQubit multiQubit, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {0, 0, 1};
    controlledRotateAroundAxis(multiQubit, controlQubit, targetQubit, angle, unitAxis);
}

void sigmaZ(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, SIGMA_Z);
}

void sGate(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, S_GATE);
} 

void tGate(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, T_GATE);
}

int validateMatrixIsUnitary(ComplexMatrix2 u){

    if ( fabs(u.r0c0.real*u.r0c0.real 
        + u.r0c0.imag*u.r0c0.imag
        + u.r1c0.real*u.r1c0.real
        + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    // check
    if ( fabs(u.r0c1.real*u.r0c1.real 
        + u.r0c1.imag*u.r0c1.imag
        + u.r1c1.real*u.r1c1.real
        + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;

    if ( fabs(u.r0c0.real*u.r0c1.real 
        + u.r0c0.imag*u.r0c1.imag
        + u.r1c0.real*u.r1c1.real
        + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;

    if ( fabs(u.r0c1.real*u.r0c0.imag
        - u.r0c0.real*u.r0c1.imag
        + u.r1c1.real*u.r1c0.imag
        - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;

    return 1;
}

int validateAlphaBeta(Complex alpha, Complex beta){
    if ( fabs(alpha.real*alpha.real 
        + alpha.imag*alpha.imag
        + beta.real*beta.real 
        + beta.imag*beta.imag - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateUnitVector(REAL ux, REAL uy, REAL uz){
    if ( fabs(sqrt(ux*ux + uy*uy + uz*uz) - 1) > REAL_EPS ) return 0;
    else return 1;
}

void QuESTSeedRandomDefault(){
    // init MT random number generator with three keys -- time, pid and a hash of hostname 
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process

    struct timeval  tv;
    gettimeofday(&tv, NULL);

    double time_in_mill = 
        (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ; // convert tv_sec & tv_usec to millisecond

    unsigned long int pid = getpid();
    unsigned long int msecs = (unsigned long int) time_in_mill;
    char hostName[MAXHOSTNAMELEN+1];
    gethostname(hostName, sizeof(hostName));
    unsigned long int hostNameInt = hashString(hostName);

    unsigned long int key[3];
    key[0] = msecs; key[1] = pid; key[2] = hostNameInt;
    init_by_array(key, 3); 
}

/** 
 * numSeeds <= 64
 */
void QuESTSeedRandom(unsigned long int *seedArray, int numSeeds){
    // init MT random number generator with user defined list of seeds
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process
    init_by_array(seedArray, numSeeds); 
}

unsigned long int hashString(char *str){
    unsigned long int hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;    
}

#ifdef __cplusplus
}
#endif
