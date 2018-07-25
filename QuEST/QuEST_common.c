// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Internal and API functions which are hardware-agnostic
 */

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "QuEST_ops.h"
# include "mt19937ar.h"

# define _BSD_SOURCE
# include <unistd.h>
# include <sys/types.h> 
# include <sys/time.h>
# include <sys/param.h>

# include <stdio.h>
# include <stdlib.h>

// error codes in QuEST_internal.h
static const char* errorMessages[] = {
	[E_INVALID_NUM_QUBITS] = "Invalid number of qubits. Must create >0.",
	[E_INVALID_TARGET_QUBIT] = "Invalid target qubit. Note qubits are zero indexed.",
	[E_INVALID_CONTROL_QUBIT] = "Invalid control qubit. Note qubits are zero indexed.",
	[E_TARGET_IS_CONTROL] = "Control qubit cannot equal target qubit.",
	[E_TARGET_IN_CONTROLS] = "Control qubits cannot include target qubit.",
	[E_INVALID_NUM_CONTROLS] = "Invalid number of control qubits. Must be >0 and <numQubits.",
	[E_NON_UNITARY_MATRIX] = "Matrix is not unitary.",
	[E_NON_UNITARY_COMPACT] = "Compact matrix is not unitary.",
	[E_SYS_TOO_BIG_TO_PRINT] = "Invalid system size. Cannot print output for systems greater than 5 qubits.",
	[E_COLLAPSE_STATE_ZERO_PROB] = "Can't collapse to state with zero probability.",
	[E_INVALID_QUBIT_OUTCOME] = "Invalid measurement outcome -- must be either 0 or 1.",
	[E_CANNOT_OPEN_FILE] = "Could not open file",
	[E_SECOND_ARG_MUST_BE_STATEVEC] = "Second argument must be a state-vector.",
	[E_MISMATCHING_REGISTER_DIMENSIONS] = "Dimensions of the qubit registers don't match.",
	[E_DEFINED_ONLY_FOR_STATEVECS] = "Valid only for state-vectors.",
	[E_DEFINED_ONLY_FOR_DENSMATRS] = "Valid only for density matrices."
};

/*
const char* errorCodes[] = {
    "Success",                                              // 0
    "Invalid target qubit. Note qubits are zero indexed.",  // 1
    "Invalid control qubit. Note qubits are zero indexed.", // 2 
    "Control qubit cannot equal target qubit.",             // 3
    "Invalid number of control qubits",                     // 4
    "Invalid unitary matrix.",                              // 5
    "as" // wot, this is actually compact non-initary       // 6
    "Invalid system size. Cannot print output for systems greater than 5 qubits.", // 7
    "Can't collapse to state with zero probability.", 		// 8
    "Invalid number of qubits.", 							// 9
    "Invalid measurement outcome -- must be either 0 or 1.",// 10
    "Could not open file.",									// 11
	"Second argument must be a pure state, not a density matrix.", // 12
	"Dimensions of the qubit registers do not match.", 		// 13
	"This operation is only defined for density matrices.",	// 14
	"This operation is only defined for two pure states.",	// 15
	"An non-unitary internal operation (phaseShift) occured.", //16
};
*/

#ifdef __cplusplus
extern "C" {
#endif

void exitWithError(ErrorCode code, const char* func){
    printf("!!!\n");
    printf("QuEST Error in function %s: %s\n", func, errorMessages[code]);
    printf("!!!\n");
    printf("exiting..\n");
    exit(code);
}

void QuESTAssert(int isValid, ErrorCode code, const char* func){
    if (!isValid) exitWithError(code, func);
}

unsigned long int hashString(char *str){
    unsigned long int hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;    
}

void seedQuESTDefault(){
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
void seedQuEST(unsigned long int *seedArray, int numSeeds){
    // init MT random number generator with user defined list of seeds
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process
    init_by_array(seedArray, numSeeds); 
}

int validateUnitComplex(Complex alpha) {
	return (absReal(1 - sqrt(alpha.real*alpha.real + alpha.imag*alpha.imag)) < REAL_EPS); 
}

int validateAlphaBeta(Complex alpha, Complex beta){
    if ( absReal(alpha.real*alpha.real 
                + alpha.imag*alpha.imag
                + beta.real*beta.real 
                + beta.imag*beta.imag - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateUnitVector(REAL ux, REAL uy, REAL uz){
    if ( absReal(sqrt(ux*ux + uy*uy + uz*uz) - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateMatrixIsUnitary(ComplexMatrix2 u){

    if ( absReal(u.r0c0.real*u.r0c0.real 
                + u.r0c0.imag*u.r0c0.imag
                + u.r1c0.real*u.r1c0.real
                + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    if ( absReal(u.r0c1.real*u.r0c1.real 
                + u.r0c1.imag*u.r0c1.imag
                + u.r1c1.real*u.r1c1.real
                + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;
    if ( absReal(u.r0c0.real*u.r0c1.real 
                + u.r0c0.imag*u.r0c1.imag
                + u.r1c0.real*u.r1c1.real
                + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;
    if ( absReal(u.r0c1.real*u.r0c0.imag
                - u.r0c0.real*u.r0c1.imag
                + u.r1c1.real*u.r1c0.imag
                - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;
    return 1;
}

REAL statevec_getProbEl(QubitRegister qureg, long long int index){
    REAL real = statevec_getRealAmpEl(qureg, index);
    REAL imag = statevec_getImagAmpEl(qureg, index);
    return real*real + imag*imag;
}

int statevec_getNumQubits(QubitRegister qureg){
    return qureg.numQubits;
}

int statevec_getNumAmps(QubitRegister qureg){
    return qureg.numAmpsPerChunk*qureg.numChunks;
}

void reportState(QubitRegister qureg){
	FILE *state;
	char filename[100];
	long long int index;
	sprintf(filename, "state_rank_%d.csv", qureg.chunkId);
	state = fopen(filename, "w");
	if (qureg.chunkId==0) fprintf(state, "real, imag\n");

	for(index=0; index<qureg.numAmpsPerChunk; index++){
		# if QuEST_PREC==1 || QuEST_PREC==2
		fprintf(state, "%.12f, %.12f\n", qureg.stateVec.real[index], qureg.stateVec.imag[index]);
		# elif QuEST_PREC == 4
		fprintf(state, "%.12Lf, %.12Lf\n", qureg.stateVec.real[index], qureg.stateVec.imag[index]);
		#endif
	}
	fclose(state);
}

void reportQubitRegisterParams(QubitRegister qureg){
    long long int numAmps = 1L << qureg.numQubits;
    long long int numAmpsPerRank = numAmps/qureg.numChunks;
    if (qureg.chunkId==0){
        printf("QUBITS:\n");
        printf("Number of qubits is %d.\n", qureg.numQubits);
        printf("Number of amps is %lld.\n", numAmps);
        printf("Number of amps per rank is %lld.\n", numAmpsPerRank);
    }
}

void statevec_phaseShift(QubitRegister qureg, const int targetQubit, REAL angle) {
	Complex term; 
	term.real = cos(angle); 
	term.imag = sin(angle);
	statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_sigmaZ(QubitRegister qureg, const int targetQubit) {
	Complex term; 
	term.real = -1;
	term.imag =  0;
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_sGate(QubitRegister qureg, const int targetQubit) {
	Complex term; 
	term.real = 0;
	term.imag = 1;
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
} 

void statevec_tGate(QubitRegister qureg, const int targetQubit) {
	Complex term; 
	term.real = 1/sqrt(2);
	term.imag = 1/sqrt(2);
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_sGateConj(QubitRegister qureg, const int targetQubit) {
	Complex term; 
	term.real =  0;
	term.imag = -1;
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
} 

void statevec_tGateConj(QubitRegister qureg, const int targetQubit) {
	Complex term; 
	term.real =  1/sqrt(2);
	term.imag = -1/sqrt(2);
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_rotateX(QubitRegister qureg, const int rotQubit, REAL angle){

    Vector unitAxis = {1, 0, 0};
    statevec_rotateAroundAxis(qureg, rotQubit, angle, unitAxis);
}

void statevec_rotateY(QubitRegister qureg, const int rotQubit, REAL angle){

    Vector unitAxis = {0, 1, 0};
    statevec_rotateAroundAxis(qureg, rotQubit, angle, unitAxis);
}

void statevec_rotateZ(QubitRegister qureg, const int rotQubit, REAL angle){

    Vector unitAxis = {0, 0, 1};
    statevec_rotateAroundAxis(qureg, rotQubit, angle, unitAxis);
}

void getAlphaBetaFromRotation(REAL angle, Vector axis, Complex* alpha, Complex* beta) {
	
    REAL mag = sqrt(pow(axis.x,2) + pow(axis.y,2) + pow(axis.z,2));
    Vector unitAxis = {axis.x/mag, axis.y/mag, axis.z/mag};

    alpha->real =   cos(angle/2.0);
    alpha->imag = - sin(angle/2.0)*unitAxis.z;	
    beta->real  =   sin(angle/2.0)*unitAxis.y;
    beta->imag  = - sin(angle/2.0)*unitAxis.x;
}

void statevec_rotateAroundAxis(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis){

    Complex alpha, beta;
    getAlphaBetaFromRotation(angle, axis, &alpha, &beta);
    statevec_compactUnitary(qureg, rotQubit, alpha, beta);
}

void statevec_rotateAroundAxisConj(QubitRegister qureg, const int rotQubit, REAL angle, Vector axis){

    Complex alpha, beta;
    getAlphaBetaFromRotation(angle, axis, &alpha, &beta);
	alpha.imag *= -1; 
	beta.imag *= -1;
    statevec_compactUnitary(qureg, rotQubit, alpha, beta);
}

void statevec_controlledRotateAroundAxis(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis){

    Complex alpha, beta;
    getAlphaBetaFromRotation(angle, axis, &alpha, &beta);
    statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledRotateAroundAxisConj(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle, Vector axis){

    Complex alpha, beta;
    getAlphaBetaFromRotation(angle, axis, &alpha, &beta);
	alpha.imag *= -1; 
	beta.imag *= -1;
    statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledRotateX(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {1, 0, 0};
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, unitAxis);
}

void statevec_controlledRotateY(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {0, 1, 0};
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, unitAxis);
}

void statevec_controlledRotateZ(QubitRegister qureg, const int controlQubit, const int targetQubit, REAL angle){

    Vector unitAxis = {0, 0, 1};
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, unitAxis);
}

Complex getConjugateScalar(Complex scalar) {
	
	Complex conjScalar;
	conjScalar.real =   scalar.real;
	conjScalar.imag = - scalar.imag;
	return conjScalar;
}

ComplexMatrix2 getConjugateMatrix(ComplexMatrix2 matrix) {
	
	ComplexMatrix2 conjMatrix;
	conjMatrix.r0c0 = getConjugateScalar(matrix.r0c0);
	conjMatrix.r0c1 = getConjugateScalar(matrix.r0c1);
	conjMatrix.r1c0 = getConjugateScalar(matrix.r1c0);
	conjMatrix.r1c1 = getConjugateScalar(matrix.r1c1);
	return conjMatrix;
}

void shiftIndices(int* indices, int numIndices, int shift) {
	for (int j=0; j < numIndices; j++)
		indices[j] += shift;
}



#ifdef __cplusplus
}
#endif