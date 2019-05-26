// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Internal and API functions which are hardware-agnostic.
 * These must never call a front-end function in QuEST.c, which would lead to 
 * duplication of e.g. QASM logging and validation. Note that though many of
 * these functions are prefixed with statevec_, they will be called multiple times
 * to effect their equivalent operation on density matrices, so the passed Qureg
 * can be assumed a statevector. Functions prefixed with densmatr_ may still
 * explicitly call statevec_ functions, but will need to manually apply the
 * conjugate qubit-shifted operations to satisfy the Choi–Jamiolkowski isomorphism
 */

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "mt19937ar.h"

# include <unistd.h>
# include <sys/types.h> 
# include <sys/time.h>
# include <stdio.h>
# include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif

/* builds a bit-string where 1 indicates a qubit is present in this list */
long long int getQubitBitMask(int* qubits, const int numQubits) {
    
    long long int mask=0; 
    for (int i=0; i<numQubits; i++)
        mask = mask | (1LL << qubits[i]);
        
    return mask;
}

/* builds a bit-string where 1 indicates control qubits conditioned on 0 ('flipped') */
long long int getControlFlipMask(int* controlQubits, int* controlState, const int numControlQubits) {
    
    long long int mask=0;
    for (int i=0; i<numControlQubits; i++)
        if (controlState[i] == 0)
            mask = mask | (1LL << controlQubits[i]);
            
    return mask;
}

void ensureIndsIncrease(int* ind1, int* ind2) {
    if (*ind1 > *ind2) {
        int copy = *ind1;
        *ind1 = *ind2;
        *ind2 = copy;
    }
}

qreal getVectorMagnitude(Vector vec) {
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

Vector getUnitVector(Vector vec) {
    
    qreal mag = getVectorMagnitude(vec);
    Vector unitVec = (Vector) {.x=vec.x/mag, .y=vec.y/mag, .z=vec.z/mag};
    return unitVec;
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

void getComplexPairFromRotation(qreal angle, Vector axis, Complex* alpha, Complex* beta) {
    
    Vector unitAxis = getUnitVector(axis);
    alpha->real =   cos(angle/2.0);
    alpha->imag = - sin(angle/2.0)*unitAxis.z;  
    beta->real  =   sin(angle/2.0)*unitAxis.y;
    beta->imag  = - sin(angle/2.0)*unitAxis.x;
}

/** maps U(alpha, beta) to Rz(rz2) Ry(ry) Rz(rz1) */
void getZYZRotAnglesFromComplexPair(Complex alpha, Complex beta, qreal* rz2, qreal* ry, qreal* rz1) {
    
    qreal alphaMag = sqrt(alpha.real*alpha.real + alpha.imag*alpha.imag);
    *ry = 2.0 * acos(alphaMag);
    
    qreal alphaPhase = atan2(alpha.imag, alpha.real);
    qreal betaPhase  = atan2(beta.imag,  beta.real);
    *rz2 = - alphaPhase + betaPhase;
    *rz1 = - alphaPhase - betaPhase;
}

/** maps U(r0c0, r0c1, r1c0, r1c1) to exp(i globalPhase) U(alpha, beta) */
void getComplexPairAndPhaseFromUnitary(ComplexMatrix2 u, Complex* alpha, Complex* beta, qreal* globalPhase) {
    
    qreal r0c0Phase = atan2(u.r0c0.imag, u.r0c0.real);
    qreal r1c1Phase = atan2(u.r1c1.imag, u.r1c1.real);
    *globalPhase = (r0c0Phase + r1c1Phase)/2.0;
    
    qreal cosPhase = cos(*globalPhase);
    qreal sinPhase = sin(*globalPhase);
    alpha->real = u.r0c0.real*cosPhase + u.r0c0.imag*sinPhase;
    alpha->imag = u.r0c0.imag*cosPhase - u.r0c0.real*sinPhase;
    beta->real = u.r1c0.real*cosPhase + u.r1c0.imag*sinPhase;
    beta->imag = u.r1c0.imag*cosPhase - u.r1c0.real*sinPhase;
}

void shiftIndices(int* indices, int numIndices, int shift) {
    for (int j=0; j < numIndices; j++)
        indices[j] += shift;
}

int generateMeasurementOutcome(qreal zeroProb, qreal *outcomeProb) {
    
    // randomly choose outcome
    int outcome;
    if (zeroProb < REAL_EPS) 
        outcome = 1;
    else if (1-zeroProb < REAL_EPS) 
        outcome = 0;
    else
        outcome = (genrand_real1() > zeroProb);
    
    // set probability of outcome
    if (outcome == 0)
        *outcomeProb = zeroProb;
    else
        *outcomeProb = 1 - zeroProb;
    
    return outcome;
}

unsigned long int hashString(char *str){
    unsigned long int hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;    
}

void getQuESTDefaultSeedKey(unsigned long int *key){
    // init MT random number generator with two keys -- time and pid
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process

    struct timeval  tv;
    gettimeofday(&tv, NULL);

    double time_in_mill =
        (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ; // convert tv_sec & tv_usec to millisecond

    unsigned long int pid = getpid();
    unsigned long int msecs = (unsigned long int) time_in_mill;

    key[0] = msecs; key[1] = pid;
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

void reportState(Qureg qureg){
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

void reportQuregParams(Qureg qureg){
    long long int numAmps = 1L << qureg.numQubitsInStateVec;
    long long int numAmpsPerRank = numAmps/qureg.numChunks;
    if (qureg.chunkId==0){
        printf("QUBITS:\n");
        printf("Number of qubits is %d.\n", qureg.numQubitsInStateVec);
        printf("Number of amps is %lld.\n", numAmps);
        printf("Number of amps per rank is %lld.\n", numAmpsPerRank);
    }
}

qreal statevec_getProbAmp(Qureg qureg, long long int index){
    qreal real = statevec_getRealAmp(qureg, index);
    qreal imag = statevec_getImagAmp(qureg, index);
    return real*real + imag*imag;
}

void statevec_phaseShift(Qureg qureg, const int targetQubit, qreal angle) {
    Complex term; 
    term.real = cos(angle); 
    term.imag = sin(angle);
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_pauliZ(Qureg qureg, const int targetQubit) {
    Complex term; 
    term.real = -1;
    term.imag =  0;
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_sGate(Qureg qureg, const int targetQubit) {
    Complex term; 
    term.real = 0;
    term.imag = 1;
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
} 

void statevec_tGate(Qureg qureg, const int targetQubit) {
    Complex term; 
    term.real = 1/sqrt(2);
    term.imag = 1/sqrt(2);
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_sGateConj(Qureg qureg, const int targetQubit) {
    Complex term; 
    term.real =  0;
    term.imag = -1;
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
} 

void statevec_tGateConj(Qureg qureg, const int targetQubit) {
    Complex term; 
    term.real =  1/sqrt(2);
    term.imag = -1/sqrt(2);
    statevec_phaseShiftByTerm(qureg, targetQubit, term);
}

void statevec_rotateX(Qureg qureg, const int rotQubit, qreal angle){

    Vector unitAxis = {1, 0, 0};
    statevec_rotateAroundAxis(qureg, rotQubit, angle, unitAxis);
}

void statevec_rotateY(Qureg qureg, const int rotQubit, qreal angle){

    Vector unitAxis = {0, 1, 0};
    statevec_rotateAroundAxis(qureg, rotQubit, angle, unitAxis);
}

void statevec_rotateZ(Qureg qureg, const int rotQubit, qreal angle){

    Vector unitAxis = {0, 0, 1};
    statevec_rotateAroundAxis(qureg, rotQubit, angle, unitAxis);
}

void statevec_rotateAroundAxis(Qureg qureg, const int rotQubit, qreal angle, Vector axis){

    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    statevec_compactUnitary(qureg, rotQubit, alpha, beta);
}

void statevec_rotateAroundAxisConj(Qureg qureg, const int rotQubit, qreal angle, Vector axis){

    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    alpha.imag *= -1; 
    beta.imag *= -1;
    statevec_compactUnitary(qureg, rotQubit, alpha, beta);
}

void statevec_controlledRotateAroundAxis(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle, Vector axis){

    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledRotateAroundAxisConj(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle, Vector axis){

    Complex alpha, beta;
    getComplexPairFromRotation(angle, axis, &alpha, &beta);
    alpha.imag *= -1; 
    beta.imag *= -1;
    statevec_controlledCompactUnitary(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledRotateX(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle){

    Vector unitAxis = {1, 0, 0};
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, unitAxis);
}

void statevec_controlledRotateY(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle){

    Vector unitAxis = {0, 1, 0};
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, unitAxis);
}

void statevec_controlledRotateZ(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle){

    Vector unitAxis = {0, 0, 1};
    statevec_controlledRotateAroundAxis(qureg, controlQubit, targetQubit, angle, unitAxis);
}

int statevec_measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb) {
    
    qreal zeroProb = statevec_calcProbOfOutcome(qureg, measureQubit, 0);
    int outcome = generateMeasurementOutcome(zeroProb, outcomeProb);
    statevec_collapseToKnownProbOutcome(qureg, measureQubit, outcome, *outcomeProb);
    return outcome;
}

int densmatr_measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb) {
    
    qreal zeroProb = densmatr_calcProbOfOutcome(qureg, measureQubit, 0);
    int outcome = generateMeasurementOutcome(zeroProb, outcomeProb);
    densmatr_collapseToKnownProbOutcome(qureg, measureQubit, outcome, *outcomeProb);
    return outcome;
}

qreal statevec_calcFidelity(Qureg qureg, Qureg pureState) {
    
    Complex innerProd = statevec_calcInnerProduct(qureg, pureState);
    qreal innerProdMag = innerProd.real*innerProd.real + innerProd.imag*innerProd.imag;
    return innerProdMag;
}

void statevec_swapGate(Qureg qureg, int qb1, int qb2) {

    statevec_controlledNot(qureg, qb1, qb2);
    statevec_controlledNot(qureg, qb2, qb1);
    statevec_controlledNot(qureg, qb1, qb2);
}

void statevec_sqrtSwapGate(Qureg qureg, int qb1, int qb2) {
    
    ComplexMatrix2 u;
    u.r0c0.real = .5; u.r0c0.imag = .5;
    u.r0c1.real = .5; u.r0c1.imag =-.5;
    u.r1c0.real = .5; u.r1c0.imag =-.5;
    u.r1c1.real = .5; u.r1c1.imag = .5;
    
    statevec_controlledNot(qureg, qb1, qb2);
    statevec_controlledUnitary(qureg, qb2, qb1, u);
    statevec_controlledNot(qureg, qb1, qb2);
}

void statevec_sqrtSwapGateConj(Qureg qureg, int qb1, int qb2) {
    
    ComplexMatrix2 u;
    u.r0c0.real = .5; u.r0c0.imag =-.5;
    u.r0c1.real = .5; u.r0c1.imag = .5;
    u.r1c0.real = .5; u.r1c0.imag = .5;
    u.r1c1.real = .5; u.r1c1.imag =-.5;
    
    statevec_controlledNot(qureg, qb1, qb2);
    statevec_controlledUnitary(qureg, qb2, qb1, u);
    statevec_controlledNot(qureg, qb1, qb2);
}

void densmatr_oneQubitPauliError(Qureg qureg, int qubit, qreal pX, qreal pY, qreal pZ) {
    
    // The accepted pX, pY, pZ do NOT need to be pre-scaled.
    // The passed probabilities below are modified so that the final state produced is
    // q + px X q X + py Y q Y + pz Z q Z
    // The *2 are due to oneQubitDephase accepting 'dephaseLevel', not a probability
    // The double statevec_compactUnitary calls are needed (with complex conjugates)
    // since we're operating on a density matrix, not a statevector, and this function
    // won't be called twice from the front-end
    
    Complex alpha, beta;
    qreal fac = 1/sqrt(2);
    int numQb = qureg.numQubitsRepresented;
    
    // add Z error
    densmatr_oneQubitDephase(qureg, qubit, 2 * pZ/(1-pX-pY));
    
    // rotate basis via Rx(pi/2)
    alpha.real = fac; alpha.imag = 0;
    beta.real = 0;    beta.imag = -fac;
    statevec_compactUnitary(qureg, qubit, alpha, beta);
    alpha.imag *= -1; beta.imag *= -1;
    statevec_compactUnitary(qureg, qubit + numQb, alpha, beta);
    
    // add Z -> Y Rx(pi/2) error 
    densmatr_oneQubitDephase(qureg, qubit, 2 * pY/(1-pX));
    
    // rotate basis by Rx(-pi/2) then Ry(pi/2) 
    alpha.real = .5; alpha.imag = -.5;
    beta.real = .5;  beta.imag = .5;
    statevec_compactUnitary(qureg, qubit, alpha, beta);
    alpha.imag *= -1; beta.imag *= -1;
    statevec_compactUnitary(qureg, qubit + numQb, alpha, beta);
    
    // add Z -> X Ry(pi/2) error
    densmatr_oneQubitDephase(qureg, qubit, 2 * pX);
    
    // restore basis
    alpha.real = fac; alpha.imag = 0;
    beta.real = -fac; beta.imag = 0;
    statevec_compactUnitary(qureg, qubit, alpha, beta);
    alpha.imag *= -1; beta.imag *= -1;
    statevec_compactUnitary(qureg, qubit + numQb, alpha, beta);
}

/** applyConj=1 will apply conjugate operation, else applyConj=0 */
void statevec_multiRotatePauli(
    Qureg qureg, int* targetQubits, int* targetPaulis, int numTargets, qreal angle,
    int applyConj
) {
    qreal fac = 1/sqrt(2);
    Complex uRxAlpha = {.real = fac, .imag = 0}; // Rx(pi/2)* rotates Z -> Y
    Complex uRxBeta = {.real = 0, .imag = (applyConj)? fac : -fac};
    Complex uRyAlpha = {.real = fac, .imag = 0}; // Ry(pi/2) rotates Z -> X
    Complex uRyBeta = {.real = fac, .imag = 0};
    
    // mask may be modified to remove superfluous Identity ops
    long long int mask = getQubitBitMask(targetQubits, numTargets);
    
    // rotate basis so that exp(Z) will effect exp(Y) and exp(X)
    for (int t=0; t < numTargets; t++) {
        if (targetPaulis[t] == 0)
            mask -= 1LL << targetPaulis[t]; // remove target from mask
        if (targetPaulis[t] == 1)
            statevec_compactUnitary(qureg, targetQubits[t], uRyAlpha, uRyBeta);
        if (targetPaulis[t] == 2)
            statevec_compactUnitary(qureg, targetQubits[t], uRxAlpha, uRxBeta);
        // (targetPaulis[t] == 3) is Z basis
    }
    
    statevec_multiRotateZ(qureg, mask, (applyConj)? -angle : angle);
    
    // undo X and Y basis rotations
    uRxBeta.imag *= -1;
    uRyBeta.real *= -1;
    for (int t=0; t < numTargets; t++) {
        if (targetPaulis[t] == 1)
            statevec_compactUnitary(qureg, targetQubits[t], uRyAlpha, uRyBeta);
        if (targetPaulis[t] == 2)
            statevec_compactUnitary(qureg, targetQubits[t], uRxAlpha, uRxBeta);
    }
}

#ifdef __cplusplus
}
#endif
