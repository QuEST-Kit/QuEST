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

#if defined(_WIN32) && ! defined(__MINGW32__)
  #include <Windows.h>
  #include <io.h>
  #include <process.h>
#else
  #include <unistd.h>
  #include <sys/time.h>
#endif

# include <sys/types.h> 
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

/* sets dest to be the complex conjugate of src */
void setConjugateMatrix2(ComplexMatrix2* dest, ComplexMatrix2 src) {
    
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++) {
            (dest->real)[i][j] =   src.real[i][j];
            (dest->imag)[i][j] = - src.imag[i][j]; // negative for conjugate
        }
    }
}

ComplexMatrix4 getConjugateMatrix4(ComplexMatrix4 u) {
    ComplexMatrix4 c = u;
    c.r0c0.imag *= -1; c.r0c1.imag *= -1; c.r0c2.imag *= -1; c.r0c3.imag *= -1;
    c.r1c0.imag *= -1; c.r1c1.imag *= -1; c.r1c2.imag *= -1; c.r1c3.imag *= -1;
    c.r2c0.imag *= -1; c.r2c1.imag *= -1; c.r2c2.imag *= -1; c.r2c3.imag *= -1;
    c.r3c0.imag *= -1; c.r3c1.imag *= -1; c.r3c2.imag *= -1; c.r3c3.imag *= -1;
    return c;
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
    
    qreal r0c0Phase = atan2(u.imag[0][0], u.real[0][0]);
    qreal r1c1Phase = atan2(u.imag[1][1], u.real[1][1]);
    *globalPhase = (r0c0Phase + r1c1Phase)/2.0;
    
    qreal cosPhase = cos(*globalPhase);
    qreal sinPhase = sin(*globalPhase);
    alpha->real = u.real[0][0]*cosPhase + u.imag[0][0]*sinPhase;
    alpha->imag = u.imag[0][0]*cosPhase - u.real[0][0]*sinPhase;
    beta->real  = u.real[1][0]*cosPhase + u.imag[1][0]*sinPhase;
    beta->imag  = u.imag[1][0]*cosPhase - u.real[1][0]*sinPhase;
}

void shiftIndices(int* indices, int numIndices, int shift) {
    for (int j=0; j < numIndices; j++)
        indices[j] += shift;
}

void conjugateMatrixN(ComplexMatrixN u) {
    for (long long int r=0; r < u.numRows; r++)
        for (long long int c=0; c < u.numRows; c++)
            u.elems[r][c].imag *= -1;
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
    *outcomeProb = (outcome==0)? zeroProb : 1-zeroProb;

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
#if defined(_WIN32) && ! defined(__MINGW32__)
  
    unsigned long int pid = (unsigned long int) _getpid();
    unsigned long int msecs = (unsigned long int) GetTickCount64();

    key[0] = msecs; key[1] = pid;
#else
    struct timeval  tv;
    gettimeofday(&tv, NULL);

    double time_in_mill =
        (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ; // convert tv_sec & tv_usec to millisecond

    unsigned long int pid = getpid();
    unsigned long int msecs = (unsigned long int) time_in_mill;

    key[0] = msecs; key[1] = pid;
#endif 
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
    long long int numAmps = 1LL << qureg.numQubitsInStateVec;
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

void statevec_sqrtSwapGate(Qureg qureg, int qb1, int qb2) {
    
    ComplexMatrix4 u = {0};
    u.r0c0.real=1;
    u.r3c3.real=1;
    u.r1c1.real = .5; u.r1c1.imag = .5;
    u.r1c2.real = .5; u.r1c2.imag =-.5;
    u.r2c1.real = .5; u.r2c1.imag =-.5;
    u.r2c2.real = .5; u.r2c2.imag = .5;
    
    statevec_twoQubitUnitary(qureg, qb1, qb2, u);
}

void statevec_sqrtSwapGateConj(Qureg qureg, int qb1, int qb2) {
    
    ComplexMatrix4 u = {0};
    u.r0c0.real=1;
    u.r3c3.real=1;
    u.r1c1.real = .5; u.r1c1.imag =-.5;
    u.r1c2.real = .5; u.r1c2.imag = .5;
    u.r2c1.real = .5; u.r2c1.imag = .5;
    u.r2c2.real = .5; u.r2c2.imag =-.5;
    
    statevec_twoQubitUnitary(qureg, qb1, qb2, u);
}

/** applyConj=1 will apply conjugate operation, else applyConj=0 */
void statevec_multiRotatePauli(
    Qureg qureg, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle,
    int applyConj
) {
    qreal fac = 1/sqrt(2);
    Complex uRxAlpha = {.real = fac, .imag = 0}; // Rx(pi/2)* rotates Z -> Y
    Complex uRxBeta = {.real = 0, .imag = (applyConj)? fac : -fac};
    Complex uRyAlpha = {.real = fac, .imag = 0}; // Ry(-pi/2) rotates Z -> X
    Complex uRyBeta = {.real = -fac, .imag = 0};
    
    // mask may be modified to remove superfluous Identity ops
    long long int mask = getQubitBitMask(targetQubits, numTargets);
    
    // rotate basis so that exp(Z) will effect exp(Y) and exp(X)
    for (int t=0; t < numTargets; t++) {
        if (targetPaulis[t] == PAULI_I)
            mask -= 1LL << targetPaulis[t]; // remove target from mask
        if (targetPaulis[t] == PAULI_X)
            statevec_compactUnitary(qureg, targetQubits[t], uRyAlpha, uRyBeta);
        if (targetPaulis[t] == PAULI_Y)
            statevec_compactUnitary(qureg, targetQubits[t], uRxAlpha, uRxBeta);
        // (targetPaulis[t] == 3) is Z basis
    }
    
    statevec_multiRotateZ(qureg, mask, (applyConj)? -angle : angle);
    
    // undo X and Y basis rotations
    uRxBeta.imag *= -1;
    uRyBeta.real *= -1;
    for (int t=0; t < numTargets; t++) {
        if (targetPaulis[t] == PAULI_X)
            statevec_compactUnitary(qureg, targetQubits[t], uRyAlpha, uRyBeta);
        if (targetPaulis[t] == PAULI_Y)
            statevec_compactUnitary(qureg, targetQubits[t], uRxAlpha, uRxBeta);
    }
}

void applyPauliProd(Qureg workspace, int* targetQubits, enum pauliOpType* pauliCodes, int numTargets) {
    
    // produces both pauli|qureg> or pauli * qureg (as a density matrix)
    for (int i=0; i < numTargets; i++) {
        // (pauliCodes[i] == PAULI_I) applies no operation
        if (pauliCodes[i] == PAULI_X)
            statevec_pauliX(workspace, targetQubits[i]);
        if (pauliCodes[i] == PAULI_Y)
            statevec_pauliY(workspace, targetQubits[i]);
        if (pauliCodes[i] == PAULI_Z)
            statevec_pauliZ(workspace, targetQubits[i]);
    }
}

// <pauli> = <qureg|pauli|qureg> = qureg . pauli(qureg)
qreal statevec_calcExpecPauliProd(Qureg qureg, int* targetQubits, enum pauliOpType* pauliCodes, int numTargets, Qureg workspace) {
    
    statevec_cloneQureg(workspace, qureg);
    applyPauliProd(workspace, targetQubits, pauliCodes, numTargets);
    
    // compute the expected value
    qreal value;
    if (qureg.isDensityMatrix)
        value = densmatr_calcTotalProb(workspace); // Trace(ops qureg)
    else
        value = statevec_calcInnerProduct(workspace, qureg).real; // <qureg|ops|qureg>
                
    return value;
}

qreal statevec_calcExpecPauliSum(Qureg qureg, enum pauliOpType* allCodes, qreal* termCoeffs, int numSumTerms, Qureg workspace) {
    
    int numQb = qureg.numQubitsRepresented;
    int targs[numQb];
    for (int q=0; q < numQb; q++)
        targs[q] = q;
        
    qreal value = 0;
    for (int t=0; t < numSumTerms; t++)
        value += termCoeffs[t] * statevec_calcExpecPauliProd(qureg, targs, &allCodes[t*numQb], numQb, workspace);
        
    return value;
}

void statevec_applyPauliSum(Qureg inQureg, enum pauliOpType* allCodes, qreal* termCoeffs, int numSumTerms, Qureg outQureg) {
    
    int numQb = inQureg.numQubitsRepresented;
    int targs[numQb];
    for (int q=0; q < numQb; q++)
        targs[q] = q;
        
    statevec_initBlankState(outQureg);
    
    for (int t=0; t < numSumTerms; t++) {
        Complex coef = (Complex) {.real=termCoeffs[t], .imag=0};
        Complex iden = (Complex) {.real=1, .imag=0};
        Complex zero = (Complex) {.real=0, .imag=0};
        
        // outQureg += coef paulis(inQureg)
        applyPauliProd(inQureg, targs, &allCodes[t*numQb], numQb);
        setWeightedQureg(coef, inQureg, iden, outQureg, zero, outQureg); 
        
        // undero paulis(inQureg), exploiting XX=YY=ZZ=I
        applyPauliProd(inQureg, targs, &allCodes[t*numQb], numQb);
    }
}

void statevec_twoQubitUnitary(Qureg qureg, const int targetQubit1, const int targetQubit2, ComplexMatrix4 u) {
    
    long long int ctrlMask = 0;
    statevec_multiControlledTwoQubitUnitary(qureg, ctrlMask, targetQubit1, targetQubit2, u);
}

void statevec_controlledTwoQubitUnitary(Qureg qureg, const int controlQubit, const int targetQubit1, const int targetQubit2, ComplexMatrix4 u) {
    
    long long int ctrlMask = 1LL << controlQubit;
    statevec_multiControlledTwoQubitUnitary(qureg, ctrlMask, targetQubit1, targetQubit2, u);
}

void statevec_multiQubitUnitary(Qureg qureg, int* targets, const int numTargets, ComplexMatrixN u) {
    
    long long int ctrlMask = 0;
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targets, numTargets, u);
}

void statevec_controlledMultiQubitUnitary(Qureg qureg, int ctrl, int* targets, const int numTargets, ComplexMatrixN u) {
    
    long long int ctrlMask = 1LL << ctrl;
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targets, numTargets, u);
}

/* returns conj(a) * b */
Complex getConjComplexProd(Complex a, Complex b) {
    
    Complex prod;
    prod.real = a.real*b.real + a.imag*b.imag;
    prod.imag = a.real*b.imag - a.imag*b.real;
    return prod;
}

/* adds conj(a)*b to dest */
void addConjComplexProd(Complex* dest, Complex a, Complex b) {
    
    Complex prod = getConjComplexProd(a, b);
    dest->real += prod.real;
    dest->imag += prod.imag;
}

/* sets the superOp pointer to be the super-operator of the list of ops */
void populateOneQubitKrausSuperoperator(ComplexMatrix4* superOp, ComplexMatrix2* ops, int numOps)
{
    const int opDim = 2;
    
    // clear the superop
    int opLen = opDim*opDim; 
    for (int r=0; r < opLen; r++) {
        for (int c=0; c < opLen; c++) {
            (superOp->real)[r][c] = 0;
            (superOp->imag)[r][c] = 0;
        }
    }
    
    for (int n=0; n<numOps; n++) {
        ComplexMatrix2 op = ops[n];
	
        // superop += conjugate(op) (x) op, where (x) is a tensor product
        qreal aRe, aIm, bRe, bIm;
        for (int i = 0; i < opDim; i++) {
            for (int j = 0; j < opDim; j++) { 
    			for (int k = 0; k < opDim; k++) { 
                    for (int l = 0; l < opDim; l++) { 
    					aRe = op.real[i][j];
    					aIm = -op.imag[i][j]; // minus -- conjugate of A
    					bRe = op.real[k][l];
    					bIm = op.imag[k][l];
                        
                        // @TODO: this won't compile until ComplexMatrix4 is refactored
                        (superOp->real)[i*opDim + k][j*opDim + l] +=  aRe*bRe - aIm*bIm;
    					(superOp->imag)[i*opDim + k][j*opDim + l] +=  aIm*bRe + aRe*bIm;
                    } 
                } 
            } 
        } 
    }
}

/* sets the dest pointer to be the super-operator of the list of ops */
void populateTwoQubitKrausSuperoperator(ComplexMatrixN* superOp, ComplexMatrix4* ops, int numOps) {
    
    // clear the superop
    int opLen = 16; 
    for (int r=0; r < opLen; r++)
        for (int c=0; c < opLen; c++)
            (superOp->elems)[r][c] = (Complex) {.real=0, .imag=0};
    
    for (int n=0; n < numOps; n++) {
        
        // unpack the Kraus map for convenience
        ComplexMatrix4 op = ops[n];
        Complex opArr[4][4] = {
            {op.r0c0, op.r0c1, op.r0c2, op.r0c3},
            {op.r1c0, op.r1c1, op.r1c2, op.r1c3},
            {op.r2c0, op.r2c1, op.r2c2, op.r2c3},
            {op.r3c0, op.r3c1, op.r3c2, op.r3c3}
        };
        
        // add conj(op) (x) op to the superoperator
        int opDim = 4;
        for (int i=0; i < opDim; i++)
            for (int j=0; j < opDim; j++)
                for (int k=0; k < opDim; k++)
                    for (int l=0; l < opDim; l++)
                        addConjComplexProd(&((superOp->elems)[i*opDim + k][j*opDim + l]), opArr[i][j], opArr[k][l]);
    }
}

void densmatr_applyKrausSuperoperator(Qureg qureg, int target, ComplexMatrix4 s) {
        
    long long int ctrlMask = 0;
    statevec_multiControlledTwoQubitUnitary(qureg, ctrlMask, target, target + qureg.numQubitsRepresented, s);
}

void densmatr_applyTwoQubitKrausSuperoperator(Qureg qureg, int target1, int target2, ComplexMatrixN s) {

    long long int ctrlMask = 0;
    int numQb = qureg.numQubitsRepresented;
    int targets[4] = {target1, target2, target1+numQb, target2+numQb};
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, targets, 4, s);
}

void densmatr_applyKrausMap(Qureg qureg, int target, ComplexMatrix2 *ops, int numOps) {
        
    ComplexMatrix4 superOp; 
    populateOneQubitKrausSuperoperator(&superOp, ops, numOps);
    densmatr_applyKrausSuperoperator(qureg, target, superOp);
}

void densmatr_applyTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps) {
        
    // create non-dynamic ComplexMatrixN instance
    ComplexMatrixN superOp;
    superOp.numQubits = 4;
    superOp.numRows = 1 << superOp.numQubits;
    
    // keep superOp.elems in stack
    Complex* matr[16];
    Complex flat[16][16];
    for (int r=0; r < superOp.numRows; r++)
        matr[r] = flat[r];
    superOp.elems = matr;

    populateTwoQubitKrausSuperoperator(&superOp, ops, numOps);
    densmatr_applyTwoQubitKrausSuperoperator(qureg, target1, target2, superOp);
}

void densmatr_oneQubitPauliError(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {
    
    // convert pauli probabilities into Kraus map
    ComplexMatrix2 ops[4];
    for (int n=0; n < 4; n++)
        ops[n] = (ComplexMatrix2) {.real={{0}}, .imag={{0}}};
    
    qreal facs[4] = {
		sqrt(1-(probX + probY + probZ)),
		sqrt(probX),
		sqrt(probY),
		sqrt(probZ)
	};
    ops[0].real[0][0] =  facs[0]; ops[0].real[1][1] =  facs[0];
    ops[1].real[0][1] =  facs[1]; ops[1].real[1][0] =  facs[1];
    ops[2].imag[0][1] = -facs[2]; ops[2].imag[1][0] =  facs[2];
    ops[3].real[0][0] =  facs[3]; ops[3].real[1][1] = -facs[3];
    
    densmatr_applyKrausMap(qureg, qubit, ops, 4);
}

#ifdef __cplusplus
}
#endif
