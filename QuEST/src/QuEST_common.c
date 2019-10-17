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
 *
 * @author Tyson Jones
 * @author Ania Brown (seeding, reporting)
 * @author Balint Koczor (Kraus maps, mixPauli, Windows compatibility)
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

#define macro_setConjugateMatrix(dest, src, dim) { \
    for (int i=0; i<dim; i++) \
        for (int j=0; j<dim; j++) { \
            dest.real[i][j] =   src.real[i][j]; \
            dest.imag[i][j] = - src.imag[i][j]; /* negative for conjugate */ \
        } \
}
ComplexMatrix2 getConjugateMatrix2(ComplexMatrix2 src) {
    ComplexMatrix2 conj;
    macro_setConjugateMatrix(conj, src, 2);
    return conj;
}
ComplexMatrix4 getConjugateMatrix4(ComplexMatrix4 src) {
    ComplexMatrix4 conj;
    macro_setConjugateMatrix(conj, src, 4);
    return conj;
}
void setConjugateMatrixN(ComplexMatrixN m) {
    int len = 1 << m.numQubits;
    macro_setConjugateMatrix(m, m, len);
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
    
    ComplexMatrix4 u = (ComplexMatrix4) {.real={{0}}, .imag={{0}}};
    u.real[0][0]=1;
    u.real[3][3]=1;
    u.real[1][1] = .5; u.imag[1][1] = .5;
    u.real[1][2] = .5; u.imag[1][2] =-.5;
    u.real[2][1] = .5; u.imag[2][1] =-.5;
    u.real[2][2] = .5; u.imag[2][2] = .5;
    
    statevec_twoQubitUnitary(qureg, qb1, qb2, u);
}

void statevec_sqrtSwapGateConj(Qureg qureg, int qb1, int qb2) {
    
    ComplexMatrix4 u = (ComplexMatrix4) {.real={{0}}, .imag={{0}}};
    u.real[0][0]=1;
    u.real[3][3]=1;
    u.real[1][1] = .5; u.imag[1][1] =-.5;
    u.real[1][2] = .5; u.imag[1][2] = .5;
    u.real[2][1] = .5; u.imag[2][1] = .5;
    u.real[2][2] = .5; u.imag[2][2] =-.5;
    
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

#define macro_populateKrausOperator(superOp, ops, numOps, opDim) { \
    /* clear the superop */ \
    for (int r=0; r < (opDim)*(opDim); r++) \
        for (int c=0; c < (opDim)*(opDim); c++) { \
            superOp->real[r][c] = 0; \
            superOp->imag[r][c] = 0; \
        } \
    /* add each op's contribution to the superop */ \
    for (int n=0; n<(numOps); n++) \
        /* superop += conjugate(op) (x) op, where (x) is a tensor product */ \
        for (int i = 0; i < (opDim); i++) \
            for (int j = 0; j < (opDim); j++) \
    			for (int k = 0; k < (opDim); k++) \
                    for (int l = 0; l < (opDim); l++) { \
                        superOp->real[i*(opDim) + k][j*(opDim) + l] += \
                            ops[n].real[i][j]*ops[n].real[k][l] + \
                            ops[n].imag[i][j]*ops[n].imag[k][l];  \
    					superOp->imag[i*(opDim) + k][j*(opDim) + l] += \
                            ops[n].real[i][j]*ops[n].imag[k][l] - \
                            ops[n].imag[i][j]*ops[n].real[k][l];  \
                    } \
}
void populateKrausSuperOperator2(ComplexMatrix4* superOp, ComplexMatrix2* ops, int numOps) {
    int opDim = 2;
    macro_populateKrausOperator(superOp, ops, numOps, opDim);
}
void populateKrausSuperOperator4(ComplexMatrixN* superOp, ComplexMatrix4* ops, int numOps) {
    int opDim = 4;
    macro_populateKrausOperator(superOp, ops, numOps, opDim);
}
void populateKrausSuperOperatorN(ComplexMatrixN* superOp, ComplexMatrixN* ops, int numOps) {
    int opDim = 1 << ops[0].numQubits;
    macro_populateKrausOperator(superOp, ops, numOps, opDim);
}

void densmatr_applyKrausSuperoperator(Qureg qureg, int target, ComplexMatrix4 superOp) {
        
    long long int ctrlMask = 0;
    statevec_multiControlledTwoQubitUnitary(qureg, ctrlMask, target, target + qureg.numQubitsRepresented, superOp);
}

void densmatr_applyTwoQubitKrausSuperoperator(Qureg qureg, int target1, int target2, ComplexMatrixN superOp) {

    long long int ctrlMask = 0;
    int numQb = qureg.numQubitsRepresented;
    int allTargets[4] = {target1, target2, target1+numQb, target2+numQb};
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, allTargets, 4, superOp);
}

void densmatr_applyMultiQubitKrausSuperoperator(Qureg qureg, int *targets, int numTargets, ComplexMatrixN superOp) {
    long long int ctrlMask = 0;
    int allTargets[2*numTargets];
    for (int t=0; t < numTargets; t++) {
        allTargets[t] = targets[t];
        allTargets[t+numTargets] = targets[t] + qureg.numQubitsRepresented;
    }
    statevec_multiControlledMultiQubitUnitary(qureg, ctrlMask, allTargets, 2*numTargets, superOp);
}

void densmatr_mixKrausMap(Qureg qureg, int target, ComplexMatrix2 *ops, int numOps) {
        
    ComplexMatrix4 superOp; 
    populateKrausSuperOperator2(&superOp, ops, numOps);
    densmatr_applyKrausSuperoperator(qureg, target, superOp);
}

ComplexMatrixN bindArraysToStackComplexMatrixN(
    int numQubits, qreal re[][1<<numQubits], qreal im[][1<<numQubits], 
    qreal** reStorage, qreal** imStorage)
{
    ComplexMatrixN m;
    m.numQubits = numQubits;
    m.real = reStorage;
    m.imag = imStorage;
    
    int len = 1<<numQubits;
    for (int i=0; i<len; i++) {
        m.real[i] = re[i];
        m.imag[i] = im[i];
    }
    return m;
}
#define macro_initialiseStackComplexMatrixN(matrix, numQubits, real, imag) { \
    /* reStorage_ and imStorage_ must not exist in calling scope */ \
    qreal* reStorage_[1<<(numQubits)]; \
    qreal* imStorage_[1<<(numQubits)]; \
    matrix = bindArraysToStackComplexMatrixN((numQubits), real, imag, reStorage_, imStorage_); \
}
#define macro_allocStackComplexMatrixN(matrix, numQubits) { \
    /* reArr_, imArr_, reStorage_, and imStorage_ must not exist in calling scope */ \
    qreal reArr_[1<<(numQubits)][1<<(numQubits)]; \
    qreal imArr_[1<<(numQubits)][1<<(numQubits)]; \
    macro_initialiseStackComplexMatrixN(matrix, (numQubits), reArr_, imArr_); \
}

void densmatr_mixTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps) {
    
    ComplexMatrixN superOp;
    macro_allocStackComplexMatrixN(superOp, 4);
    populateKrausSuperOperator4(&superOp, ops, numOps);
    densmatr_applyTwoQubitKrausSuperoperator(qureg, target1, target2, superOp);
}

void densmatr_mixMultiQubitKrausMap(Qureg qureg, int* targets, int numTargets, ComplexMatrixN* ops, int numOps) {

    ComplexMatrixN superOp;
    macro_allocStackComplexMatrixN(superOp, 2*numTargets);
    populateKrausSuperOperatorN(&superOp, ops, numOps);
    densmatr_applyMultiQubitKrausSuperoperator(qureg, targets, numTargets, superOp);
}

void densmatr_mixPauli(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ) {
    
    // convert pauli probabilities into Kraus map
    const int numOps = 4;
    ComplexMatrix2 ops[numOps];
    for (int n=0; n < numOps; n++)
        ops[n] = (ComplexMatrix2) {.real={{0}}, .imag={{0}}};
    
    qreal facs[4] = { // literal numOps=4 for valid initialisation
		sqrt(1-(probX + probY + probZ)),
		sqrt(probX),
		sqrt(probY),
		sqrt(probZ)
	};
    ops[0].real[0][0] =  facs[0]; ops[0].real[1][1] =  facs[0];
    ops[1].real[0][1] =  facs[1]; ops[1].real[1][0] =  facs[1];
    ops[2].imag[0][1] = -facs[2]; ops[2].imag[1][0] =  facs[2];
    ops[3].real[0][0] =  facs[3]; ops[3].real[1][1] = -facs[3];
    
    densmatr_mixKrausMap(qureg, qubit, ops, numOps);
}

#ifdef __cplusplus
}
#endif
