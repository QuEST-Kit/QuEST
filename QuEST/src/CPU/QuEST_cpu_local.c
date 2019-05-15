// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the pure backend in ../QuEST_ops_pure.h for a local (non-MPI, non-GPU) environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 */

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "mt19937ar.h"

# include "QuEST_cpu_internal.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <sys/types.h>

# ifdef _OPENMP
# include <omp.h>
# endif


void densmatr_oneQubitDepolarise(Qureg qureg, const int targetQubit, qreal depolLevel) {
    if (depolLevel == 0)
        return;

    densmatr_oneQubitDepolariseLocal(qureg, targetQubit, depolLevel);
}

void densmatr_oneQubitDamping(Qureg qureg, const int targetQubit, qreal damping) {
    if (damping == 0)
        return;

    densmatr_oneQubitDampingLocal(qureg, targetQubit, damping);
}

void densmatr_twoQubitDepolarise(Qureg qureg, int qubit1, int qubit2, qreal depolLevel){
    if (depolLevel == 0)
        return;
    qreal eta = 2/depolLevel;
    qreal delta = eta - 1 - sqrt( (eta-1)*(eta-1) - 1 );
    qreal gamma = 1+delta;
    // TODO -- test delta too small

    gamma = 1/(gamma*gamma*gamma);
    densmatr_twoQubitDephase(qureg, qubit1, qubit2, depolLevel);
    densmatr_twoQubitDepolariseLocal(qureg, qubit1, qubit2, delta, gamma);
}


qreal densmatr_calcPurity(Qureg qureg) {
    return densmatr_calcPurityLocal(qureg);
}

qreal densmatr_calcFidelity(Qureg qureg, Qureg pureState) {
    
    // save pointers to qureg's pair state
    qreal* quregPairRePtr = qureg.pairStateVec.real;
    qreal* quregPairImPtr = qureg.pairStateVec.imag;
    
    // populate qureg pair state with pure state (by repointing)
    qureg.pairStateVec.real = pureState.stateVec.real;
    qureg.pairStateVec.imag = pureState.stateVec.imag;
    
    // calculate fidelity using pairState
    qreal fid = densmatr_calcFidelityLocal(qureg, pureState);
    
    // restore pointers
    qureg.pairStateVec.real = quregPairRePtr;
    qureg.pairStateVec.imag = quregPairImPtr;
    
    return fid;
}

void densmatr_initPureState(Qureg qureg, Qureg pureState) {
    
    // save pointers to qureg's pair state
    qreal* quregPairRePtr = qureg.pairStateVec.real;
    qreal* quregPairImPtr = qureg.pairStateVec.imag;
    
    // populate qureg pair state with pure state (by repointing)
    qureg.pairStateVec.real = pureState.stateVec.real;
    qureg.pairStateVec.imag = pureState.stateVec.imag;

    // populate density matrix via it's pairState
    densmatr_initPureStateLocal(qureg, pureState);
    
    // restore pointers
    qureg.pairStateVec.real = quregPairRePtr;
    qureg.pairStateVec.imag = quregPairImPtr;
}

Complex statevec_calcInnerProduct(Qureg bra, Qureg ket) {
    return statevec_calcInnerProductLocal(bra, ket);
}

qreal densmatr_calcTotalProb(Qureg qureg) {
    
    // computes the trace using Kahan summation
    qreal pTotal=0;
    qreal y, t, c;
    c = 0;
    
    long long int numCols = 1LL << qureg.numQubitsRepresented;
    long long diagIndex;
    
    for (int col=0; col< numCols; col++) {
        diagIndex = col*(numCols + 1);
        y = qureg.stateVec.real[diagIndex] - c;
        t = pTotal + y;
        c = ( t - pTotal ) - y; // brackets are important
        pTotal = t;
    }
    
    // @TODO should maybe do a cheap test that imaginary components are ~0
    
    return pTotal;
}

qreal statevec_calcTotalProb(Qureg qureg){
    // implemented using Kahan summation for greater accuracy at a slight floating
    // point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    qreal pTotal=0; 
    qreal y, t, c;
    long long int index;
    long long int numAmpsPerRank = qureg.numAmpsPerChunk;
    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){ 
        // Perform pTotal+=qureg.stateVec.real[index]*qureg.stateVec.real[index]; by Kahan

        y = qureg.stateVec.real[index]*qureg.stateVec.real[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;

        // Perform pTotal+=qureg.stateVec.imag[index]*qureg.stateVec.imag[index]; by Kahan

        y = qureg.stateVec.imag[index]*qureg.stateVec.imag[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;


    } 
    return pTotal;
}


QuESTEnv createQuESTEnv(void) {
    // init MPI environment
    
    QuESTEnv env;
    env.rank=0;
    env.numRanks=1;
    
    seedQuESTDefault();
    
    return env;
}

void syncQuESTEnv(QuESTEnv env){
    // MPI Barrier goes here in MPI version. 
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void destroyQuESTEnv(QuESTEnv env){
    // MPI finalize goes here in MPI version. Call this function anyway for consistency
}

void reportQuESTEnv(QuESTEnv env){
    printf("EXECUTION ENVIRONMENT:\n");
    printf("Running locally on one node\n");
    printf("Number of ranks is %d\n", env.numRanks);
# ifdef _OPENMP
    printf("OpenMP enabled\n");
    printf("Number of threads available is %d\n", omp_get_max_threads());
# else
    printf("OpenMP disabled\n");
# endif
    printf("Precision: size of qreal is %ld bytes\n", sizeof(qreal));
}

qreal statevec_getRealAmp(Qureg qureg, long long int index){
    return qureg.stateVec.real[index];
}

qreal statevec_getImagAmp(Qureg qureg, long long int index){
    return qureg.stateVec.imag[index];
}

void statevec_compactUnitary(Qureg qureg, const int targetQubit, Complex alpha, Complex beta) 
{
    statevec_compactUnitaryLocal(qureg, targetQubit, alpha, beta);
}

void statevec_unitary(Qureg qureg, const int targetQubit, ComplexMatrix2 u) 
{
    statevec_unitaryLocal(qureg, targetQubit, u);
}

void statevec_controlledCompactUnitary(Qureg qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    statevec_controlledCompactUnitaryLocal(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledUnitary(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u) 
{
    statevec_controlledUnitaryLocal(qureg, controlQubit, targetQubit, u);
}

void statevec_multiControlledUnitary(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) 
{
    long long int mask=0; 
    for (int i=0; i<numControlQubits; i++)
        mask = mask | (1LL<<controlQubits[i]);

    statevec_multiControlledUnitaryLocal(qureg, targetQubit, mask, u);
}

void statevec_pauliX(Qureg qureg, const int targetQubit) 
{
    statevec_pauliXLocal(qureg, targetQubit);
}

void statevec_pauliY(Qureg qureg, const int targetQubit) 
{
    int conjFac = 1;
    statevec_pauliYLocal(qureg, targetQubit, conjFac);
}

void statevec_pauliYConj(Qureg qureg, const int targetQubit) 
{
    int conjFac = -1;
    statevec_pauliYLocal(qureg, targetQubit, conjFac);
}

void statevec_controlledPauliY(Qureg qureg, const int controlQubit, const int targetQubit)
{
    int conjFac = 1;
    statevec_controlledPauliYLocal(qureg, controlQubit, targetQubit, conjFac);
}

void statevec_controlledPauliYConj(Qureg qureg, const int controlQubit, const int targetQubit)
{
    int conjFac = -1;
    statevec_controlledPauliYLocal(qureg, controlQubit, targetQubit, conjFac);
}

void statevec_hadamard(Qureg qureg, const int targetQubit) 
{
    statevec_hadamardLocal(qureg, targetQubit);
}

void statevec_controlledNot(Qureg qureg, const int controlQubit, const int targetQubit) 
{
    statevec_controlledNotLocal(qureg, controlQubit, targetQubit);
}

qreal statevec_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome)
{
    qreal stateProb=0;
    stateProb = statevec_findProbabilityOfZeroLocal(qureg, measureQubit);
    if (outcome==1) stateProb = 1.0 - stateProb;
    return stateProb;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome) {
    
    qreal outcomeProb = densmatr_findProbabilityOfZeroLocal(qureg, measureQubit);
    if (outcome == 1)
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

void statevec_collapseToKnownProbOutcome(Qureg qureg, const int measureQubit, int outcome, qreal stateProb)
{
    statevec_collapseToKnownProbOutcomeLocal(qureg, measureQubit, outcome, stateProb);
}

void seedQuESTDefault(){
    // init MT random number generator with three keys -- time and pid
    // for the MPI version, it is ok that all procs will get the same seed as random numbers will only be 
    // used by the master process

    unsigned long int key[2];
    getQuESTDefaultSeedKey(key);
    init_by_array(key, 2);
}
