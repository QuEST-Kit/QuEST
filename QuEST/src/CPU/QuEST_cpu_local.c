// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the pure backend in ../QuEST_ops_pure.h for a local (non-MPI, non-GPU) environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 *
 * @author Ania Brown
 * @author Tyson Jones
 * @author Balint Koczor
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


void densmatr_mixDepolarising(Qureg qureg, int targetQubit, qreal depolLevel) {
    if (depolLevel == 0)
        return;

    densmatr_mixDepolarisingLocal(qureg, targetQubit, depolLevel);
}

void densmatr_mixDamping(Qureg qureg, int targetQubit, qreal damping) {
    if (damping == 0)
        return;

    densmatr_mixDampingLocal(qureg, targetQubit, damping);
}

void densmatr_mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal depolLevel){
    if (depolLevel == 0)
        return;
    qreal eta = 2/depolLevel;
    qreal delta = eta - 1 - sqrt( (eta-1)*(eta-1) - 1 );
    qreal gamma = 1+delta;
    // TODO -- test delta too small

    gamma = 1/(gamma*gamma*gamma);
    densmatr_mixTwoQubitDephasing(qureg, qubit1, qubit2, depolLevel);
    densmatr_mixTwoQubitDepolarisingLocal(qureg, qubit1, qubit2, delta, gamma);
}

qreal densmatr_calcPurity(Qureg qureg) {
    
    return densmatr_calcPurityLocal(qureg);
}

qreal densmatr_calcHilbertSchmidtDistance(Qureg a, Qureg b) {
    
    qreal distSquared = densmatr_calcHilbertSchmidtDistanceSquaredLocal(a, b);
    qreal dist = sqrt(distSquared);
    return dist;
}

qreal densmatr_calcInnerProduct(Qureg a, Qureg b) {
    
    qreal scalar = densmatr_calcInnerProductLocal(a, b);
    return scalar;
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
    
    // does not check imaginary component, by design
        
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
    
    env.seeds = NULL;
    env.numSeeds = 0;
    seedQuESTDefault(&env);
    
    return env;
}

void syncQuESTEnv(QuESTEnv env){
    // MPI Barrier goes here in MPI version. 
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void destroyQuESTEnv(QuESTEnv env){
    free(env.seeds);
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

void getEnvironmentString(QuESTEnv env, char str[200]){
    int ompStatus=0;
    int numThreads=1;
# ifdef _OPENMP
    ompStatus=1;
    numThreads=omp_get_max_threads(); 
# endif
    sprintf(str, "CUDA=0 OpenMP=%d MPI=0 threads=%d ranks=1", ompStatus, numThreads);
}

qreal statevec_getRealAmp(Qureg qureg, long long int index){
    return qureg.stateVec.real[index];
}

qreal statevec_getImagAmp(Qureg qureg, long long int index){
    return qureg.stateVec.imag[index];
}

void statevec_compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta) 
{
    statevec_compactUnitaryLocal(qureg, targetQubit, alpha, beta);
}

void statevec_unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u) 
{
    statevec_unitaryLocal(qureg, targetQubit, u);
}

void statevec_controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta) 
{
    statevec_controlledCompactUnitaryLocal(qureg, controlQubit, targetQubit, alpha, beta);
}

void statevec_controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u) 
{
    statevec_controlledUnitaryLocal(qureg, controlQubit, targetQubit, u);
}

void statevec_multiControlledUnitary(Qureg qureg, long long int ctrlQubitsMask, long long int ctrlFlipMask, int targetQubit, ComplexMatrix2 u) 
{
    statevec_multiControlledUnitaryLocal(qureg, targetQubit, ctrlQubitsMask, ctrlFlipMask, u);
}

void statevec_pauliX(Qureg qureg, int targetQubit) 
{
    statevec_pauliXLocal(qureg, targetQubit);
}

void statevec_pauliY(Qureg qureg, int targetQubit) 
{
    int conjFac = 1;
    statevec_pauliYLocal(qureg, targetQubit, conjFac);
}

void statevec_pauliYConj(Qureg qureg, int targetQubit) 
{
    int conjFac = -1;
    statevec_pauliYLocal(qureg, targetQubit, conjFac);
}

void statevec_controlledPauliY(Qureg qureg, int controlQubit, int targetQubit)
{
    int conjFac = 1;
    statevec_controlledPauliYLocal(qureg, controlQubit, targetQubit, conjFac);
}

void statevec_controlledPauliYConj(Qureg qureg, int controlQubit, int targetQubit)
{
    int conjFac = -1;
    statevec_controlledPauliYLocal(qureg, controlQubit, targetQubit, conjFac);
}

void statevec_hadamard(Qureg qureg, int targetQubit) 
{
    statevec_hadamardLocal(qureg, targetQubit);
}

void statevec_controlledNot(Qureg qureg, int controlQubit, int targetQubit) 
{
    statevec_controlledNotLocal(qureg, controlQubit, targetQubit);
}

void statevec_multiControlledMultiQubitNot(Qureg qureg, int ctrlMask, int targMask) {
    
    statevec_multiControlledMultiQubitNotLocal(qureg, ctrlMask, targMask);
}

qreal statevec_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome)
{
    qreal outcomeProb = statevec_findProbabilityOfZeroLocal(qureg, measureQubit);
    if (outcome==1)
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

qreal densmatr_calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome) {
    
    qreal outcomeProb = densmatr_findProbabilityOfZeroLocal(qureg, measureQubit);
    if (outcome == 1)
        outcomeProb = 1.0 - outcomeProb;
    return outcomeProb;
}

void statevec_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits) {
    
    statevec_calcProbOfAllOutcomesLocal(retProbs, qureg, qubits, numQubits);
}

void densmatr_calcProbOfAllOutcomes(qreal* retProbs, Qureg qureg, int* qubits, int numQubits) {
    
    densmatr_calcProbOfAllOutcomesLocal(retProbs, qureg, qubits, numQubits);
}

void statevec_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal stateProb)
{
    statevec_collapseToKnownProbOutcomeLocal(qureg, measureQubit, outcome, stateProb);
}

void seedQuEST(QuESTEnv *env, unsigned long int *seedArray, int numSeeds) {

    // free existing seed array, if exists
    if (env->seeds != NULL)
        free(env->seeds);
        
    // record keys in permanent heap
    env->seeds = malloc(numSeeds * sizeof *(env->seeds));
    for (int i=0; i<numSeeds; i++)
        (env->seeds)[i] = seedArray[i];
    env->numSeeds = numSeeds;
    
    // pass keys to Mersenne Twister seeder
    init_by_array(seedArray, numSeeds); 
}

void statevec_multiControlledTwoQubitUnitary(Qureg qureg, long long int ctrlMask, int q1, int q2, ComplexMatrix4 u)
{
    statevec_multiControlledTwoQubitUnitaryLocal(qureg, ctrlMask, q1, q2, u);
}

void statevec_multiControlledMultiQubitUnitary(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u)
{
    statevec_multiControlledMultiQubitUnitaryLocal(qureg, ctrlMask, targs, numTargs, u);
}

void statevec_swapQubitAmps(Qureg qureg, int qb1, int qb2) 
{
    statevec_swapQubitAmpsLocal(qureg, qb1, qb2);
}

void densmatr_applyDiagonalOp(Qureg qureg, DiagonalOp op) {

    // we must preload qureg.pairStateVec with the elements of op.
    // instead of needless cloning, we'll just temporarily swap the pointers
    qreal* rePtr = qureg.pairStateVec.real;
    qreal* imPtr = qureg.pairStateVec.imag;
    qureg.pairStateVec.real = op.real;
    qureg.pairStateVec.imag = op.imag;
    
    densmatr_applyDiagonalOpLocal(qureg, op);
    
    qureg.pairStateVec.real = rePtr;
    qureg.pairStateVec.imag = imPtr;
}

Complex statevec_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {

    return statevec_calcExpecDiagonalOpLocal(qureg, op);
}

Complex densmatr_calcExpecDiagonalOp(Qureg qureg, DiagonalOp op) {
    
    return densmatr_calcExpecDiagonalOpLocal(qureg, op);
}
