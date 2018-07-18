// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * An implementation of the pure backend in ../QuEST_ops_pure.h for a local (non-MPI, non-GPU) environment.
 * Mostly pure-state wrappers for the local/distributed functions implemented in QuEST_cpu
 */

# include "../QuEST.h"
# include "../QuEST_internal.h"
# include "../QuEST_precision.h"
# include "../mt19937ar.h"

# include "QuEST_cpu_internal.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <sys/types.h>

# ifdef _OPENMP
# include <omp.h>
# endif

void initQuESTEnv(QuESTEnv *env){
    // init MPI environment
    env->rank=0;
    env->numRanks=1;
	
	seedQuESTDefault();
}

void syncQuESTEnv(QuESTEnv env){
    // MPI Barrier goes here in MPI version. 
} 

int syncQuESTSuccess(int successCode){
    return successCode;
}

void closeQuESTEnv(QuESTEnv env){
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
    printf("Precision: size of REAL is %ld bytes\n", sizeof(REAL));
}

void reportNodeList(QuESTEnv env){
    printf("Hostname unknown: running locally\n");
}

REAL pure_calcTotalProbability(MultiQubit multiQubit){
    // implemented using Kahan summation for greater accuracy at a slight floating
    // point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    REAL pTotal=0; 
    REAL y, t, c;
    long long int index;
    long long int numAmpsPerRank = multiQubit.numAmpsPerChunk;
    c = 0.0;
    for (index=0; index<numAmpsPerRank; index++){ 
        // Perform pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index]; by Kahan

        y = multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;

        // Perform pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index]; by Kahan

        y = multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index] - c;
        t = pTotal + y;
        // Don't change the bracketing on the following line
        c = ( t - pTotal ) - y;
        pTotal = t;


    } 
    return pTotal;
}

REAL pure_getRealAmpEl(MultiQubit multiQubit, long long int index){
    return multiQubit.stateVec.real[index];
}

REAL pure_getImagAmpEl(MultiQubit multiQubit, long long int index){
    return multiQubit.stateVec.imag[index];
}

void pure_compactUnitary(MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);

    // all values required to update state vector lie in this rank
    compactUnitaryLocal(multiQubit, targetQubit, alpha, beta);
}

void pure_unitary(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    // all values required to update state vector lie in this rank
    unitaryLocal(multiQubit, targetQubit, u);
}

void pure_controlledCompactUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, Complex alpha, Complex beta) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateAlphaBeta(alpha, beta), 6, __func__);


    controlledCompactUnitaryLocal(multiQubit, controlQubit, targetQubit, alpha, beta);
}

void pure_controlledUnitary(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    controlledUnitaryLocal(multiQubit, controlQubit, targetQubit, u);
}

void pure_multiControlledUnitary(MultiQubit multiQubit, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(numControlQubits > 0 && numControlQubits <= multiQubit.numQubits, 4, __func__);
    QuESTAssert(validateMatrixIsUnitary(u), 5, __func__);

    long long int mask=0; 
    for (int i=0; i<numControlQubits; i++) mask = mask | (1LL<<controlQubits[i]);
    QuESTAssert(mask >=0 && mask <= (1LL<<multiQubit.numQubits)-1, 2, __func__);
    QuESTAssert((mask & (1LL<<targetQubit)) != (1LL<<targetQubit), 3, __func__);

    multiControlledUnitaryLocal(multiQubit, targetQubit, mask, u);
}

void pure_sigmaX(MultiQubit multiQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    sigmaXLocal(multiQubit, targetQubit);
}

void pure_sigmaY(MultiQubit multiQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    sigmaYLocal(multiQubit, targetQubit);
}

void pure_phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type)
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    phaseGateLocal(multiQubit, targetQubit, type);
}

void pure_hadamard(MultiQubit multiQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    hadamardLocal(multiQubit, targetQubit);
}

void pure_controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit) 
{
    QuESTAssert(targetQubit >= 0 && targetQubit < multiQubit.numQubits, 1, __func__);
    QuESTAssert(controlQubit >= 0 && controlQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert(controlQubit != targetQubit, 3, __func__);
    controlledNotLocal(multiQubit, controlQubit, targetQubit);
}

REAL pure_findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    REAL stateProb=0;
    stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
    if (outcome==1) stateProb = 1.0 - stateProb;
    return stateProb;
}

REAL pure_collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome)
{
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    QuESTAssert((outcome==0 || outcome==1), 10, __func__);
    REAL stateProb;
    stateProb = pure_findProbabilityOfOutcome(multiQubit, measureQubit, outcome);
    QuESTAssert(absReal(stateProb)>REAL_EPS, 8, __func__);
    collapseToOutcomeLocal(multiQubit, measureQubit, stateProb, outcome);
    return stateProb;
}

int pure_measureWithStats(MultiQubit multiQubit, int measureQubit, REAL *stateProb){
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);

    int outcome;
    // find probability of qubit being in state 1
    REAL stateProbInternal = pure_findProbabilityOfOutcome(multiQubit, measureQubit, 1);

    // we can't collapse to a state that has a probability too close to zero
    if (stateProbInternal<REAL_EPS) outcome=0;
    else if (1-stateProbInternal<REAL_EPS) outcome=1;
    else {
        // ok. both P(0) and P(1) are large enough to resolve
        // generate random float on [0,1]
        float randNum = genrand_real1();
        if (randNum<=stateProbInternal) outcome = 1;
        else outcome = 0;
    } 
    if (outcome==0) stateProbInternal = 1-stateProbInternal;
    collapseToOutcomeLocal(multiQubit, measureQubit, stateProbInternal, outcome);
    *stateProb = stateProbInternal;
    return outcome;
}

int pure_measure(MultiQubit multiQubit, int measureQubit){
    QuESTAssert(measureQubit >= 0 && measureQubit < multiQubit.numQubits, 2, __func__);
    REAL stateProb;
    return pure_measureWithStats(multiQubit, measureQubit, &stateProb);
}
