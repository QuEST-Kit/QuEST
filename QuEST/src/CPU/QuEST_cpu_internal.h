// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Internal functions used to implement the pure backend in ../QuEST_ops_pure.h. Do not call these functions
 * directly. In general, qubits_cpu_local.c and qubits_cpu_mpi.c will implement the pure backend by choosing
 * the correct function or combination of functions to use from those included here, which are defined in QuEST_cpu.c
 *
 * @author Ania Brown
 * @author Tyson Jones
 */

# ifndef QUEST_CPU_INTERNAL_H
# define QUEST_CPU_INTERNAL_H

# include "QuEST_precision.h"


/*
* Bit twiddling functions are defined seperately here in the CPU backend, 
* since the GPU backend  needs a device-specific redefinition to be callable 
* from GPU kernels. These are called in both QuEST_cpu and QuEST_cpu_distributed 
* and defined in here since public inline methods in C must go in the header
*/

inline int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber) {
    return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

inline long long int flipBit(long long int number, int bitInd) {
    return (number ^ (1LL << bitInd));
}

inline int maskContainsBit(long long int mask, int bitInd) {
    return mask & (1LL << bitInd);
}

inline int isOddParity(long long int number, int qb1, int qb2) {
    return extractBit(qb1, number) != extractBit(qb2, number);
}

inline long long int insertZeroBit(long long int number, int index) {
    long long int left, right;
    left = (number >> index) << index;
    right = number - left;
    return (left << 1) ^ right;
}

inline long long int insertTwoZeroBits(long long int number, int bit1, int bit2) {
    int small = (bit1 < bit2)? bit1 : bit2;
    int big = (bit1 < bit2)? bit2 : bit1;
    return insertZeroBit(insertZeroBit(number, small), big);
}


/*
 * density matrix operations
 */

qreal densmatr_calcPurityLocal(Qureg qureg);

void densmatr_initPureStateLocal(Qureg targetQureg, Qureg copyQureg);

qreal densmatr_calcFidelityLocal(Qureg qureg, Qureg pureState);

qreal densmatr_calcHilbertSchmidtDistanceSquaredLocal(Qureg a, Qureg b);

qreal densmatr_calcInnerProductLocal(Qureg a, Qureg b);

qreal densmatr_findProbabilityOfZeroLocal(Qureg qureg, const int measureQubit);

void densmatr_mixDepolarisingLocal(Qureg qureg, const int targetQubit, qreal depolLevel);

void densmatr_mixDepolarisingDistributed(Qureg qureg, const int targetQubit, qreal depolLevel);

void densmatr_mixDampingLocal(Qureg qureg, const int targetQubit, qreal damping);

void densmatr_mixDampingDistributed(Qureg qureg, const int targetQubit, qreal damping);

void densmatr_mixTwoQubitDepolarisingLocal(Qureg qureg, int qubit1, int qubit2, qreal delta, qreal gamma);

void densmatr_mixTwoQubitDepolarisingLocalPart1(Qureg qureg, int qubit1, int qubit2, qreal delta);

void densmatr_mixTwoQubitDepolarisingDistributed(Qureg qureg, const int targetQubit,
                const int qubit2, qreal delta, qreal gamma);

void densmatr_mixTwoQubitDepolarisingQ1LocalQ2DistributedPart3(Qureg qureg, const int targetQubit,
                const int qubit2, qreal delta, qreal gamma);


/*
 * state vector operations
 */

Complex statevec_calcInnerProductLocal(Qureg bra, Qureg ket);

void statevec_compactUnitaryLocal (Qureg qureg, const int targetQubit, Complex alpha, Complex beta);

void statevec_compactUnitaryDistributed (Qureg qureg,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_unitaryLocal(Qureg qureg, const int targetQubit, ComplexMatrix2 u);

void statevec_unitaryDistributed (Qureg qureg,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_controlledCompactUnitaryLocal (Qureg qureg, const int controlQubit, const int targetQubit,
        Complex alpha, Complex beta);

void statevec_controlledCompactUnitaryDistributed (Qureg qureg, const int controlQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_controlledUnitaryLocal(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void statevec_controlledUnitaryDistributed (Qureg qureg, const int controlQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_multiControlledUnitaryLocal(Qureg qureg, const int targetQubit,
        long long int ctrlQubitsMask, long long int ctrlFlipMask, ComplexMatrix2 u);

void statevec_multiControlledUnitaryDistributed (Qureg qureg,
        const int targetQubit,
        long long int ctrlQubitsMask, long long int ctrlFlipMask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_pauliXLocal(Qureg qureg, const int targetQubit);

void statevec_pauliXDistributed (Qureg qureg,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

void statevec_pauliYLocal(Qureg qureg, const int targetQubit, const int conjFac);

void statevec_pauliYDistributed(Qureg qureg,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut,
        int updateUpper, const int conjFac);

void statevec_controlledPauliYLocal(Qureg qureg, const int controlQubit, const int targetQubit, const int conjFactor);

void statevec_controlledPauliYDistributed(Qureg qureg, const int controlQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, const int conjFactor);
        
void statevec_hadamardLocal (Qureg qureg, const int targetQubit);

void statevec_hadamardDistributed (Qureg qureg,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut, int updateUpper);

void statevec_controlledNotLocal(Qureg qureg, const int controlQubit, const int targetQubit);

void statevec_controlledNotDistributed (Qureg qureg, const int controlQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

qreal statevec_findProbabilityOfZeroLocal (Qureg qureg, const int measureQubit);

qreal statevec_findProbabilityOfZeroDistributed (Qureg qureg);

void statevec_collapseToKnownProbOutcomeLocal(Qureg qureg, int measureQubit, int outcome, qreal totalProbability);

void statevec_collapseToKnownProbOutcomeDistributedRenorm (Qureg qureg, const int measureQubit, const qreal totalProbability);

void statevec_collapseToOutcomeDistributedSetZero(Qureg qureg);

void statevec_swapQubitAmpsLocal(Qureg qureg, int qb1, int qb2);

void statevec_swapQubitAmpsDistributed(Qureg qureg, int pairRank, int qb1, int qb2);

void statevec_multiControlledTwoQubitUnitaryLocal(Qureg qureg, long long int ctrlMask, const int q1, const int q2, ComplexMatrix4 u);

void statevec_multiControlledMultiQubitUnitaryLocal(Qureg qureg, long long int ctrlMask, int* targs, const int numTargs, ComplexMatrixN u);


# endif // QUEST_CPU_INTERNAL_H
