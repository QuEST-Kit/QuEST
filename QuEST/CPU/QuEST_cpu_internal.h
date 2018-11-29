// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Internal functions used to implement the pure backend in ../QuEST_ops_pure.h. Do not call these functions
 * directly. In general, qubits_cpu_local.c and qubits_cpu_mpi.c will implement the pure backend by choosing
 * the correct function or combination of functions to use from those included here, which are defined in QuEST_cpu.c
 */

# ifndef QuEST_CPU_INTERNAL
# define QuEST_CPU_INTERNAL

# include "../QuEST_precision.h"

REAL densmatr_calcPurityLocal(Qureg qureg);

void densmatr_initPureStateLocal(Qureg targetQureg, Qureg copyQureg);

REAL densmatr_calcFidelityLocal(Qureg qureg, Qureg pureState);

REAL densmatr_findProbabilityOfZeroLocal(Qureg qureg, const int measureQubit);

void densmatr_oneQubitDepolariseLocal(Qureg qureg, const int targetQubit, REAL depolLevel);

void densmatr_oneQubitDepolariseDistributed(Qureg qureg, const int targetQubit, REAL depolLevel);

void densmatr_twoQubitDepolariseLocal(Qureg qureg, int qubit1, int qubit2, REAL delta, REAL gamma);

void densmatr_twoQubitDepolariseLocalPart1(Qureg qureg, int qubit1, int qubit2, REAL delta);

void densmatr_twoQubitDepolariseDistributed(Qureg qureg, const int targetQubit,
                const int qubit2, REAL delta, REAL gamma);

void densmatr_twoQubitDepolariseQ1LocalQ2DistributedPart3(Qureg qureg, const int targetQubit,
                const int qubit2, REAL delta, REAL gamma);

Complex statevec_calcInnerProductLocal(Qureg bra, Qureg ket);

void statevec_compactUnitaryLocal (Qureg qureg, const int targetQubit, Complex alpha, Complex beta);

void statevec_compactUnitaryDistributed (Qureg qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_unitaryLocal(Qureg qureg, const int targetQubit, ComplexMatrix2 u);

void statevec_unitaryDistributed (Qureg qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_controlledCompactUnitaryLocal (Qureg qureg, const int controlQubit, const int targetQubit,
        Complex alpha, Complex beta);

void statevec_controlledCompactUnitaryDistributed (Qureg qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_controlledUnitaryLocal(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void statevec_controlledUnitaryDistributed (Qureg qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_multiControlledUnitaryLocal(Qureg qureg, const int targetQubit,
        long long int mask, ComplexMatrix2 u);

void statevec_multiControlledUnitaryDistributed (Qureg qureg,
        const int targetQubit,
        long long int mask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_pauliXLocal(Qureg qureg, const int targetQubit);

void statevec_pauliXDistributed (Qureg qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

void statevec_pauliYLocal(Qureg qureg, const int targetQubit, int conjFac);

void statevec_pauliYDistributed(Qureg qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut,
        int updateUpper, int conjFac);

void statevec_controlledPauliYLocal(Qureg qureg, const int controlQubit, const int targetQubit, const int conjFactor);

void statevec_controlledPauliYDistributed(Qureg qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, const int conjFactor);
        
void statevec_hadamardLocal (Qureg qureg, const int targetQubit);

void statevec_hadamardDistributed (Qureg qureg, const int targetQubit,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut, int updateUpper);

void statevec_controlledNotLocal(Qureg qureg, const int controlQubit, const int targetQubit);

void statevec_controlledNotDistributed (Qureg qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

REAL statevec_findProbabilityOfZeroLocal (Qureg qureg, const int measureQubit);

REAL statevec_findProbabilityOfZeroDistributed (Qureg qureg, const int measureQubit);

void statevec_collapseToKnownProbOutcomeLocal(Qureg qureg, int measureQubit, int outcome, REAL totalProbability);

void statevec_collapseToKnownProbOutcomeDistributedRenorm (Qureg qureg, const int measureQubit, const REAL totalProbability);

void statevec_collapseToOutcomeDistributedSetZero(Qureg qureg);




# endif
