// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Internal functions used to implement the pure backend in ../QuEST_ops_pure.h. Do not call these functions
 * directly. In general, qubits_cpu_local.c and qubits_cpu_mpi.c will implement the pure backend by choosing
 * the correct function or combination of functions to use from those included here, which are defined in QuEST_cpu.c
 */

# ifndef QuEST_CPU_INTERNAL
# define QuEST_CPU_INTERNAL

# include "../QuEST_precision.h"

REAL densmatr_calcPurityLocal(QubitRegister qureg);

void densmatr_initPureStateLocal(QubitRegister targetQureg, QubitRegister copyQureg);

REAL densmatr_calcFidelityLocal(QubitRegister qureg, QubitRegister pureState);

REAL densmatr_findProbabilityOfZeroLocal(QubitRegister qureg, const int measureQubit);

void densmatr_oneQubitDepolariseLocal(QubitRegister qureg, const int targetQubit, REAL depolLevel);

void densmatr_oneQubitDepolariseDistributed(QubitRegister qureg, const int targetQubit, REAL depolLevel);

void densmatr_twoQubitDepolariseLocal(QubitRegister qureg, int qubit1, int qubit2, REAL delta, REAL gamma);

void densmatr_twoQubitDepolariseLocalPart1(QubitRegister qureg, int qubit1, int qubit2, REAL delta);

void densmatr_twoQubitDepolariseDistributed(QubitRegister qureg, const int targetQubit,
                const int qubit2, REAL delta, REAL gamma);

void densmatr_twoQubitDepolariseQ1LocalQ2DistributedPart3(QubitRegister qureg, const int targetQubit,
                const int qubit2, REAL delta, REAL gamma);

Complex statevec_calcInnerProductLocal(QubitRegister bra, QubitRegister ket);

void statevec_compactUnitaryLocal (QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta);

void statevec_compactUnitaryDistributed (QubitRegister qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_unitaryLocal(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u);

void statevec_unitaryDistributed (QubitRegister qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_controlledCompactUnitaryLocal (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex alpha, Complex beta);

void statevec_controlledCompactUnitaryDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_controlledUnitaryLocal(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void statevec_controlledUnitaryDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_multiControlledUnitaryLocal(QubitRegister qureg, const int targetQubit,
        long long int mask, ComplexMatrix2 u);

void statevec_multiControlledUnitaryDistributed (QubitRegister qureg,
        const int targetQubit,
        long long int mask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void statevec_pauliXLocal(QubitRegister qureg, const int targetQubit);

void statevec_pauliXDistributed (QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

void statevec_pauliYLocal(QubitRegister qureg, const int targetQubit, int conjFac);

void statevec_pauliYDistributed(QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut,
        int updateUpper, int conjFac);

void statevec_controlledPauliYLocal(QubitRegister qureg, const int controlQubit, const int targetQubit, const int conjFactor);

void statevec_controlledPauliYDistributed(QubitRegister qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, const int conjFactor);
        
void statevec_hadamardLocal (QubitRegister qureg, const int targetQubit);

void statevec_hadamardDistributed (QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut, int updateUpper);

void statevec_controlledNotLocal(QubitRegister qureg, const int controlQubit, const int targetQubit);

void statevec_controlledNotDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

REAL statevec_findProbabilityOfZeroLocal (QubitRegister qureg, const int measureQubit);

REAL statevec_findProbabilityOfZeroDistributed (QubitRegister qureg, const int measureQubit);

void statevec_collapseToKnownProbOutcomeLocal(QubitRegister qureg, int measureQubit, int outcome, REAL totalProbability);

void statevec_collapseToKnownProbOutcomeDistributedRenorm (QubitRegister qureg, const int measureQubit, const REAL totalProbability);

void statevec_collapseToOutcomeDistributedSetZero(QubitRegister qureg);




# endif
