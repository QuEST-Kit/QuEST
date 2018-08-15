// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Internal functions used to implement the pure backend in ../QuEST_ops_pure.h. Do not call these functions
 * directly. In general, qubits_cpu_local.c and qubits_cpu_mpi.c will implement the pure backend by choosing
 * the correct function or combination of functions to use from those included here, which are defined in QuEST_cpu.c
 */

# ifndef QuEST_CPU_INTERNAL
# define QuEST_CPU_INTERNAL

# include "../QuEST_precision.h"

void densmatr_initPureStateLocal(QubitRegister targetQureg, QubitRegister copyQureg);

REAL densmatr_calcFidelityLocal(QubitRegister qureg, QubitRegister pureState);

REAL densmatr_findProbabilityOfZeroLocal(QubitRegister qureg, const int measureQubit);

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

void statevec_sigmaXLocal(QubitRegister qureg, const int targetQubit);

void statevec_sigmaXDistributed (QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

void statevec_sigmaYLocal(QubitRegister qureg, const int targetQubit, int conjFac);

void statevec_sigmaYDistributed(QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut,
        int updateUpper, int conjFac);

void statevec_controlledSigmaYLocal(QubitRegister qureg, const int controlQubit, const int targetQubit, const int conjFactor);

void statevec_controlledSigmaYDistributed(QubitRegister qureg, const int controlQubit, const int targetQubit,
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
