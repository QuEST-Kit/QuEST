// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Internal functions used to implement the pure backend in ../QuEST_ops_pure.h. Do not call these functions
 * directly. In general, qubits_cpu_local.c and qubits_cpu_mpi.c will implement the pure backend by choosing
 * the correct function or combination of functions to use from those included here, which are defined in QuEST_cpu.c
 */

# ifndef QuEST_CPU_INTERNAL
# define QuEST_CPU_INTERNAL

# include "../QuEST_precision.h"

// @TODO
void mixed_initPureStateLocal(QubitRegister targetQureg, QubitRegister copyQureg);

// @TODO
void mixed_initPureStateDistributed(QubitRegister targetQureg, QubitRegister copyQureg);


void pure_compactUnitaryLocal (QubitRegister qureg, const int targetQubit, Complex alpha, Complex beta);

void pure_compactUnitaryDistributed (QubitRegister qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void pure_unitaryLocal(QubitRegister qureg, const int targetQubit, ComplexMatrix2 u);

void pure_unitaryDistributed (QubitRegister qureg, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void pure_controlledCompactUnitaryLocal (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex alpha, Complex beta);

void pure_controlledCompactUnitaryDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void pure_controlledUnitaryLocal(QubitRegister qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void pure_controlledUnitaryDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void pure_multiControlledUnitaryLocal(QubitRegister qureg, const int targetQubit,
        long long int mask, ComplexMatrix2 u);

void pure_multiControlledUnitaryDistributed (QubitRegister qureg,
        const int targetQubit,
        long long int mask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void pure_sigmaXLocal(QubitRegister qureg, const int targetQubit);

void pure_sigmaXDistributed (QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

void pure_sigmaYLocal(QubitRegister qureg, const int targetQubit);

void pure_sigmaYDistributed(QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut,
        int updateUpper);

void pure_hadamardLocal (QubitRegister qureg, const int targetQubit);

void pure_hadamardDistributed (QubitRegister qureg, const int targetQubit,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut, int updateUpper);


void pure_specialPhaseGateLocal(QubitRegister qureg, const int targetQubit, enum phaseGateType type);

void pure_specialPhaseGateDistributed(QubitRegister qureg, const int targetQubit, enum phaseGateType type);

void pure_controlledNotLocal(QubitRegister qureg, const int controlQubit, const int targetQubit);

void pure_controlledNotDistributed (QubitRegister qureg, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

REAL pure_findProbabilityOfZeroLocal (QubitRegister qureg,
        const int measureQubit);

REAL pure_findProbabilityOfZeroDistributed (QubitRegister qureg,
        const int measureQubit);

void pure_collapseToOutcomeLocal(QubitRegister qureg, int measureQubit, REAL totalProbability, int outcome);

REAL pure_collapseToOutcomeDistributedRenorm (QubitRegister qureg, const int measureQubit, const REAL totalProbability);

void pure_collapseToOutcomeDistributedSetZero(QubitRegister qureg, const int measureQubit);




# endif
