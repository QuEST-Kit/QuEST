// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * Internal functions used to implement the pure backend in ../QuEST_ops_pure.h. Do not call these functions
 * directly. In general, qubits_cpu_local.c and qubits_cpu_mpi.c will implement the pure backend by choosing
 * the correct function or combination of functions to use from those included here, which are defined in QuEST_cpu.c
 */

# ifndef QuEST_CPU_INTERNAL
# define QuEST_CPU_INTERNAL

# include "../QuEST_precision.h"

void compactUnitaryLocal (MultiQubit multiQubit, const int targetQubit, Complex alpha, Complex beta);

void compactUnitaryDistributed (MultiQubit multiQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void unitaryLocal(MultiQubit multiQubit, const int targetQubit, ComplexMatrix2 u);

void unitaryDistributed (MultiQubit multiQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void controlledCompactUnitaryLocal (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        Complex alpha, Complex beta);

void controlledCompactUnitaryDistributed (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void controlledUnitaryLocal(MultiQubit multiQubit, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

void controlledUnitaryDistributed (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void multiControlledUnitaryLocal(MultiQubit multiQubit, const int targetQubit,
        long long int mask, ComplexMatrix2 u);

void multiControlledUnitaryDistributed (MultiQubit multiQubit,
        const int targetQubit,
        long long int mask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut);

void sigmaXLocal(MultiQubit multiQubit, const int targetQubit);

void sigmaXDistributed (MultiQubit multiQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

void sigmaYLocal(MultiQubit multiQubit, const int targetQubit);

void sigmaYDistributed(MultiQubit multiQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut,
        int updateUpper);

void hadamardLocal (MultiQubit multiQubit, const int targetQubit);

void hadamardDistributed (MultiQubit multiQubit, const int targetQubit,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut, int updateUpper);


void phaseGateLocal(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);

void phaseGateDistributed(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);

void controlledNotLocal(MultiQubit multiQubit, const int controlQubit, const int targetQubit);

void controlledNotDistributed (MultiQubit multiQubit, const int controlQubit, const int targetQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut);

REAL findProbabilityOfZeroLocal (MultiQubit multiQubit,
        const int measureQubit);

REAL findProbabilityOfZeroDistributed (MultiQubit multiQubit,
        const int measureQubit);

void collapseToOutcomeLocal(MultiQubit multiQubit, int measureQubit, REAL totalProbability, int outcome);

REAL collapseToOutcomeDistributedRenorm (MultiQubit multiQubit, const int measureQubit, const REAL totalProbability);

void collapseToOutcomeDistributedSetZero(MultiQubit multiQubit, const int measureQubit);




# endif
