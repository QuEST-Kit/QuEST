# ifndef QUBITS_INTERNAL
# define QUBITS_INTERNAL

# include "precision.h"

/** @file
 * Internal functions used to implement the public facing API in qubits.h. Do not call these functions
 * directly. In general, qubits_env_local.c and qubits_env_mpi.c will implement the public API by choosing
 * the correct function or combination of functions to use from those included here.  
 */

// functions in qubits.c

extern const char* errorCodes[];

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

// Validation

int validateMatrixIsUnitary(ComplexMatrix2 u);

int validateAlphaBeta(Complex alpha, Complex beta);

int validateUnitVector(REAL ux, REAL uy, REAL uz);

// Helper functions in qubits_env_local.c and qubits_env_mpi.c that aren't part of the public API

/**
 * Rotate a single qubit by {{1,0},{0,p}} where p is a phase term determined by the type argument
 * @param[in,out] multiQubit object representing the set of qubits
 * @param[in] targetQubit qubit to rotate
 * @param[in] type the type of phase gate to apply -- one of {SIGMA_Z, S_GATE, T_GATE}
 * */
void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);

void exitWithError(int errorCode, const char *func);

void QuESTAssert(int isValid, int errorCode, const char *func);

# endif
