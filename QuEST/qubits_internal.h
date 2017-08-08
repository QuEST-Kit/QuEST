# ifndef QUBITS_INTERNAL
# define QUBITS_INTERNAL

# include "precision.h"

/** @file
 * Internal functions used to implement the public facing API in qubits.h. Do not call these functions
 * directly. In general, qubits_env_local.c and qubits_env_mpi.c will implement the public API by choosing
 * the correct function or combination of functions to use from those included here.  
 */

void rotateQubitLocal (MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta);

void rotateQubitDistributed (MultiQubit multiQubit, const int rotQubit,
                Complex rot1, Complex rot2,
                ComplexArray stateVecUp,
                ComplexArray stateVecLo,
                ComplexArray stateVecOut);

void controlRotateQubitLocal (MultiQubit multiQubit, const int rotQubit, const int controlQubit,
		Complex alpha, Complex beta);

void controlRotateQubitDistributed (MultiQubit multiQubit, const int rotQubit, const int controlQubit,
                Complex rot1, Complex rot2,
                ComplexArray stateVecUp,
                ComplexArray stateVecLo,
                ComplexArray stateVecOut);

void sigmaXLocal(MultiQubit multiQubit, const int rotQubit);

void sigmaXDistributed (MultiQubit multiQubit, const int rotQubit,
                ComplexArray stateVecIn,
                ComplexArray stateVecOut);

void sigmaYLocal(MultiQubit multiQubit, const int rotQubit);

void sigmaYDistributed(MultiQubit multiQubit, const int rotQubit,
                ComplexArray stateVecIn,
                ComplexArray stateVecOut,
                int updateUpper);

void hadamardLocal (MultiQubit multiQubit, const int rotQubit);

void hadamardDistributed (MultiQubit multiQubit, const int rotQubit,
                ComplexArray stateVecUp,
                ComplexArray stateVecLo,
                ComplexArray stateVecOut, int updateUpper);


void phaseGateLocal(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type);

void phaseGateDistributed(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type);

void controlNotLocal(MultiQubit multiQubit, const int targetQubit, const int controlQubit);

void controlNotDistributed (MultiQubit multiQubit, const int targetQubit, const int controlQubit,
		        ComplexArray stateVecIn,
			ComplexArray stateVecOut);

REAL findProbabilityOfZeroLocal (MultiQubit multiQubit,
                const int measureQubit);

REAL findProbabilityOfZeroDistributed (MultiQubit multiQubit,
                const int measureQubit);

void measureInStateLocal(MultiQubit multiQubit, int measureQubit, REAL totalProbability, int outcome);

REAL measureInStateDistributedRenorm (MultiQubit multiQubit, const int measureQubit, const REAL totalProbability);

void measureInStateDistributedSetZero(MultiQubit multiQubit, const int measureQubit);

void filterOut111Local(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3,
	const REAL probOfFilter);

REAL probOfFilterOut111Local(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3);



# endif
