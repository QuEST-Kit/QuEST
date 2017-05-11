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

REAL findProbabilityOfZeroLocal (MultiQubit multiQubit,
                const int measureQubit);

REAL findProbabilityOfZeroDistributed (MultiQubit multiQubit,
                const int measureQubit);

void measureInZeroLocal(MultiQubit multiQubit, int measureQubit, REAL totalProbability);

REAL measureInZeroDistributedRenorm (MultiQubit multiQubit, const int measureQubit, const REAL totalProbability);

void measureInZeroDistributedSetZero(MultiQubit multiQubit, const int measureQubit);

void filterOut111Local(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3,
	const REAL probOfFilter);

REAL probOfFilterOut111Local(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3);



# endif
