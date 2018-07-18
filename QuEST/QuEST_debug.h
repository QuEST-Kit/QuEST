// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Developer functions used for unit testing and debugging, which are not part of the public API. 
 * May contain functions that are incomplete or untested.
 */

# ifndef QuEST_DEBUG
# define QuEST_DEBUG

# include "QuEST_precision.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialise the state vector of probability amplitudes such that one qubit is set to 'outcome' and all other qubits are in an equal superposition of zero and one.
 * @param[in,out] qureg object representing the set of qubits to be initialised
 * @param[in] qubitId id of qubit to set to state 'outcome'
 * @param[in] value of qubit 'qubitId'
 */
void initStateOfSingleQubit(QubitRegister *qureg, int qubitId, int outcome);

/**
 * Initialise the state vector of probability amplitudes to an (unphysical) state with
 * each component of each probability amplitude a unique floating point value. For debugging processes
 * @param[in,out] qureg object representing the set of qubits to be initialised
 */
void initStateDebug(QubitRegister qureg);

/** Initialises the wavefunction amplitudes according to those specified in a file.
 * For debugging purpsoses 
 */
void initializeStateFromSingleFile(QubitRegister *qureg, char filename[200], QuESTEnv env);

/** Return whether two given wavefunctions are equivalent within a given precision
 * Global phase included in equivalence check. For debugging purposes.
 */
int compareStates(QubitRegister mq1, QubitRegister mq2, REAL precision);

/** Report a list of CPU hostnames and the rank that is running on each if running with MPI enabled and an 
 * error message otherwise. For debugging purposes. 
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * */
void reportNodeList(QuESTEnv env);

#ifdef __cplusplus
}
#endif

# endif
