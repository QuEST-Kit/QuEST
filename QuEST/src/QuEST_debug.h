// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Developer functions used for unit testing and debugging, which are not part of the public API. 
 * May contain functions that are incomplete or untested.
 *
 * @author Ania Brown
 */

# ifndef QUEST_DEBUG_H
# define QUEST_DEBUG_H

# include "QuEST_precision.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Initialise the state vector of probability amplitudes such that one qubit is set to 'outcome' and all other qubits are in an equal superposition of zero and one.
 * @param[in,out] qureg object representing the set of qubits to be initialised
 * @param[in] qubitId id of qubit to set to state 'outcome'
 * @param[in] outcome value of qubit 'qubitId' to set
 */
void initStateOfSingleQubit(Qureg *qureg, int qubitId, int outcome);

/**
 * Initialise the state vector of probability amplitudes to an (unphysical) state with
 * each component of each probability amplitude a unique floating point value. For debugging processes
 * @param[in,out] qureg object representing the set of qubits to be initialised
 */
void initStateDebug(Qureg qureg);

/** Initialises the wavefunction amplitudes according to those specified in a file.
 * For debugging purpsoses 
 */
void initStateFromSingleFile(Qureg *qureg, char filename[200], QuESTEnv env);

/** Return whether two given wavefunctions are equivalent within a given precision
 * Global phase included in equivalence check. For debugging purposes.
 */
int compareStates(Qureg mq1, Qureg mq2, qreal precision);

/** Set elements in the underlying state vector represenation of a density matrix. Not exposed in the public
 * API as this requires an understanding of how the state vector is used to represent a density matrix.
 * Currently can only be used to set all amps. 
 */
void setDensityAmps(Qureg qureg, qreal* reals, qreal* imags);



/** Return the precision of qreal for use in testing 
 * */
  
int QuESTPrecision(void);
  
#ifdef __cplusplus
}
#endif

# endif // QUEST_DEBUG_H
