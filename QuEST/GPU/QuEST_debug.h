// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

# ifndef QuEST_DEBUG
# define QuEST_DEBUG

# include "QuEST_precision.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @file
 * Developer functions used for unit testing and debugging. Not part of the public API. May contain
 * functions that are incomplete or untested.
 */

void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome);

void initStateDebug(MultiQubit *multiQubit);

void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env);

int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision);

/** Report a list of CPU hostnames and the rank that is running on each if running with MPI enabled and an 
 * error message otherwise. For debugging purposes. 
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * */
void reportNodeList(QuESTEnv env);

#ifdef __cplusplus
}
#endif

# endif
