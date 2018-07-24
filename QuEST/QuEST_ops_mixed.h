// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file
 * Functions operating upon mixed states which are provided by a hardware-specific backend
 */

# ifndef QuEST_OPS_MIXED
# define QuEST_OPS_MIXED

# include "QuEST.h"
# include "QuEST_precision.h"

# ifdef __cplusplus
extern "C" {
# endif

void mixed_initStatePlus(QubitRegister targetQureg);

void mixed_initClassicalState(QubitRegister qureg, long long int stateInd);

void mixed_initPureState(QubitRegister targetQureg, QubitRegister copyQureg);

REAL mixed_calcTotalProbability(QubitRegister qureg);

REAL mixed_findProbabilityOfOutcome(QubitRegister qureg, const int measureQubit, int outcome);


# ifdef __cplusplus
}
# endif

# endif //QuEST_OPS_MIXED