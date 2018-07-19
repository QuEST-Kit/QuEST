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
	

void mixed_initPureState(QubitRegister targetQureg, QubitRegister copyQureg);


	
# ifdef __cplusplus
}
# endif

# endif //QuEST_OPS_MIXED