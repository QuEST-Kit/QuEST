// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * GPU backend routines which are invoked by both the cuQuantum and bespoke-kernel backends.
 * This excludes backend-agnostic API implementations which do not need a header declaration.
 *
 * @author Tyson Jones
 */

# ifndef QUEST_GPU_COMMON_H
# define QUEST_GPU_COMMON_H

# ifdef __cplusplus
extern "C" {
# endif



int GPUExists(void);



#ifdef __cplusplus
}
#endif

# endif // QUEST_GPU_COMMON_H