/** @file
 * API signatures for debugging QuEST behaviour, or controlling
 * input validation.
 */

#ifndef DEBUG_H
#define DEBUG_H

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



// TODO:
// is this really the right place to put CUQUANTUM_MEM_POOL_BYTES??????


// sets the size at which the CUDA memory pool will 
// automatically deallocate temporary memory. Below this
// size, temporary memory structures (like a CompMatr)
// will persist in GPU memory to save time. This is
// only relevant to GPU-mode with cuQuantum enabled,
// and is effected at createQuESTEnv().
int CUQUANTUM_MEM_POOL_BYTES = 16*(1<<15); // 1 MiB ~ 8 qubit complex<double> matrix


void invalidQuESTInputError(const char* msg, const char* func);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // DEBUG_H