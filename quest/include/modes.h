/** @file
 * API signatures for users to specify deployment modes.
 */

#ifndef MODES_H
#define MODES_H



// ensure all mode flags are defined

#ifndef ENABLE_DISTRIBUTION
    #error "Compiler must define ENABLE_DISTRIBUTION"
#endif

#ifndef ENABLE_MULTITHREADING
    #error "Compiler must define ENABLE_MULTITHREADING"
#endif

#ifndef ENABLE_GPU_ACCELERATION
    #error "Compiler must define ENABLE_GPU_ACCELERATION"
#endif



// ensure all mode flags are valid values

#if ! (ENABLE_DISTRIBUTION == 0 || ENABLE_DISTRIBUTION == 1)
    #error "Macro ENABLE_DISTRIBUTION must have value 0 or 1"
#endif

#if ! (ENABLE_MULTITHREADING == 0 || ENABLE_MULTITHREADING == 1)
    #error "Macro ENABLE_MULTITHREADING must have value 0 or 1"
#endif

#if ! (ENABLE_GPU_ACCELERATION == 0 || ENABLE_GPU_ACCELERATION == 1)
    #error "Macro ENABLE_GPU_ACCELERATION must have value 0 or 1"
#endif



// user flags for choosing automatic deployment; only accessible by C++ 
// backend and C++ users; C users must hardcode -1 

#ifdef __cplusplus

namespace modeflag { 

    extern int USE_AUTO;
}

#endif // __cplusplus



#endif // MODES_H