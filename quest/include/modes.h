/** @file
 * API signatures for users to specify deployment modes.
 */

#ifndef MODES_H
#define MODES_H



// ensure all mode flags are defined

#ifndef COMPILE_MPI
    #error "Compiler must define COMPILE_MPI"
#endif

#ifndef COMPILE_OPENMP
    #error "Compiler must define COMPILE_OPENMP"
#endif

#ifndef COMPILE_CUDA
    #error "Compiler must define COMPILE_CUDA"
#endif

#ifndef COMPILE_CUQUANTUM
    #error "Compiler must define COMPILE_CUQUANTUM"
#endif



// ensure all mode flags are valid values

#if ! (COMPILE_MPI == 0 || COMPILE_MPI == 1)
    #error "Macro COMPILE_MPI must have value 0 or 1"
#endif

#if ! (COMPILE_OPENMP == 0 || COMPILE_OPENMP == 1)
    #error "Macro COMPILE_OPENMP must have value 0 or 1"
#endif

#if ! (COMPILE_CUDA == 0 || COMPILE_CUDA == 1)
    #error "Macro COMPILE_CUDA must have value 0 or 1"
#endif

#if ! (COMPILE_CUQUANTUM == 0 || COMPILE_CUQUANTUM == 1)
    #error "Macro COMPILE_CUQUANTUM must have value 0 or 1"
#endif



// ensure mode flags are compatible

#if COMPILE_CUQUANTUM && ! COMPILE_CUDA
    #error "Cannot enable cuQuantum without simultaneously enabling GPU-acceleration"
#endif



// user flags for choosing automatic deployment; only accessible by C++ 
// backend and C++ users; C users must hardcode -1 

#ifdef __cplusplus

namespace modeflag { 

    extern int USE_AUTO;
}

#endif // __cplusplus



#endif // MODES_H