/** @file
 * Utility definitions for querying CPU hardware.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/errors.hpp"


#if ENABLE_MULTITHREADING && !defined(_OPENMP)
    #error "Attempted to compile in multithreaded mode without enabling OpenMP."
#endif


#if ENABLE_MULTITHREADING
    #include <omp.h>
#endif



/*
 * ENABLE OPENMP REDUCTION OF qcomp (except on MSVC compilers)
 */

#if defined(ENABLE_MULTITHREADING) && !defined(_MSC_VER)
     #pragma omp declare reduction(+ : qcomp : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#endif



/*
 * OPENMP CONFIG
 */


bool cpu_isOpenmpCompiled() {
    return (bool) ENABLE_MULTITHREADING;
}


int cpu_getCurrentNumThreads() {
#if ENABLE_MULTITHREADING
    int n = -1;

    #pragma omp parallel shared(n)
    n = omp_get_num_threads();

    return n;
#else
    error_cpuThreadsQueriedButEnvNotMultithreaded();
    return -1;
#endif
}


int cpu_getNumOpenmpProcessors() {
#if ENABLE_MULTITHREADING
    return omp_get_num_procs();
#else
    error_cpuThreadsQueriedButEnvNotMultithreaded();
    return -1;
#endif
}
