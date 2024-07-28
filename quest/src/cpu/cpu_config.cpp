/** @file
 * Utility definitions for querying CPU hardware.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/errors.hpp"


#if COMPILE_OPENMP && !defined(_OPENMP)
    #error "Attempted to compile in multithreaded mode without enabling OpenMP in the compiler flags."
#endif


#if COMPILE_OPENMP
    #include <omp.h>
#endif



/*
 * ENABLE OPENMP REDUCTION OF qcomp (except on MSVC compilers)
 */

#if defined(COMPILE_OPENMP) && !defined(_MSC_VER)
     #pragma omp declare reduction(+ : qcomp : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#endif



/*
 * OPENMP CONFIG
 */


bool cpu_isOpenmpCompiled() {
    return (bool) COMPILE_OPENMP;
}


int cpu_getCurrentNumThreads() {
#if COMPILE_OPENMP
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
#if COMPILE_OPENMP
    return omp_get_num_procs();
#else
    error_cpuThreadsQueriedButEnvNotMultithreaded();
    return -1;
#endif
}



/*
 * MEMORY MANAGEMENT
 */


qcomp* cpu_allocAmps(qindex numLocalAmps) {

    // we call calloc over malloc in order to fail immediately if mem isn't available;
    // caller must handle NULL result
    return (qcomp*) calloc(numLocalAmps, sizeof(qcomp));
}


void cpu_deallocAmps(qcomp* amps) {

    free(amps);
}