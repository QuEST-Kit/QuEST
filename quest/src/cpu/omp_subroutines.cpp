/** @file
 * CPU OpenMP-accelerated definitions of the subroutines called by
 * accelerator.cpp. 
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"

#if ENABLE_MULTITHREADING
    #include <omp.h>
#endif


// inform OpenMP how to reduce qcomp instances (except on MSVC compilers)
#if defined(ENABLE_MULTITHREADING) && !defined(_MSC_VER)
     #pragma omp declare reduction(+ : qcomp : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#endif