/** @file
 * Functions for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"

#if ENABLE_DISTRIBUTION
    #include <mpi.h>
#endif



/*
 * MPI COMPLEX TYPE FLAG
 */

#if ENABLE_DISTRIBUTION

    #if (FLOAT_PRECISION == 1)
        #define MPI_QCOMP MPI_CXX_FLOAT_COMPLEX

    #elif (FLOAT_PRECISION == 2)
        #define MPI_QCOMP MPI_CXX_DOUBLE_COMPLEX

    // sometimes 'MPI_CXX_LONG_DOUBLE_COMPLEX' isn't defined
    #elif (FLOAT_PRECISION == 4) && defined(MPI_CXX_LONG_DOUBLE_COMPLEX)
        #define MPI_QCOMP MPI_CXX_LONG_DOUBLE_COMPLEX

    // in that case, fall back to the C type (identical memory layout)
    #else
        #define MPI_QCOMP MPI_C_LONG_DOUBLE_COMPLEX
    #endif

#endif