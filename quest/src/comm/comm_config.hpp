/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, attemptedly agnostically to
 * the implementation (like OpenMPI vs MPICH).
 */

#ifndef COMM_CONFIG_HPP
#define COMM_CONFIG_HPP



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



/*
 * MPI ENVIRONMENT MANAGEMENT
 */

bool comm_isMpiCompiled();

bool comm_isMpiGpuAware();

bool comm_isInit();

void comm_init();

void comm_end();

int comm_getRank();

int comm_getNumNodes();

void comm_sync();



#endif // COMM_CONFIG_HPP