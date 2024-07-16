/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, attemptedly agnostically to
 * the implementation (like OpenMPI vs MPICH).
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/errors.hpp"

#if ENABLE_DISTRIBUTION
    #include <mpi.h>
#endif



/*
 * WARN ABOUT CUDA-AWARENESS
 */

#if ENABLE_DISTRIBUTION && ENABLE_GPU_ACCELERATION

    // this check is OpenMPI specific
    #ifdef OPEN_MPI
        #include <mpi-ext.h>

        #ifndef MPIX_CUDA_AWARE_SUPPORT
            #warning "Could not ascertain whether MPI is CUDA-aware, so we will assume it is not. This means inter-GPU communication will be slowly routed through the CPU/RAM."
        #elif !MPIX_CUDA_AWARE_SUPPORT
            #warning "MPI compiler is not CUDA-aware, so inter-GPU communication will be slowly routed through the CPU/RAM"
        #endif
    #endif

    // TODO: check whether MPICH is CUDA-aware

#endif



/*
 * MPI ENVIRONMENT MANAGEMENT
 * all of which is safely callable in non-distributed mode
 */


bool comm_isMpiCompiled() {
    return (bool) ENABLE_DISTRIBUTION;
}


bool comm_isMpiGpuAware() {

    // TODO: these checks may be OpenMPI specific, so that
    // non-OpenMPI MPI compilers are always dismissed as
    // not being CUDA-aware. Check e.g. MPICH method!

    // definitely not GPU-aware if compiler declares it is not
    #if defined(MPIX_CUDA_AWARE_SUPPORT) && ! MPIX_CUDA_AWARE_SUPPORT
        return false;
    #endif

    // check CUDA-awareness at run-time if we know it's principally supported
    #if defined(MPIX_CUDA_AWARE_SUPPORT)
        return (bool) MPIX_Query_cuda_support();
    #endif

    // if we can't ascertain CUDA-awareness, just assume no to avoid seg-fault
    return false;
}


bool comm_isInit() {
#if ENABLE_DISTRIBUTION

    // safely callable before MPI initialisation, but NOT after comm_end()
    int isInit;
    MPI_Initialized(&isInit);
    return (bool) isInit;

#else

    // obviously MPI is never initialised if not even compiled
    return false;
#endif
}


void comm_init() {
#if ENABLE_DISTRIBUTION

    // error if attempting re-initialisation
    if (comm_isInit())
        error_commAlreadyInit();
    
    MPI_Init(NULL, NULL);
#endif
}


void comm_end() {
#if ENABLE_DISTRIBUTION

    // gracefully permit comm_end() before comm_init(), as input validation can trigger
    if (!comm_isInit())
        return;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

#endif
}


int comm_getRank() {
#if ENABLE_DISTRIBUTION

    // if MPI not yet setup (e.g. QuESTEnv creation error'd), return main rank
    if (!comm_isInit())
        return 0;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;

#else

    // if MPI isn't compiled, we're definitely non-distributed; return main rank 
    return 0;
#endif
}


bool comm_isRootNode(int rank) {
    return rank == 0;
}
bool comm_isRootNode() {
    return comm_isRootNode(comm_getRank());
}


int comm_getNumNodes() {
#if ENABLE_DISTRIBUTION

    // if MPI not yet setup, error; else we may misreport later MPI env
    if (!comm_isInit())
        error_commNotInit();

    int numNodes;
    MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
    return numNodes;

#else

    // if MPI isn't compiled, we're definitely non-distributed; return single node
    return 1;
#endif
}


void comm_sync() {
#if ENABLE_DISTRIBUTION

    // gracefully handle when not distributed, needed by pre-MPI-setup validation 
    if (!comm_isInit())
        return;

    MPI_Barrier(MPI_COMM_WORLD);
#endif
}