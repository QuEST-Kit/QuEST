/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, agnostically to the specific
 * implementation (like OpenMPI vs MPICH). These functions
 * are callable even when MPI has not been compiled/linked.
 * 
 * Note that even when COMPILE_MPI=1, the user may have
 * disabled distribution when creating the QuEST environment
 * at runtime. Ergo we use comm_isInit() to determine whether
 * functions should invoke the MPI API.
 * 
 * @author Tyson Jones
 */

#include "quest/include/config.h"
#include "quest/include/types.h"

#include "quest/src/comm/comm_config.hpp"
#include "quest/src/core/errors.hpp"

#if COMPILE_MPI
    #include <mpi.h>

    static MPI_Comm mpiCommQuest = MPI_COMM_NULL;
#endif



/*
 * WARN ABOUT CUDA-AWARENESS
 */

#if COMPILE_MPI && COMPILE_CUDA

    // this check is OpenMPI specific
    #ifdef OPEN_MPI
        #include <mpi-ext.h>

        // #warning command is always recognised (OpenMPI is not Windows compatible)
        #ifndef MPIX_CUDA_AWARE_SUPPORT
            #warning "Could not ascertain whether MPI is CUDA-aware, so we will assume it is not. This means inter-GPU communication will be slowly routed through the CPU/RAM."
        #elif !MPIX_CUDA_AWARE_SUPPORT
            #warning "MPI compiler is not CUDA-aware, so inter-GPU communication will be slowly routed through the CPU/RAM"
        #endif
    #endif

    /// @todo check whether MPICH is CUDA-aware
    /// beware MSVC cannot parse #warning, and
    /// Intel MPI would crash (but not MSMPI?)

#endif



/*
 * MPI ENVIRONMENT MANAGEMENT
 * all of which is safely callable in non-distributed mode
 */


bool comm_isMpiCompiled() {
    return (bool) COMPILE_MPI;
}

bool comm_isMpiSubCommunicatorCompiled() {
    return (bool) COMPILE_SUBCOMM;
}

bool comm_isMpiGpuAware() {

    /// @todo these checks may be OpenMPI specific, so that
    /// non-OpenMPI MPI compilers are always dismissed as
    /// not being CUDA-aware. Check e.g. MPICH method!

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
#if COMPILE_MPI

    // safely callable before MPI initialisation, but NOT after comm_end()
    int isInit;
    MPI_Initialized(&isInit);
    return (bool) isInit;

#else

    // obviously MPI is never initialised if not even compiled
    return false;
#endif
}


void comm_init(int userOwnsMpi) {
#if COMPILE_MPI

    // error if attempting re-initialisation
    if (!userOwnsMpi && comm_isInit())
        error_commAlreadyInit();

    // error if user owns MPI but has not initialised
    if (userOwnsMpi && !comm_isInit()) {
        error_commNotInit();
    }
   
    // QuEST must initialise MPI if the user does not own it
    if (!userOwnsMpi)
        MPI_Init(NULL, NULL);
    
    // If user is setting their own comm, mpiCommQuest will be NOT MPI_COMM_NULL, 
    // and we should not touch it.
    // If user is NOT setting their own comm, mpiCommQuest will be MPI_COMM_NULL,
    // and we should set it to MPI_COMM_WORLD.
    if (mpiCommQuest == MPI_COMM_NULL)
        MPI_Comm_dup(MPI_COMM_WORLD, &mpiCommQuest);

#endif
}


void comm_end(int userOwnsMpi) {
#if COMPILE_MPI

    // gracefully permit comm_end() before comm_init(), as input validation can trigger
    if (!comm_isInit())
        return;

    MPI_Barrier(mpiCommQuest);
    MPI_Comm_free(&mpiCommQuest);
    
    // QuEST must finalise MPI if the user does not own it
    if (!userOwnsMpi)
        MPI_Finalize();

#endif
}


int comm_getRank() {
#if COMPILE_MPI

    // if distribution was not runtime enabled (or a validation error was 
    // triggered), every node (if many MPI processes were launched)
    // believes it is the root rank
    if (!comm_isInit())
        return ROOT_RANK;

    int rank;
    MPI_Comm_rank(mpiCommQuest, &rank);
    return rank;

#else

    // if MPI isn't compiled, we're definitely non-distributed; return main rank 
    return ROOT_RANK;
#endif
}


bool comm_isRootNode(int rank) {
    return rank == ROOT_RANK;
}
bool comm_isRootNode() {
    return comm_isRootNode(comm_getRank());
}


int comm_getNumNodes() {
#if COMPILE_MPI

    // if distribution was not runtime enabled (or a validation error was 
    // triggered), every node (if many MPI processes were launched)
    // believes it is the one and only node
    if (!comm_isInit())
        return 1;

    int numNodes;
    MPI_Comm_size(mpiCommQuest, &numNodes);
    return numNodes;

#else

    // if MPI isn't compiled, we're definitely non-distributed; return single node
    return 1;
#endif
}


void comm_sync() {
#if COMPILE_MPI

    // gracefully handle when not distributed, needed by e.g. pre-MPI-setup validation 
    if (!comm_isInit())
        return;

    MPI_Barrier(mpiCommQuest);
#endif
}

#if COMPILE_MPI
    MPI_Comm comm_getMpiComm() {
        return mpiCommQuest;
    }

    #if COMPILE_SUBCOMM
        void comm_setMpiComm(MPI_Comm newComm) {
            // error if mpiCommQuEST is already set!
            if (mpiCommQuest != MPI_COMM_NULL) {
                error_commDoubleSetMpiComm();
            }

            MPI_Comm_dup(newComm, &mpiCommQuest);

            return;
        }
    #endif
#endif
