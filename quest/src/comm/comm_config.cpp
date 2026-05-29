/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, agnostically to the specific
 * implementation (like OpenMPI vs MPICH). These functions
 * are callable even when MPI has not been compiled/linked.
 * 
 * Note that even when QUEST_COMPILE_MPI=1, the user may have
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

#if QUEST_COMPILE_MPI
    #include <mpi.h>

    static MPI_Comm global_mpiComm = MPI_COMM_NULL;
#endif



/*
 * WARN ABOUT CUDA-AWARENESS
 */

#if QUEST_COMPILE_MPI && QUEST_COMPILE_CUDA

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
    return (bool) QUEST_COMPILE_MPI;
}

bool comm_isMpiSubCommunicatorCompiled() {
    return (bool) QUEST_COMPILE_SUBCOMM;
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
#if QUEST_COMPILE_MPI

    // safely callable before MPI initialisation, but NOT after comm_end()
    int isInit;
    MPI_Initialized(&isInit);
    return (bool) isInit;

#else

    // obviously MPI is never initialised if not even compiled
    return false;
#endif
}


void comm_init(bool userOwnsMpi) {
#if QUEST_COMPILE_MPI

    // re-assert prior user-validations for robustness
    if (userOwnsMpi && !comm_isInit())
        error_commNotInit();
    if (!userOwnsMpi && comm_isInit())
        error_commAlreadyInit();
   
    // init MPI only when it's not the user's responsibility
    if (!userOwnsMpi)
        MPI_Init(NULL, NULL);

    // choose communicator only when the user hasn't 
    if (global_mpiComm == MPI_COMM_NULL)
        MPI_Comm_dup(MPI_COMM_WORLD, &global_mpiComm);

#endif
}



void comm_end(bool userOwnsMpi) {
#if QUEST_COMPILE_MPI

    // gracefully permit comm_end() before comm_init(), as input validation can trigger
    if (!comm_isInit())
        return;

    // ungracefully handle when the communicator is still NULL, because comm_end() may be
    // triggered by "bad MPI init" validation, during which, the communicator may not yet
    // have been set. We choose NOT to divert to MPI_COMM_WORLD, which is likely just to
    // stall at MPI_Barrier, and instead let the user's communicator live on; then crash!
    if (global_mpiComm == MPI_COMM_NULL)
        return;

    MPI_Barrier(global_mpiComm);
    MPI_Comm_free(&global_mpiComm);
    
    // QuEST must finalise MPI if the user does not own it
    if (!userOwnsMpi)
        MPI_Finalize();

#endif
}


int comm_getRank() {
#if QUEST_COMPILE_MPI

    // if distribution was not runtime enabled (or a validation error was 
    // triggered), every node (if many MPI processes were launched)
    // believes it is the root rank
    if (!comm_isInit())
        return ROOT_RANK;

    // consult the (potentially sub-) communicator for rank; if it is still
    // NULL, as can only validly happen during failed MPI status validation (the
    // error msg is attemptedly printed on only the root process), fallback to
    // using WORLD (and pray the user hasn't silenced world-root std-out!). We
    // COULD safely return ROOT_RANK instead, letting all processes believe they
    // are root, but this grossly duplicates the output across ALL processes
    MPI_Comm comm = (global_mpiComm == MPI_COMM_NULL)? MPI_COMM_WORLD : global_mpiComm;

    int rank;
    MPI_Comm_rank(comm, &rank);
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
#if QUEST_COMPILE_MPI

    // if distribution was not runtime enabled (or a validation error was 
    // triggered), every node (if many MPI processes were launched)
    // believes it is the one and only node
    if (!comm_isInit())
        return 1;

    int numNodes;
    MPI_Comm_size(global_mpiComm, &numNodes);
    return numNodes;

#else

    // if MPI isn't compiled, we're definitely non-distributed; return single node
    return 1;
#endif
}


void comm_sync() {
#if QUEST_COMPILE_MPI

    // gracefully handle when not distributed, needed by e.g. pre-MPI-setup validation 
    if (!comm_isInit())
        return;

    // gracefully handle when the communicator is still NULL, because comm_sync() is
    // triggered by "bad MPI init" validation (during the error message printing)
    // during which, the communicator may not yet have been overriden
    if (global_mpiComm == MPI_COMM_NULL)
        return;

    MPI_Barrier(global_mpiComm);
#endif
}

#if QUEST_COMPILE_MPI
    MPI_Comm comm_getMpiComm() {

        if (global_mpiComm == MPI_COMM_NULL)
            error_commMpiCommIsNull();

        return global_mpiComm;
    }

    #if QUEST_COMPILE_SUBCOMM
        void comm_setMpiComm(MPI_Comm newComm) {

            // error if global_mpiComm is already set!
            if (global_mpiComm != MPI_COMM_NULL) {
                MPI_Barrier(global_mpiComm);
                MPI_Comm_free(&global_mpiComm);
                error_commDoubleSetMpiComm();
            }

            int mpi_err = MPI_Comm_dup(newComm, &global_mpiComm);
            if (mpi_err != MPI_SUCCESS) {
                error_commInvalidMpiComm();
            }

            return;
        }
    #endif
#endif
