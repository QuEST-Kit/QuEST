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
 *
 * all of which is safely callable in non-distributed mode
 */


bool comm_isMpiCompiled() {
    return (bool) QUEST_COMPILE_MPI;
}

bool comm_isMpiSubCommCompiled() {
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

    // gracefully handle when the communicator is still NULL, because comm_end() may be
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

    // Consult the (potentially sub-) communicator for rank; if it is still
    // NULL, as can only validly happen during failed QuESTEnv init validation
    // (which triggers root-only error printing and ergo this function), we
    // fall back to every process believing it is root and so attempting to
    // print. This safely avoids consulting a potentially bugged MPI communicator
    // and losing the message. We once tried to fallback to MPI_COMM_WORLD here,
    // to avoid duplicate output, but it is not worth the risk of msg loss!
    if (global_mpiComm == MPI_COMM_NULL)
        return ROOT_RANK;

    int rank;
    MPI_Comm_rank(global_mpiComm, &rank);
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



/*
 * MPI COMMUNICATOR MANAGEMENT
 *
 * some of which requires exposing MPI_Comm in external-facing signatures.
 * In lieu of leaking these into comm_config.hpp, callers must extern them.
 */

bool comm_isMpiCommSet() {
#if QUEST_COMPILE_MPI

    // once comm_init() or comm_setMpiComm() overwrite
    // the communicator, is can never return to NULL  
    return (global_mpiComm != MPI_COMM_NULL);
# else
    return false;
#endif
}

#if QUEST_COMPILE_MPI

MPI_Comm comm_getMpiComm() {

    if (global_mpiComm == MPI_COMM_NULL)
        error_commMpiCommIsNull();

    return global_mpiComm;
}

bool comm_setMpiComm(MPI_Comm newComm) {

    // this is called prior to QuEST initialisation,
    // and merely seeks to overwrite global_mpiComm 

    if (global_mpiComm != MPI_COMM_NULL)
        error_commAlreadyHasSetMpiComm();
    if (newComm == MPI_COMM_NULL)
        error_commMpiCommIsNull();

    auto status = MPI_Comm_dup(newComm, &global_mpiComm);
    return status == MPI_SUCCESS;
}

#endif // QUEST_COMPILE_MPI
