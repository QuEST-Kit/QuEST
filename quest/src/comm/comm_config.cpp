/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, agnostically to the specific
 * implementation (like OpenMPI vs MPICH). These functions
 * are callable even when MPI has not been compiled/linked.
 * 
 * Note that even when QUEST_COMPILE_MPI=1, the user may have
 * disabled distribution when creating the QuEST environment
 * at runtime - even despite they themselves initialising and
 * using MPI. So we must be careful about consulting MPI status!
 * Furthermore, all routines here will only ever consult/affect
 * the QuEST communicator, never the entire MPI environment,
 * the latter of which may contain non-participating processes.
 * 
 * @author Tyson Jones
 */

#include "quest/include/config.h"
#include "quest/include/types.h"

#include "quest/src/comm/comm_config.hpp"
#include "quest/src/core/errors.hpp"

#include <string>

#if QUEST_COMPILE_MPI
    #include <mpi.h>
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
 * COMMUNICATOR MANAGEMENT
 *
 * QuEST will only ever use the overridable global_mpiComm communicator,
 * so that superusers can dedicate external MPI processes to other tasks.
 * Beware that it's valid for QuEST to be compiled with MPI, but have
 * distribution runtime-disabled, while the user is themselves using
 * (and ergo have initialised) MPI. In that scenario, we must not touch
 * MPI, hence why comm_isActive() below is distinct from comm_isMpiInit().
 */


// We must record whether the user owns MPI, so that we do not ever attempt
// to kill it when gracefully exiting, or due to a validation error
static bool global_isMpiUserOwned = false;


// Guarded since MPI_Comm cannot be exposed when not compiling MPI. This
// communicator is overridden from NULL either BEFORE or DURING comm_init()
#if QUEST_COMPILE_MPI
    static MPI_Comm global_mpiComm = MPI_COMM_NULL;
#endif


bool comm_isActive() {
#if QUEST_COMPILE_MPI

    // comm_init(), or potentially comm_setMpiComm() before it, will only
    // ever override mpiComm with non-NULL, indicating active comm. Note
    // it's principally for mpiComm to later return to NULL, via comm_end(),
    // and for QuEST execution to continue (though not supported presently).
    // if comm_isActive() is true, then it is guaranteed MPI is initialised
    return global_mpiComm != MPI_COMM_NULL;

    // note it is legal for QuEST distribution to be disabled (and ergo
    // mpiComm never initialised) even when the user is themselves accessing
    // MPI, hence this function is semantically distinct from comm_isMpiInit()
#else

    // QuEST communication is obviously never active if
    // not even MPI is compiled; though this does not
    // imply at all the user isn't themselves using MPI!
    return false;

#endif
}


// Hide MPI_Comm from signatures when MPI is not compiled. Beware that
// these are not exposed in comm_config.hpp; callers must 'extern' them!
#if QUEST_COMPILE_MPI


MPI_Comm comm_getMpiComm() {

    // illegal to call before communicator has been overridden
    if (global_mpiComm == MPI_COMM_NULL)
        error_commMpiCommIsNull();

    return global_mpiComm;
}


bool comm_setMpiComm(MPI_Comm newComm, bool userOwnsMpi) {

    // illegal to re-set, or set to null
    if (global_mpiComm != MPI_COMM_NULL)
        error_commAlreadyHasSetMpiComm();
    if (newComm == MPI_COMM_NULL)
        error_commNewMpiCommIsNull();

    // detect bad communicator, and inform validation
    auto status = MPI_Comm_dup(newComm, &global_mpiComm);
    if (status != MPI_SUCCESS)
        return false;

    // record ownership as soon as QuEST communication becomes active, so
    // validation errors during env initialisation never kill user-owned MPI
    global_isMpiUserOwned = userOwnsMpi;
    return true;
}


#endif // QUEST_COMPILE_MPI



/*
 * MPI ENVIRONMENT MANAGEMENT
 *
 * which queries MPI itself (as may be user-activated), rather
 * than QuEST's (possibly more limited) MPI environment
 */


bool comm_isMpiCompiled() {
    return (bool) QUEST_COMPILE_MPI;
}

bool comm_isMpiSubCommCompiled() {
    return (bool) QUEST_COMPILE_SUBCOMM;
}


bool comm_isMpiGpuAware() {

    // well duh
    if (!comm_isMpiCompiled())
        return false;

    // definitely not GPU-aware if compiler declares it is not
    #if defined(MPIX_CUDA_AWARE_SUPPORT) && ! MPIX_CUDA_AWARE_SUPPORT
        return false;
    #endif

    // check CUDA-awareness at run-time if we know it's principally supported
    #if defined(MPIX_CUDA_AWARE_SUPPORT)
        return (bool) MPIX_Query_cuda_support();
    #endif

    // check whether an MPICH env-var indicates support (we assume it never lies!)
    static const auto var = std::getenv("MPICH_GPU_SUPPORT_ENABLED");
    if (var && std::string(var) == "1") // ill-formed vars = 0
        return true;

    // if we can't ascertain CUDA-awareness, just assume no to avoid seg-fault
    return false;
}


bool comm_isMpiInit() {
#if QUEST_COMPILE_MPI

    // safely callable before MPI initialisation, but NOT after comm_end()
    int isInit;
    MPI_Initialized(&isInit);

    // when MPI is not initialised, it is guaranteed that QuEST's communicator
    // is inactive, which we double check here so callers can be absolutely sure
    if (!isInit && comm_isActive())
        error_commActiveButMpiNotInit();

    return (bool) isInit;

#else

    // obviously MPI is never initialised if not even compiled
    return false;

#endif
}


bool comm_isMpiUserOwned() {

    // this isn't presently used by the code base; I'm just naughtily silencing
    // "unused var" warning when compiling without MPI :^)
    return global_isMpiUserOwned;
}



/*
 * QUEST COMMUNICATION MANAGEMENT
 *
 * which interacts only with QuEST's MPI environment,
 * which may be smaller than the user-controlled MPI env
 */


void comm_init(bool userOwnsMpi) {
#if QUEST_COMPILE_MPI

    // re-assert prior user-validations for clarity
    if (userOwnsMpi && !comm_isMpiInit())
        error_commNotInit();
    if (!userOwnsMpi && comm_isMpiInit())
        error_commAlreadyInit();
   
    // init MPI only when it's not the user's responsibility
    if (!userOwnsMpi)
        MPI_Init(NULL, NULL);

    // choose communicator only when the user hasn't already
    // (via comm_setMpiComm, during custom env initialisation)
    if (global_mpiComm == MPI_COMM_NULL)
        comm_setMpiComm(MPI_COMM_WORLD, userOwnsMpi);

#endif
}


void comm_end() {
#if QUEST_COMPILE_MPI

    // If QuEST isn't using distribution, regardless of whether the user is using MPI,
    // then we gracefully exit. We do NOT attempt to end MPI on the user's behalf (as we
    // may be tempted to do during validation failure to avoid their MPI-crash), because
    // it's possible/legal that not all processes are participating in this comm_end()
    // call, in which case so MPI_Finalize() could just cause a hang.
    if (!comm_isActive())
        return;

    // Syncing is not strictly necessary, but it ensures that finalizeQuESTEnv() never
    // completes on one process while another process is still performing simulation
    // (though that'd be weird), and so may avoid a silly user benchmarking pitfall
    MPI_Barrier(global_mpiComm);
    MPI_Comm_free(&global_mpiComm);
    
    // Do NOT close MPI if the user owns; they may still wish to use it after QuEST!
    if (!global_isMpiUserOwned)
        MPI_Finalize();

    // Presently, comm_end() is only ever called during QuESTEnv destruction (either
    // deliberately, or because of failed validation during QuESTEnv initialisation).
    // This means any comm_*() call hereafter is invalid/illegal and will be prevented
    // by validation. However, we can imagine a future where distribution gets runtime
    // disabled while QuEST execution continues (e.g. initQuESTEnv automatically
    // disabled distribution), and so we must indicate that communication is no longer
    // active by overwriting comm to NULL. BEWARE that this is "hacky"; we have
    // updated mpiComm here without MPI_Comm_dup(), but that's fine, because hereafter
    // MPI will never be used again (illegal to re-init both MPI, and QuEST!)
    global_mpiComm = MPI_COMM_NULL;
    global_isMpiUserOwned = false;

#endif
}


int comm_getRank() {
#if QUEST_COMPILE_MPI

    // if distribution was not runtime enabled (or a validation error was 
    // triggered during distributed initialisation), every process believes
    // it is the root rank; this may lead to unavoidable error msg spam!
    if (!comm_isActive())
        return ROOT_RANK;

    // obtain the process rank within the QuEST communicator, which can
    // differ from the global MPI process rank when users own MPI
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
    // triggered during distributed initialisation), every process is told
    // it is the one and only node; this may lead to error msg spam, but
    // appears unavoidable!
    if (!comm_isActive())
        return 1;

    // obtain the number of processes within the QuEST communicator, which
    // can be smaller than global MPI process count when users own MPI
    int numNodes;
    MPI_Comm_size(global_mpiComm, &numNodes);
    return numNodes;

#else

    // if MPI isn't compiled, QuEST is definitely non-distributed and
    // each process only knows itself (though users may own MPI and
    // actually have many processes; that's none of our business!)
    return 1;

#endif
}


void comm_sync() {
#if QUEST_COMPILE_MPI

    // gracefully handle when not distributed, needed by e.g. pre-MPI-setup validation
    if (!comm_isActive())
        return;

    MPI_Barrier(global_mpiComm);

#endif

    // do nothing at all when MPI is not compiled (user owned MPI processes go unsynced)
}
