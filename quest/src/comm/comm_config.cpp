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

#include <string>

#if QUEST_COMPILE_MPI
    #include <mpi.h>

    static MPI_Comm mpiCommQuest = MPI_COMM_NULL;
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


void comm_init(int useDistrib, bool userOwnsMpi) {
#if QUEST_COMPILE_MPI

    // error if user owns MPI but has not initialised
    if (userOwnsMpi && !comm_isInit()) {
        error_commNotInit();
    }
   
    // Overall mpiCommQuest should be set in the following ways
    // however only useDistrib = 1 and userOwnsMpi = false
    // and useDistrib = 0 and userOwnsMpi = true 
    // require action here
    //
    // | useDistrib | userOwnsMpi |  mpiCommQuest  |
    // | ---------- | ----------- | -------------- |
    // |     0      |    false    | MPI_COMM_NULL  |
    // | ---------- | ----------- | -------------- |
    // |     1      |    false    | MPI_COMM_WORLD |
    // | ---------- | ----------- | -------------- |
    // |     0      |    true     | MPI_COMM_SELF  |
    // | ---------- | ----------- | -------------- |
    // |            |             | MPI_COMM_WORLD |
    // |     1      |    true     |      or        |
    // |            |             | userQuestComm  |
    // | ---------- | ----------- | -------------- |
    

    if (useDistrib && !userOwnsMpi) {
        // error if attempting re-initialisation
        if (comm_isInit()) {
            error_commAlreadyInit();
        } else {
            MPI_Init(NULL, NULL);
            // The user wants MPI and is leaving it to QuEST
            MPI_Comm_dup(MPI_COMM_WORLD, &mpiCommQuest);
        }
    } else if (!useDistrib && userOwnsMpi) {
        // The user has initialised MPI but wants QuEST to ignore it
        MPI_Comm_dup(MPI_COMM_SELF, &mpiCommQuest);
    } else if (useDistrib && userOwnsMpi) {
        // if mpiCommQuEST is still MPI_COMM_NULL the user is not 
        // providing their own MPI_Comm and we should set mpiCommQuest
        // to MPI_COMM_WORLD
        if (mpiCommQuest == MPI_COMM_NULL)
            MPI_Comm_dup(MPI_COMM_WORLD, &mpiCommQuest);
    }

#endif
    return;
}



void comm_end(bool userOwnsMpi) {
#if QUEST_COMPILE_MPI

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
#if QUEST_COMPILE_MPI

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
#if QUEST_COMPILE_MPI

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
#if QUEST_COMPILE_MPI

    // gracefully handle when not distributed, needed by e.g. pre-MPI-setup validation 
    if (!comm_isInit())
        return;

    MPI_Barrier(mpiCommQuest);
#endif
}

#if QUEST_COMPILE_MPI
    MPI_Comm comm_getMpiComm() {
        return mpiCommQuest;
    }

    #if QUEST_COMPILE_SUBCOMM
        void comm_setMpiComm(MPI_Comm newComm) {

            // error if mpiCommQuEST is already set!
            if (mpiCommQuest != MPI_COMM_NULL) {
                MPI_Barrier(mpiCommQuest);
                MPI_Comm_free(&mpiCommQuest);
                error_commDoubleSetMpiComm();
            }

            int mpi_err = MPI_Comm_dup(newComm, &mpiCommQuest);
            if (mpi_err != MPI_SUCCESS) {
                error_commInvalidMpiComm();
            }

            return;
        }
    #endif
#endif
