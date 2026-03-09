#include "quest/include/config.h"
#include "quest/include/environment.h"
#include "quest/include/subcommunicator.h"

#include "quest/src/comm/comm_config.hpp"

#if COMPILE_MPI && COMPILE_SUBCOMM

#include <mpi.h>

void initCustomMpiCommQuESTEnv(MPI_Comm userQuestComm, int useGpuAccel, int useMultithread) {
    // useDistrib and userOwnsMpi are implied by the user of this initialiser
    const int USE_DISTRIB = 1;
    const int USER_OWNS_MPI = 1;

    // set mpiCommQuest to user provided communicator
    comm_setMpiComm(userQuestComm);

    // initialise QuEST around that communicator
    initCustomMpiQuESTEnv(USE_DISTRIB, USER_OWNS_MPI, useGpuAccel, useMultithread);

    return;
}

#endif
