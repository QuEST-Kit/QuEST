#include "quest/include/config.h"
#include "quest/include/environment.h"
#include "quest/include/subcommunicator.h"

#include "quest/src/comm/comm_config.hpp"
#include "quest/src/core/errors.hpp"

#if QUEST_COMPILE_MPI && QUEST_COMPILE_SUBCOMM

#include <stdbool.h>
#include <mpi.h>

void initCustomMpiCommQuESTEnv(MPI_Comm userQuestComm, int useGpuAccel, int useMultithread) {
    // useDistrib and userOwnsMpi are implied by the user of this initialiser
    const int useDistrib = 1;
    const bool userOwnsMpi = true;

    // set mpiCommQuest to user provided communicator
    if (comm_isInit()) {
        comm_setMpiComm(userQuestComm);
    } else {
        error_commNotInit();
    }

    // initialise QuEST around that communicator
    initCustomMpiQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread);

    return;
}

#endif
