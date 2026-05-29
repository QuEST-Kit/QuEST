#include "quest/include/config.h"
#include "quest/include/environment.h"
#include "quest/include/subcommunicator.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"

#if QUEST_COMPILE_MPI && QUEST_COMPILE_SUBCOMM

#include <mpi.h>


// TODO:
// We must resolve this communicator function which contains an MPI type
// and ergo should not be leaked outside comm_config.cpp. For now, we cheat! 
extern void comm_setMpiComm(MPI_Comm newComm);


// TODO:
// We must resolve this inner function of QuEST initialisation, but which is
// private to api/environment.cpp, and so cannot be exposed in the user-facing
// include/environment.hpp. Grr! For now, we here just cheekily extern it c:
extern void validateAndInitCustomQuESTEnv(
    int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread, const char* caller);


void initCustomMpiCommQuESTEnv(MPI_Comm userQuestComm, int useGpuAccel, int useMultithread) {

    // useDistrib and userOwnsMpi are implied by the user of this initialiser
    const int useDistrib = 1;
    const bool userOwnsMpi = true;

    // pre-validate that we are able to set the MPI communicator
    validate_mpiInitStatus(useDistrib, userOwnsMpi, __func__);

    // avoid re-setting the MPI comm (to avoid an internal error), which happens
    // if a user illegally re-calls this function, which will be subsequently
    // caught by the validation in validateAndInitCustomQuESTEnv() below
    if (comm_isMpiCommSet())
        comm_setMpiComm(userQuestComm);

    // perform remaining validation (some is harmlessly repeated) and init QuEST env
    validateAndInitCustomQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread, __func__);
}

#endif
