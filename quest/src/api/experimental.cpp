/** @file
 * Experimental functions which are liable to
 * API breaks within QuEST minor version releases.
 * Some optional functions require compiling this
 * file against MPI, despite being outside of /comm/, 
 * and so require opt-in macros (QUEST_COMPILE_SUBCOMM)
 * 
 * @author Oliver Brown
 */

#include "quest/include/config.h"
#include "quest/include/environment.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"

#if QUEST_COMPILE_SUBCOMM && ! QUEST_COMPILE_MPI
    #error "Macro QUEST_COMPILE_SUBCOMM was true, but QUEST_COMPILE_MPI was illegally false."
#endif

#if QUEST_COMPILE_SUBCOMM
    #include <mpi.h>
#endif



/*
 * EXTERNAL FUNCTIONS
 *
 * which we here regretfully 'extern' because we are either
 * unsure which header should expose them, or because they
 * contain deployment-specific types (like MPI_Comm) which
 * we do not wish to expose within internal headers 
 */


extern void validateAndInitCustomQuESTEnv(
    int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread, const char* caller);


#if QUEST_COMPILE_SUBCOMM // hide MPI_Comm
    extern bool comm_setMpiComm(MPI_Comm newComm, bool userOwnsMpi);
#endif



/*
 * API FUNCTIONS
 */


// enable invocation by both C and C++ binaries
extern "C" {


void initCustomMpiQuESTEnv(int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread) {
    validateAndInitCustomQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread, __func__);
}


#if QUEST_COMPILE_SUBCOMM // hide MPI_Comm

void initCustomMpiCommQuESTEnv(MPI_Comm userQuestComm, int useGpuAccel, int useMultithread) {

    // useDistrib and userOwnsMpi are implied by the user of this initialiser
    const int useDistrib = 1;
    const bool userOwnsMpi = true;

    // pre-validate that we are able to set the MPI communicator
    validate_mpiInitStatus(useDistrib, userOwnsMpi, __func__);
    validate_mpiSubCommIsNonNull(userQuestComm != MPI_COMM_NULL, __func__);

    // avoid re-setting the MPI comm (to avoid an internal error), which happens
    // if a user illegally re-calls this function, which will be subsequently
    // caught by the validation in validateAndInitCustomQuESTEnv() below
    if (!comm_isActive()) {
        bool success = comm_setMpiComm(userQuestComm, userOwnsMpi);
        validate_mpiSubCommSetSucceeded(success, __func__);
    }

    // perform remaining validation (some is harmlessly repeated) and init QuEST env
    validateAndInitCustomQuESTEnv(useDistrib, userOwnsMpi, useGpuAccel, useMultithread, __func__);
}

#endif // QUEST_COMPILE_SUBCOMM


// end de-mangler
}
