#ifndef SUBCOMMUNICATOR_H
#define SUBCOMMUNICATOR_H

#include "quest/include/config.h" 

#if QUEST_COMPILE_MPI && QUEST_COMPILE_SUBCOMM

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @notyetdoced
 *  Advanced initialiser which allows the user to provide an MPI communicator for QuEST to use.
 *  Use of this initialiser implies userOwnsMpi = true, (exposed by initCustomMpiQuESTEnv) and 
 *  therefore that they have already initialised MPI, and they will call MPI_Finalize at the 
 *  appropriate time.
 *
 *  The user-provided MPI communicator undergoes the same validation procedure as any that QuEST
 *  would use, and so must contain a power-of-2 number of processes.
 */
void initCustomMpiCommQuESTEnv(MPI_Comm questComm, int useGpuAccel, int useMultithread);

#ifdef __cplusplus
}
#endif

#endif

#endif
