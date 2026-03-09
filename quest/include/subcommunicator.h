#ifndef SUBCOMMUNICATOR_H
#define SUBCOMMUNICATOR_H

#include "quest/include/config.h" 

#if COMPILE_MPI && COMPILE_SUBCOMM

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void initCustomMpiCommQuESTEnv(MPI_Comm questComm, int useGpuAccel, int useMultithread);

#ifdef __cplusplus
}
#endif

#endif

#endif
