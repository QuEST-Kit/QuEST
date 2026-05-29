/** @file
 * 
 * TODO
 * 
 * @author Oliver Brown
 */

#include "quest.h"
#include <stdio.h>

// This example requires linking with MPI, which the CMake
// build only enables when QUEST_ENABLE_SUBCOMM is ON, which
// results in quest.h defining QUEST_COMPILE_SUBCOMM
#if ! QUEST_COMPILE_SUBCOMM

int main(void)
{    
    printf("Example skipped since MPI is not linked.\n");
    return 0;
}

#else 

#include <mpi.h>

int main(void)
{
    const int  USE_DISTRIB = 1;
    const bool USER_MPI    = 1;
    const int  USE_OPENMP  = 1;
    const int  USE_GPU     = 0;

    MPI_Init(NULL, NULL);
    initCustomMpiQuESTEnv(USE_DISTRIB, USER_MPI, USE_GPU, USE_OPENMP);
    reportQuESTEnv();
    finalizeQuESTEnv();
    MPI_Finalize();

    return 0;
}

#endif // QUEST_COMPILE_SUBCOMM
