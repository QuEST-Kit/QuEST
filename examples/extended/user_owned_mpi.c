/** @file
 * 
 * TODO
 * 
 * @author Oliver Brown
 */

#include <mpi.h>
#include "quest.h"


    // TODO:
    // this file will only receive mpi.h from CMakeLists.txt if
    // we are also compiling with QUEST_ENABLE_SUBCOMM. Fix this!


int main (void)
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
