/** @file
 * 
 * An example of using QuEST's experimental
 * initCustomMpiQuESTEnv() function to
 * initialise QuEST in an environment where
 * MPI is owned and controlled by the user.
 * 
 * @author Oliver Brown
 * @author Tyson Jones (doc)
 */

#include "quest.h"
#include <cstdio>


// This example requires linking with MPI, which the CMake
// build only enables when QUEST_ENABLE_SUBCOMM is ON, which
// results in quest.h defining QUEST_COMPILE_SUBCOMM. To
// enable this example to always be compilable (like during
// our CI), we guard against when QUEST_ENABLE_SUBCOMM is OFF.
#if ! QUEST_COMPILE_SUBCOMM
int main(void)
{    
    std::printf("Example skipped since MPI is not linked.\n");
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
