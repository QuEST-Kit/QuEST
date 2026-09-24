/** @file
 * 
 * An example of using QuEST's experimental
 * initCustomMpiCommQuESTEnv() function to
 * dedicate only some user-owned MPI processes
 * to QuEST, and dedicate the remainder to
 * other tasks.
 * 
 * @author Oliver Brown
 * @author Tyson Jones (doc)
 */

#include "quest.h"
#include <stdio.h>


// This example requires linking with MPI, which the CMake
// build only enables when QUEST_ENABLE_SUBCOMM is ON, which
// results in quest.h defining QUEST_COMPILE_SUBCOMM. To
// enable this example to always be compilable (like during
// our CI), we guard against when QUEST_ENABLE_SUBCOMM is OFF.
#if ! QUEST_COMPILE_SUBCOMM
int main()
{    
    printf("Example skipped since MPI is not linked.\n");
    return 0;
}
#else 


#include <mpi.h>

int main (void)
{
    int nprocs, quest_nprocs, world_rank, quest_rank;
    MPI_Comm comm_split, comm_quantum, comm_classical;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int I_AM_QUANTUM = world_rank % 2;

    printf("[%d] Hello from rank %d of %d in MPI_COMM_WORLD.\n", world_rank, world_rank, nprocs);

    MPI_Comm_split(MPI_COMM_WORLD, I_AM_QUANTUM, world_rank, &comm_split);

    if (I_AM_QUANTUM) {
        MPI_Comm_dup(comm_split, &comm_quantum);
        MPI_Comm_size(comm_quantum, &quest_nprocs);
        MPI_Comm_rank(comm_quantum, &quest_rank);
        printf("[%d] Hello from rank %d of %d in comm_quantum.\n", world_rank, quest_rank, quest_nprocs);
    } else {
        MPI_Comm_dup(comm_split, &comm_classical);
        quest_rank = -1;
        quest_nprocs = -1;
    }

    // only procs in quantum comm initialise QuEST
    if (I_AM_QUANTUM) {
        printf("[%d] Initialising QuEST.\n", world_rank);
        initCustomMpiCommQuESTEnv(comm_quantum, -1, -1); // -1 = auto-deployments

        reportQuESTEnv();

        printf("[%d] Finalising QuEST.\n", world_rank);
        finalizeQuESTEnv();
    }

    MPI_Comm_free(&comm_split);
    if (I_AM_QUANTUM) {
        MPI_Comm_free(&comm_quantum);
    } else {
        MPI_Comm_free(&comm_classical);
    }

    MPI_Finalize();

    return 0;
}


#endif // QUEST_COMPILE_SUBCOMM
