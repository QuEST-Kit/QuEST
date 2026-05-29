/** @file
 * 
 * TODO
 * 
 * @author Oliver Brown
 */

#include <cstdio>
#include <mpi.h>
#include <quest.h>


    // TODO:
    // this file will only receive mpi.h from CMakeLists.txt if
    // we are also compiling with QUEST_ENABLE_SUBCOMM. Fix this!


int main (void)
{
    int nprocs, quest_nprocs, world_rank, quest_rank;
    MPI_Comm comm_split, comm_quantum, comm_classical;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int I_AM_QUANTUM = world_rank % 2;

    std::printf("[%d] Hello from rank %d of %d in MPI_COMM_WORLD.\n", world_rank, world_rank, nprocs);

    MPI_Comm_split(MPI_COMM_WORLD, I_AM_QUANTUM, world_rank, &comm_split);

    if (I_AM_QUANTUM) {
        MPI_Comm_dup(comm_split, &comm_quantum);
        MPI_Comm_size(comm_quantum, &quest_nprocs);
        MPI_Comm_rank(comm_quantum, &quest_rank);
        std::printf("[%d] Hello from rank %d of %d in comm_quantum.\n", world_rank, quest_rank, quest_nprocs);
    } else {
        MPI_Comm_dup(comm_split, &comm_classical);
        quest_rank = -1;
        quest_nprocs = -1;
    }

    // only procs in quantum comm initialise QuEST
    if (I_AM_QUANTUM) {
        std::printf("[%d] Initialising QuEST.\n", world_rank);
        initCustomMpiCommQuESTEnv(comm_quantum, modeflag::USE_AUTO, modeflag::USE_AUTO);

        reportQuESTEnv();

        std::printf("[%d] Finalising QuEST.\n", world_rank);
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
