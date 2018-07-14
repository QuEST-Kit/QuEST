#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1

# load MPI compiler
module purge
module load mvapich2

# needed on Oxford's ARCUS-B
. enable_arcus-b_mpi.sh

# compile QuEST
make clean
make

# run 8 threads on each of 4 nodes
export OMP_NUM_THREADS=8
mpirun ./myExecutable