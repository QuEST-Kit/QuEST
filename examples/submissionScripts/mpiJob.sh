#!/bin/bash

# Multiple node job using - one process with 16 threads - per node

# ----------------  EDIT ---------------- 
# select one node 
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1

# set max wallclock time
#SBATCH --time=01:00:00

# set name of job
#SBATCH --job-name QuEST_JOB

# set the program executable and arguments
NUM_QUBITS=33
EXE=demo
# ----------------------------------------


# set up OMP environment
export OMP_NUM_THREADS=16

module purge
# load compiler  
module load mvapich2

# set up mpi on the arcus system
. enable_arcus-b_mpi.sh

# compile program. Comment out these lines if program already built 
make clean
make

# run program 
mpirun $MPI_HOSTS ./$EXE $NUM_QUBITS


