#!/bin/bash

# Single node using 8 processes, each with 16 threads, on the 6Tb node

# ----------------  EDIT ---------------- 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# set max wallclock time
#SBATCH --time=01:00:00

# set name of job
#SBATCH --job-name QuEST_JOB
# ----------------------------------------


# select 6TB node
#SBATCH --partition=mem6T

# run one job at a time on the 6TB node
#SBATCH --exclusive


# ----------------  EDIT ---------------- 
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

