#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# set max wallclock time
#SBATCH --time=01:30:00

# set name of job
#SBATCH --job-name QUEST_AB

# set queue
#SBATCH --partition=mem6T
#SBATCH --exclusive

NUM_QUBITS=33
EXE=demo
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=128

module purge
module load intel-compilers

make clean
make

./$EXE $NUM_QUBITS
