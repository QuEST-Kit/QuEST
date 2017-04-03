#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2

#SBATCH --mem=50Gb
# uncomment if NUM_QUBITS - log2(NUM_NODES) > 30
####SBATCH --mem=100Gb

# set max wallclock time
#SBATCH --time=00:10:00

# set name of job
#SBATCH --job-name QUEST

# set queue
#SBATCH --partition=compute

NUM_QUBITS=20
EXE=demo
export OMP_NUM_THREADS=8

module purge
module load mvapich2

make clean
make

mpirun $MPI_HOSTS ./$EXE $NUM_QUBITS
