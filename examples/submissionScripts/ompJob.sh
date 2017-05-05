#!/bin/bash

# set the number of nodes and processes per node. We are running one process on a single node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# Choose the right amount of memory to ask for
#SBATCH --mem=50Gb

# uncomment if NUM_QUBITS > 31
####SBATCH --mem=100Gb
# uncomment if NUM_QUBITS > 32
####SBATCH --mem=200Gb

# set max wallclock time
#SBATCH --time=00:10:00

# set name of job
#SBATCH --job-name QUEST

# set queue
#SBATCH --partition=compute

NUM_QUBITS=16
EXE=demo
export OMP_NUM_THREADS=16

module purge
module load intel-compilers/2016

make clean
make

./$EXE $NUM_QUBITS
