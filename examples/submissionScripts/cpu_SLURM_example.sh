#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# load C/C++ compiler
module purge
module load gcc

# compile QuEST into myExecutable
make clean
make

# run with 8 threads
export OMP_NUM_THREADS=8
./myExecutable