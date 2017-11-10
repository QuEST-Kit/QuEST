#!/bin/bash

# Single node job using all 16 threads of a 64Gb node
 
# ----------------  EDIT ---------------- 
# select one node 
#SBATCH --nodes=1

# set max wallclock time
#SBATCH --time=01:00:00

# set name of job
#SBATCH --job-name QuEST_JOB

# set the program executable and arguments
NUM_QUBITS=33
EXE=demo
# ----------------------------------------

 
# set up OMP environment
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=16

module purge
# load compiler  
# comment out this line if using gcc 
module load intel-compilers

# compile program. Comment out these lines if program already built 
make clean
make

# run program 
./$EXE $NUM_QUBITS
