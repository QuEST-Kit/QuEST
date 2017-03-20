#!/bin/bash
# set the number of nodes and processes per node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

##SBATCH --mem=200Gb

# set max wallclock time
#SBATCH --time=00:10:00

# set name of job
#SBATCH --job-name qubitJobRefactor

# set queue
#SBATCH --partition=devel

export OMP_NUM_THREADS=16

module purge

module load intel-compilers

make clean
make

./demo 16
