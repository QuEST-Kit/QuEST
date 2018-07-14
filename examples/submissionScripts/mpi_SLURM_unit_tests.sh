#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --output=results.txt

# load MPI compiler
module purge
module load mvapich2

# needed on Oxford's ARCUS-B
. enable_arcus-b_mpi.sh

# requires ../../makefile has DISTRIBUTED=1
cd ../../tests
./runTests.sh