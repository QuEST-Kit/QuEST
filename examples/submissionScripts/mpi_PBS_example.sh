#!/bin/bash

#PBS -l select=4:ncpus=8

# compile QuEST into myExecutable
make clean
make

# run 8 threads on each of 4 nodes
export OMP_NUM_THREADS=8
aprun -n 4 -d 8 -cc numa_node ./myExecutable