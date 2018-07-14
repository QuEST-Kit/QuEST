#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1 

#SBATCH --partition=gpu    ## name may vary

# load CUDA compilers
module purge
module load cuda  ## name may vary

# compile QuEST into myExecutable
make clean
make

# run on GPU
./myExecutable