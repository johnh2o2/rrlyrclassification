#!/bin/bash
# parallel job using 48 cores. and runs for 4 hours (max)
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=12
#SBATCH -t 00:05:00
# sends mail when process begins, and 
# when it ends. Make sure you define your email 
# address.
#XSBATCH --mail-type=begin
#XSBATCH --mail-type=end
#XSBATCH --mail-user=jah5@princeton.edu
module load openmpi
module load anaconda
module load mpi4py

srun python test.py
