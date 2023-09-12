#!/bin/bash
#SBATCH -J mlg1.8         # job name
#SBATCH -o myMPI.o%j       # output and error file name (%j expands to jobID)
#SBATCH -p g6
#SBATCH -w n019,n020,n022,n023,n025
#SBATCH --ntasks=4            # total number of nodes
#SBATCH --cpus-per-task=16             # total number of tasks
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## HPC ENVIRONMENT 
. /etc/profile.d/TMI.sh
##

mpirun ./pes > stdout.log
