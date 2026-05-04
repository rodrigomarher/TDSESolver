#!/bin/bash

#SBATCH --job-name="TDSE"
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --nodelist=nodo11
##SBATCH --account=usal74
##SBATCH --qos=gp_resa

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
srun python3 launch.py
