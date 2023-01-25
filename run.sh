#!/bin/bash

#SBATCH --job-name="TDSE"
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
srun python3 launch.py
