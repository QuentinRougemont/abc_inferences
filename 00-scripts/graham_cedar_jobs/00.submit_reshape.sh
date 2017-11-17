#!/bin/bash
#SBATCH --account=def-blouis #ihv-653-aa
#SBATCH --time=01:35:00
#SBATCH --job-name=abc
#SBATCH --output=abc-%J.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

./00-scripts/00.reshape.sh
