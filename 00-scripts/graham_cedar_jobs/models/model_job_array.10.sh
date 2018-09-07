#!/bin/bash
#SBATCH --account=def-blouis #ihv-653-aa
#SBATCH --time=03:45:00
#SBATCH --job-name=abc
#SBATCH --output=abc-%J.out
#SBATCH --array=1-33%33
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
##SBATCH --gres=cpu:32
##SBATCH --mail-user=quentinrougemont@orange.fr
##PBS -l nodes=1:ppn=8
##SBATCH --mail-type=EA 
source /home/quentin/software/abcinf/temp/bin/activate
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Folder to run simulations
MODEL=./00-scripts/models/theta/model.10.sh
FOLDER=./results/im.heterom.heteron.$SLURM_ARRAY_TASK_ID


for i in $(seq 32)
do
    sleep 0 # $(echo $RANDOM | cut -c -1)
    ./"$MODEL" "$FOLDER"_"$i"  &
done

# Wait for all simulations to finish
wait
sleep 30
