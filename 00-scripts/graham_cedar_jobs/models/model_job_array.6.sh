#!/bin/bash
#SBATCH --account=def-blouis #ihv-653-aa
#SBATCH --time=02:45:00
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

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

sims_type=$1 #either 1) "theta" which corresponds to the options -t of ms 
             #2) "single_snps" to generate snp with option -s
             #3) "single_snps_ascertainment" to generate snps with ascertainment bias (for SNPs chips)
			               
if [[ -z "$sims_type" ]]
  then
      echo "Error: need simulation type (eg: theta)"
      exit
fi

# Folder to run simulations
MODEL=./00-scripts/models/"$sims_type"/model.6.sh
FOLDER=./results/am.homom.homon.$SLURM_ARRAY_TASK_ID


for i in $(seq 32)
do
    sleep 0 # $(echo $RANDOM | cut -c -1)
    ./"$MODEL" "$FOLDER"_"$i"  &
done

# Wait for all simulations to finish
wait
sleep 30
