#!/bin/bash
#PBS -A ihv-653-aa
#PBS -N OutputTest
##PBS -o OutTest.out
##PBS -e OutTest.err
#PBS -l walltime=02:45:00
#PBS -l nodes=1:ppn=8
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
#PBS -t [1-10]%10

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

# Folder to run simulations
MODEL=./models/theta/model.0.sh
FOLDER=./results/pan.homom.heteron.$SLURM_ARRAY_TASK_ID

for i in $(seq 32)
do
    sleep 0 # $(echo $RANDOM | cut -c -1)
    ./"$MODEL" "$FOLDER"_"$i"  &
done

# Wait for all simulations to finish
wait
sleep 30
