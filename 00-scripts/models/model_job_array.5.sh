#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N OutputTest
##PBS -o OutTest.out
##PBS -e OutTest.err
#PBS -l walltime=09:25:00
#PBS -l nodes=1:ppn=8
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
#PBS -t [1-10]%10

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

# Folder to run simulations
MODEL=./00-scripts/models/model.5.sh
FOLDER=./results/am.homom.heteron.$MOAB_JOBARRAYINDEX


for i in $(seq 8)
do
    sleep 0 # $(echo $RANDOM | cut -c -1)
    ./"$MODEL" "$FOLDER"_"$i"  &
done

# Wait for all simulations to finish
wait
sleep 30
