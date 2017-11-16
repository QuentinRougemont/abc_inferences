#!/bin/bash
#PBS -A ihv-653-aa
#PBS -N OutputTest
##PBS -o OutTest.out
##PBS -e OutTest.err
#PBS -l walltime=00:10:00
#PBS -l nodes=1:ppn=8
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
#PBS -t [1-40]%40

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

# Folder to run simulations
MODEL=./models/model.17.sh
FOLDER=./results/2P.bot.sc.heterom.heteron.$MOAB_JOBARRAYINDEX


for i in $(seq 8)
do
    sleep 0 # $(echo $RANDOM | cut -c -1)
    ./"$MODEL" "$FOLDER"_"$i"  &
done

# Wait for all simulations to finish
wait
sleep 30
