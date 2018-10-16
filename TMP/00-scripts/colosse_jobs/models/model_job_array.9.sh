#!/bin/bash
#PBS -A ihv-653-aa
#PBS -N OutputTest
##PBS -o OutTest.out
##PBS -e OutTest.err
#PBS -l walltime=38:00:00
#PBS -l nodes=1:ppn=8
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
#PBS -t [1-40]%40

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

sims_type="theta" #either 1) "theta" which corresponds to the options -t of ms 
             #2) "single_snps" to generate snp with option -s
             #3) "single_snps_ascertainment" to generate snps with ascertainment bias (for SNPs chips)
			               
if [[ -z "$sims_type" ]]
  then
      echo "Error: need simulation type (eg: theta)"
      exit
fi

# Folder to run simulations
MODEL=./00-scripts/models/"$sims_type"/model.9.sh
FOLDER=./results/sc.homom.heteron.$MOAB_JOBARRAYINDEX


for i in $(seq 8)
do
    sleep 0 # $(echo $RANDOM | cut -c -1)
    ./"$MODEL" "$FOLDER"_"$i"  &
done

# Wait for all simulations to finish
wait
sleep 30
