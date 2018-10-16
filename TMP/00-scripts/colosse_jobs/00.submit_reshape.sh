#!/bin/bash
#Header for colosse calcul quebec
#PBS -A ihv-653-ab
#PBS -N reshape
##PBS -o reshape.out
##PBS -e reshape.err
#PBS -l walltime=00:25:00
#PBS -l nodes=1:ppn=8
##PBS -M quentinrougemont@orange.fr
##PBS -m bea

# Move to directory where job was submitted
cd $PBS_O_WORKDIR

./00-scripts/00.reshape.sh
