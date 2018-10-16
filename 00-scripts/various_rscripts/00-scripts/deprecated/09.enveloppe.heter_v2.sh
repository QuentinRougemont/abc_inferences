#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N neut.fst
#PBS -o neut.fst.hetero_v2.out
#PBS -e neut.fst.hetero_v2.err
#PBS -l walltime=00:59:00
#PBS -l nodes=1:ppn=8  -l mem=18gb
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
##PBS -t [1-10]%10
module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"

#Launch rscript to generate 5e5 prior
folder=05.sim.fst.he.hetero_v2

cd "$folder"    

empstat=../03.fst.he.data/empirical.stats.txt 

cp "$empstat" .
empstat2=empirical.stats.txt
neutstat=neutral.stats.txt

Rscript ../00-scripts/rscript/09.enveloppe_hetero.R $neutstat $empstat2 
