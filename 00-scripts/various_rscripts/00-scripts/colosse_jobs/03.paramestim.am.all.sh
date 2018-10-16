#!/00-scripts/bash
#PBS -A ihv-653-ab
#PBS -N param.am.het.het
#PBS -o ParamAM.het.het.out
#PBS -e ParamAM.het.het.err
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=8  -l mem=23gb
#PBS -M quentinrougemont@orange.fr
#PBS -m ea 
##PBS -t [1-10]%10

# Move to directory where job was submitted
cd PBS_O_WORKDIR

module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

cd ./results
Rscript ../00-scripts/rscript/03.ParamEstimAM.all.R

wait
sleep 30
