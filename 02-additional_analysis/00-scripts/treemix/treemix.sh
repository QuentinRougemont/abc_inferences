#!/bin/bash

#SBATCH -D ./
#SBATCH -J treemix
#SBATCH -o treemix.%j.out
#SBATCH -c 1
#SBATCH -p ibismax
#SBATCH -A ibismax
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR@MAIL
#SBATCH --time=1-00:00
#SBATCH --mem=10G

# Load the software with module if applicable:
module load treemix
# Type your command line here

#cd $SLURM_SUBMIT_DIR

#Global variable
i="-i treemix.frq.gz" #name of input file
#path="/usr/local/bioinfo/src/treemix/treemix-1.12/bin/" #path
m="-m" #migration
g="-g out_stem.vertices.gz out_stem.edges.gz "
o="-o out_stem_mig"
#b="-bootstrap"
k="-k 100"

#run treemix
treemix $i  -o out_stem

for K in $(seq 20)  ; do treemix $i  $m $K $b $k $g $o.$K ; done 
