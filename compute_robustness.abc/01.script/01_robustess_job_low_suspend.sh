#!/bin/bash
#SBATCH -J "ABCRobust"
#SBATCH -o log_%j
#SBATCH -c 10
#SBATCH -p low-suspend
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=asdf
#SBATCH --time=10-00:00
#SBATCH --mem=200G

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

# Get parameters
target=$1
if [[ -z "$target" ]]
then
    echo "Error: need model name (eg: AM.1)"
    exit
fi

# Launching jobs in parallel
echo running model: "$target"
cat lower_bounds | parallel -j 10 01.script/02_robustess.sh {} "$target"
