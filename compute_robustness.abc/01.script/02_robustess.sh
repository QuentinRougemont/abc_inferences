#!/bin/bash
# Verify if 100 lines are done, if not run 03_robustness.R

# Global variables
lowbnd=$1
target=$2

# Launch R script
if [[ ! -f 02.done/"$target"/"$lowbnd" ]]
then
	Rscript 01.script/03_robustess.R "$lowbnd" "$target" && touch 02.done/"$target"/"$lowbnd"
fi
