#!/bin/bash

input=$1 #genotypic matrice obtained from the script to prepare ABCdata

sed '1,2'd "$input" | cut -f 2- > geno
sed -n '1,2p' geno > header
sed 's/\(.\)/\1 /g;s/ $//' geno | sed 's/ \+ /\t/g' > geno.tmp
