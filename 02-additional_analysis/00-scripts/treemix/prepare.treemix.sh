#!/bin/bash

#Script pour construire ped files et treemix files
path=02-additional_analysis/00-scripts/treemix #your path to the script
datafolder=02-additional_analysis/02-treemix #path to folder

awk '{print $1}' "$datafolder"/cluster.dat > "$datafolder"/colomn1
cut -d " " -f2-20000 "$datafolder"/salmon.ordered.V2.full.ped > "$datafolder"/colomn2.4656

paste "$datafolder"/colomn1 "$datafolder"/colomn2.4656 > "$datafolder"/salmon.ordered.V2.tmp 

sed -e 's/N A/0 0/g' "$datafolder"/salmon.ordered.V2.tmp  -e "s/\ \ */\ /g" -e 's/\t/ /g' > "$datafolder"/salmon.ordered.V2.ped

#use plink1.9
plink --file "$datafolder"/salmon.ordered.V2 --noweb --missing --freq --double-id --allow-extra-chr --chr-set 29 --within "$datafolder"/cluster.dat --out "$datafolder"/plink

rm "$datafolder"/*tmp

gzip "$datafolder"/plink.frq.strat
$path/plink2treemix.py "$datafolder"/plink.frq.strat.gz "$datafolder"/treemix.frq.gz

#faire le f3test:
threepop -i "$datafolder"/treemix.frq.gz -k 100 >> "$datafolder"/resultat.f3

#analyser le results de threepop
mkdir "$datafolder"/f3.test
mv "$datafolder"/resultat.f3 "$datafolder"/f3.test/
cd "$datafolder"/f3.test/
../../00-scripts/treemix/analyse.threepop.sh resultat.f3
