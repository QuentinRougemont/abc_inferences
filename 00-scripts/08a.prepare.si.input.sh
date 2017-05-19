#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N neutral.fst
#PBS -o neutral.fst.estim.out
#PBS -e neutral.fst.esim.err
#PBS -l walltime=00:09:00
#PBS -l nodes=1:ppn=8  -l mem=8gb
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
##PBS -t [1-10]%10
cd "${PBS_O_WORKDIR}"

threshold=95
folder=11.isol.selection."$threshold"
#data=../*ABC.tmp

mkdir "$folder"
cd "$folder"
#cp "$data" .
cp ../*ABC.tmp ABC 

sed -e '1,2d' ABC |cut -f2-20000 | sed -e 's/NA/--/g' -r -e 's/([^ ])/\1 /g'  -e 's/\t/ /g' > geno.tmp

sed -n '1,2p' ABC > ABC.tmp

module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

#Launch Rscript 
Rscript ../00-scripts/rscript/10.prepare.data.si.R ABC.tmp $threshold

input="*."$threshold".txt"
i=$( echo $input |sed -e "s/."$threshold".txt//g" )
z=$(echo "$i."$threshold".txt ${i:0:$(expr ${#i} / 2)} ${i:$(expr ${#i} / 2)} " )

../bin/inputGen_modifquentin.py $z

cp spinput.txt spinput.txt.2 
sed -i -e 's/myfifo/locus.ms/g' -e 's/100000/1/g' spinput.txt 

../bin/mscalc spinput.txt 

mv ABCstat.txt OBS.ABC.stat.txt 
#renommer les spinput.txt2 en spinput.txt
mv spinput.txt.2 spinput.txt 
#vire les spoutput, error.txt, etc
rm error.txt  spoutput.txt  
  
