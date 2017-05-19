#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N neut.fst
#PBS -o neut.fst.hetero_v2.out
#PBS -e neut.fst.hetero_v2.err
#PBS -l walltime=11:59:00
#PBS -l nodes=1:ppn=8  -l mem=18gb
#PBS -M quentinrougemont@orange.fr
##PBS -m ea 
##PBS -t [1-10]%10
module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"
cd ./results

#Launch rscript to generate 5e5 prior
Rscript ../00-scripts/rscript/05.generate.prior.from.post.best.het.het.R
cd ..

nloc=1
nsp1=$(grep -v ^$ spinput.txt | sed -n 2~3p | sort -nr  |sed -n 2p )
nsp2=$(grep -v ^$ spinput.txt | sed -n 3~3p | sort -nr  |sed -n 1p )
#nsp1=$(grep -v ^$ spinput.txt | sed -n 2p )
#nsp2=$(grep -v ^$ spinput.txt | sed -n 3p )
ntot=$(( $nsp1 + $nsp2 ))
folder=05.sim.fst.he.hetero_v2

mkdir "$folder"
cd "$folder"    

cp ../results/prior.from.post.het.het_v2 .
mv prior.from.post.het.het_v2 prior.from.post 
sed -i '/^$/d' prior.from.post

echo $nloc > spinput.txt
for i in $(seq $nloc) ; do echo -e "$nsp1\n$nsp2\n1 " >> spinput.txt ; done
echo -e "$nloc\ntmp.ms" >> spinput.txt

while read line; do 
        echo "$line" > tmp
        col1=$(awk '{print $10}' tmp) #$10=M1
        col2=$(awk '{print $11}' tmp) #$11=M2
        col3=$(awk '{print $1}' tmp)  #$1=N1
        col4=$(awk '{print $2}' tmp)  #$2=N2
        col5=$(awk '{print $5}' tmp)  #$5=TSC
        col6=$(awk '{print $4}' tmp)  #$6=Tsplit
        col7=$(awk '{print $4}' tmp)  #$6=Tsplit
        col8=$(awk '{print $3}' tmp)  #$3=Nanc
        ../bin/msnsam $ntot 1 -s 1 -I 2 $nsp1 $nsp2 0 -m 1 2 $col1 -m 2 1 $col2 -n 1 $col3 -n 2 $col4 -eM $col5 0 -ej $col6 2 1 -eN $col7 $col8 >> tmp.ms
        
done < prior.from.post 

wait 

grep -A $ntot 'positions' tmp.ms |sed -e 's/--//g' |grep -v 'positions' > ms.3

Rscript ../00-scripts/rscript/06.neutral.fst.he.R ms.3 $nsp1 $nsp2

wait
sleep 30

rm error.txt spoutput.txt seedms spinput.txt prior.from.post

empstat=../03.fst.he.data/empirical.stats.txt 

cp "$empstat" .
empstat2=empirical.stats.txt
neutstat=neutral.stats.txt

Rscript ../00-scripts/rscript/09.enveloppe_hetero.R $neutstat $empstat2  
