#!/bin/bash
##PBS -A youracount
##PBS -N neut.fst
##PBS -o neut.fst.homo.out
##PBS -e neut.fst.homo.err
##PBS -l walltime=10:59:00
##PBS -l nodes=1:ppn=8  -l mem=18gb
##PBS -M quentinrougemont@orange.fr
##PBS -m ea 
##PBS -t [1-10]%10
#module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

# Move to directory where job was submitted
#cd "${PBS_O_WORKDIR}"
#cd ./results

#Launch rscript to generate 5e5 prior
Rscript ./00-scripts/rscript/05.eq.generate.prior.from.post.best.homo.R

nloc=1
nsp1=$(grep -v ^$ spinput.txt | sed -n 2p )
nsp2=$(grep -v ^$ spinput.txt | sed -n 3p )
ntot=$(( $nsp1 + $nsp2 ))
folder=04.sim.eq.fst.he.homo

mkdir "$folder"
cd "$folder"    

cp ../results/eq.prior.from.post.homo .
mv eq.prior.from.post.homo prior.from.post 
sed -i '/^$/d' prior.from.post
#dos2unix prior.from.post

echo $nloc >> spinput.txt
for i in $(seq $nloc) ; do echo -e "$nsp1\n$nsp2\n1 " >> spinput.txt ; done
echo -e "$nloc\ntmp.ms" >> spinput.txt

while read line; do 
        echo "$line" > tmp
        col1=$(awk '{print $5}' tmp) #$10=M1
        col2=$(awk '{print $6}' tmp) #$11=M2
        col3=$(awk '{print $1}' tmp)  #$1=N1
        col4=$(awk '{print $2}' tmp)  #$2=N2
        col5=$(awk '{print $4}' tmp)  #$5=Tsplit
        col6=$(awk '{print $3}' tmp)  #$3=Nanc
        ../bin/msnsam $ntot 1 -s 1 \
        -I 2 $nsp1 $nsp2 0 \
        -m 1 2 $col1 -m 2 1 $col2 \
        -n 1 $col3 -n 2 $col4 \
        -ej $col5 2 1 -eN $col5 $col6 > tmp.ms
#        ../bin/mscalc spinput.txt
       
done < prior.from.post 

wait 
grep -A 100 'positions' tmp.ms |sed -e 's/--//g' |grep -v 'positions' > ms.3

Rscript ../00-scripts/rscript/06.neutral.fst.he.R ms.3 $nsp1 $nsp2
wait
sleep 30
rm error.txt spoutput.txt tmp.ms seedms spinput.txt prior.from.post

empstat=../03.fst.he.data/empirical.stats.txt 

cp "$empstat" .
empstat2=empirical.stats.txt
neutstat=neutral.stats.txt
Rscript ../00-scripts/rscript/09.enveloppe.R $neutstat $empstat2  
