#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N gfit.homo
#PBS -o gfit.ho.out
#PBS -e gfit.ho.err
#PBS -l walltime=01:20:00
#PBS -l nodes=1:ppn=8 -l mem=16gb
#PBS -M quentinrougemont@orange.fr
#PBS -m ea 
##PBS -t [1-10]%10

#load module for computations
module load compilers/gcc/4.8.5 apps/mugqic_pipeline/2.1.1 mugqic/mugqic_R_packages/0.1_R_3.2.0

# Move to directory where job was submitted
cd "${PBS_O_WORKDIR}"
cd ./results

#Generate 1000 prior for the posterior predictive checks
Rscript ../00-scripts/rscript/03.prior.from.post.R
cd ..

#prepare input for ppc
folder=02.gfit.homo

mkdir "$folder"
cd "$folder"   

#copy inputfile
cp ../bpfile . 
cp ../spinput.txt .
cp ../results/eq1.prior.from.post .
sed -i '/^$/d' eq1.prior.from.post

#global variable
outfile=gfit.outfile
nlines=$(wc -l eq1.prior.from.post |sed -e 's/eq1.prior.from.post//g') 
nreps=1
nsp1=$(grep -v ^$ spinput.txt | sed -n 2~3p | sort -nr  |sed -n 2p )
nsp2=$(grep -v ^$ spinput.txt | sed -n 3~3p | sort -nr  |sed -n 1p )
ntot=$(( $nsp1 + $nsp2 ))
nlocus=$(grep -v ^$ spinput.txt | head -1 )
sed -i 's/myfifo/tmp.ms/g' spinput.txt
sed -i "s/100000/$nlines/g" spinput.txt

#launch msnsam and compute stats
while read line; do 
        echo "$line" > tmp
        col1=$(awk '{print $3}' tmp) #$10=M1
        col2=$(awk '{print $4}' tmp) #$11=M2
        col3=$(awk '{print $1}' tmp)  #$1=N1
        col4=$(awk '{print $2}' tmp)  #$2=N2
        #col5=$(awk '{print $5}' tmp)  #$5=TSC
        #col6=$(awk '{print $4}' tmp)  #$6=Tsplit
        #col7=$(awk '{print $4}' tmp)  #$6=Tsplit
        #col8=$(awk '{print $3}' tmp)  #$3=Nanc
        ../bin/msnsam $ntot $(( $nreps * $nlocus )) -t 0.16 -r 0.08 80 \
        -I 2 $nsp1 $nsp2 0 \
        -m 1 2 $col1 -m 2 1 $col2 -n 1 $col3 -n 2 $col4 > tmp.ms

        ../bin/mscalc spinput.txt

        less ABCstat.txt | grep -v "dataset" >> $outfile
done < sc1.prior.from.post 

wait 

sed -i 's/^$/d' $outfile
#Plot the goodness of fit and compute the p-value
Rscript ../00-scripts/rscript/04.gfit.R 
rm error.txt spoutput.txt tmp tmp.ms seedms bpfile spinput  
#Une fois les goodness of fit ok on peut lancer les enveloppe neutre
#Il faut vérifier à la main les p-valeurs et inspecter les graphes
