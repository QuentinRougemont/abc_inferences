#!/bin/bash

#script to reshape the data
#to be use carefully
target=1000000 #set the number of wanted simulations

ls -d results/*/ | sed -e 's/\([0-9]*\)//g' -e 's/._\///g'  -e 's/results\///g' |uniq > list

for j in $(cat list) ; do mkdir results/$j.glob ; done 
for j in  $(cat list) ; do mv results/$j.* results/$j.glob ; done 
for i in $(cat list) ; do for k in $(find results/$i.glob -name ABCstat.txt) ; do cat $k |grep -v dataset >> results/$i.ABC.stat.txt ; done ; done  
for i in $(cat list) ; do for k in $(find results/$i.glob -name priorfile) ; do cat $k |grep -v N1 >> results/$i.priorfile.txt ; done ; done 
sed -i '/^$/d' results/*.h*.ABC.stat.txt 

for j in  *.glob ; do rm results/$j/*/myfifo  results/$j/*/error.txt results/$j/*/sp*txt results/$j/*/bpfile results/$j/*/seedms ; done
for i in results ; do for j  in  $i/*.glob  ; do  tar -zcvf $j.tar.gz $j ; done ;done
for j in $(ls -d results/*glob/ ) ; do  rm  -r $j ; done
#rm *.err *.out

for i in $(cat list) ; do wc -l results/$i.ABC.stat.txt  |awk '{print $1}' >>  stat.tmp ; done
for i in $(cat list) ; do wc -l results/$i.priorfile.txt |awk '{print $1}' >> prior.tmp ; done 

mv OBS.ABC.stat.txt results/

sum_stat=$(awk '{total += $1} END {print total/NR}' stat.tmp )
if [[ $sum_stat != $target ]]  ; 
then 
	echo "you loose" 
	echo "error! Number of simulation in one of the model is not = "$target" !!"
	echo "please check models"
#	exit
fi

sum_prio=$(awk '{total += $1} END {print total/NR}' prior.tmp )
if [[ $sum_prio != $target ]]  ; 
then 
	echo "you loose" 
	echo "error! Number of simulation in one of the model is not = "$target" !!"
	echo "please check models"
#	exit
fi

rm list *.tmp 
