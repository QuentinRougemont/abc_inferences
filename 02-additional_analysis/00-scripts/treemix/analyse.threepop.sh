#!/bin/bash 

INPUT=$1
path=../../00-scripts/treemix

grep -v Estimating "$INPUT" |grep -v total |grep -v npop |sed -e 's/;/ /g' |sed -e 's/,/ /g' | LC_ALL=C sort -k 6 -n > f3.100.sorted

numline=$(awk '{print $6}' f3.100.sorted |grep -Eo "\-[0-9]+(\.[0-9]+)?" |wc -l)

head -n $numline f3.100.sorted > f3.100.sorted.significantly.admixed

sort f3.100.sorted.significantly.admixed  |awk '{print $1}' |uniq >> list.pop.significantly.admixed
sort f3.100.sorted.significantly.admixed  |awk '{print $2}' |sort |uniq >> major.source.2
sort f3.100.sorted.significantly.admixed  |awk '{print $3}' |sort |uniq >> major.source.3

awk '{print $1 }' f3.100.sorted.significantly.admixed |sort |uniq -c > f3.sorted.significantly.admixed.count
#sort f3.100.sorted.significantly.admixed  |awk '{print $1}' |grep lab |wc -l >> count.lab.admixed

for i in $(cat list.pop.significantly.admixed) 
do 
    sort f3.100.sorted.significantly.admixed |awk '{print $1}' |grep $i |wc -l >> count.$i.admixed ; 
done

for i in $(cat major.source.2)
do
    sort f3.100.sorted.significantly.admixed |awk '{print $2}' |grep $i |wc -l >> count.$i.source
done

for i in $(cat major.source.3)
do
    sort f3.100.sorted.significantly.admixed |awk '{print $3}' |grep $i |wc -l >> count.$i.source.3
done

for i in $(cat major.source.2)
do
    sed -i "1i $i" count.$i.source
done
for i in $(cat major.source.3)
do 
    sed -i "1i $i" count.$i.source.3
done

for i in $(cat list.pop.significantly.admixed)
do 
    sed -i "1i $i" count.$i.admixed
done

paste count.*.source > source.major.1
paste count.*.source.3 > source.major.2
paste count.*.admixed > admixed.major

$path/transpose_tab source.major.1 | sort -k 2 -n -r  >  source.major.1.tr
$path/transpose_tab source.major.2 | sort -k 2 -n -r  > source.major.2.tr
$path/transpose_tab admixed.major  | sort -k 2 -n -r  > admixed.major.tr


awk '{print $1}' source.major.2.tr > maj.2.tmp
awk '{print $1}' source.major.1.tr > maj.1.tmp
 
cat maj.1.tmp maj.2.tmp |sort |uniq > maj.all.tmp

mkdir tmp
mv count* tmp
mv admixed.major tmp/
mv source.major.1 source.major.2 tmp/
mv major.source.2 major.source.3 tmp/

#rm -r tmp
