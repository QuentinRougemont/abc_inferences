#!/bin/bash
# WARNING! All scripts must been edited prior to submitting and not changed until run!
# HINT: Give meaningful job names in submission scripts to distinguish them

# Global variables
SCRIPTPATH="./00-scripts"

echo "    DO NOT USE WITHOUT APPROPRIATE ABC BACKGROUND"

# Launch the models :
p01=$(msub  "$SCRIPTPATH"/models/model_job_array.1.sh  | tail -n 1 |awk  '{print $1}')  #si n.hetero
p02=$(msub  "$SCRIPTPATH"/models/model_job_array.2.sh  | tail -n 1 |awk  '{print $1}')  #si n.homo
p03=$(msub  "$SCRIPTPATH"/models/model_job_array.3.sh  | tail -n 1 |awk  '{print $1}')  #am n.hetero.m.hetero
p04=$(msub  "$SCRIPTPATH"/models/model_job_array.4.sh  | tail -n 1 |awk  '{print $1}')  #am n.homo.m.hetero
p05=$(msub  "$SCRIPTPATH"/models/model_job_array.5.sh  | tail -n 1 |awk  '{print $1}')  #am n.hetero.m.homo
p06=$(msub  "$SCRIPTPATH"/models/model_job_array.6.sh  | tail -n 1 |awk  '{print $1}')  #am n.homo.m.homo
p07=$(msub  "$SCRIPTPATH"/models/model_job_array.7.sh  | tail -n 1 |awk  '{print $1}')  #sc.heterom.heteron
p08=$(msub  "$SCRIPTPATH"/models/model_job_array.8.sh  | tail -n 1 |awk  '{print $1}')  #sc.heterom.homon
p09=$(msub  "$SCRIPTPATH"/models/model_job_array.9.sh  | tail -n 1 |awk  '{print $1}')  #sc.homom.heteron
p10=$(msub  "$SCRIPTPATH"/models/model_job_array.10.sh  | tail -n 1 |awk  '{print $1}') #im.heterom.heteron
p11=$(msub  "$SCRIPTPATH"/models/model_job_array.11.sh  | tail -n 1 |awk  '{print $1}') #im.heterom.homon
p12=$(msub  "$SCRIPTPATH"/models/model_job_array.12.sh  | tail -n 1 |awk  '{print $1}') #im.homom.heteron
p13=$(msub  "$SCRIPTPATH"/models/model_job_array.13.sh  | tail -n 1 |awk  '{print $1}') #sc.homom.homon
p14=$(msub  "$SCRIPTPATH"/models/model_job_array.14.sh  | tail -n 1 |awk  '{print $1}') #im.homom.homon
#p15=$(msub  "$SCRIPTPATH"/models/model_job_array.15.sh  | tail -n 1 |awk  '{print $1}') #bot.sc.heterom.heteron
#p16=$(msub  "$SCRIPTPATH"/models/model_job_array.16.sh  | tail -n 1 |awk  '{print $1}') #2P.sc.heterom.heteron
#p17=$(msub  "$SCRIPTPATH"/models/model_job_array.17.sh  | tail -n 1 |awk  '{print $1}') #2P.bot.sc.heterom.heteron
#p18=$(msub  "$SCRIPTPATH"/models/model_job_array.18.sh  | tail -n 1 |awk  '{print $1}') #pan heteron
#p19=$(msub  "$SCRIPTPATH"/models/model_job_array.19.sh  | tail -n 1 |awk  '{print $1}') #pan homon
#p20=$(msub  "$SCRIPTPATH"/models/model_job_array.20.sh  | tail -n 1 |awk  '{print $1}') #bot heteron
#p21=$(msub  "$SCRIPTPATH"/models/model_job_array.21.sh  | tail -n 1 |awk  '{print $1}') #bot homom

# reshape the simulation and start the abc model choice procedure
#Warning you have to chose after which model to start this script! Depending on the compared model
s01=$(msub -l depend=$p01:$p02:$p03:$p04:$p05:$p06:$p07:$p08:$p09:$p10:$p11:$p12:$p13:$p14 "$SCRIPTPATH"/00.reshape.sh		| tail -n 1 |awk  '{print $1}')
s02=$(msub -l depend=$s01 "$SCRIPTPATH"/01.model.selection.sh	| tail -n 1 |awk  '{print $1}')
s03=$(msub -l depend=$s01 "$SCRIPTPATH"/03.paramestim.im.all.sh	| tail -n 1 |awk  '{print $1}')
s04=$(msub -l depend=$s01 "$SCRIPTPATH"/03.paramestim.sc.all.sh | tail -n 1 |awk  '{print $1}')
s05=$(msub -l depend=$s01 "$SCRIPTPATH"/03.paramestim.am.all.sh	| tail -n 1 |awk  '{print $1}')
s06=$(msub -l depend=$s01 "$SCRIPTPATH"/03.paramestim.si.all.sh	| tail -n 1 |awk  '{print $1}')
#s7=$(msub -l depend=$s6 "$SCRIPTPATH"/06.goodness_fit.sh	| tail -n 1 |awk  '{print $1}')

# Drawing neutral enveloppe, computing Fst and He
#s10=$(msub -l depend=$s7 "$SCRIPTPATH"/07.fst.he.empirical.sh	| tail -n 1 |awk  '{print $1}')

#confirm job submissions
echo "All jobs submitted"
