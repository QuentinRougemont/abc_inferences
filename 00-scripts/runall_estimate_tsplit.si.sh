#!/bin/bash
# WARNING! All scripts must been edited prior to submitting and not changed until run!
# HINT: Give meaningful job names in submission scripts to distinguish them

# Global variables
SCRIPTPATH="./00-scripts"

echo "    DO NOT USE WITHOUT APPROPRIATE ABC BACKGROUND"

# Estimate Tsplit under SI
s01=$(msub  "$SCRIPTPATH"/10.prepare.si.input.sh  | tail -n 1 | awk '{print $1}') #this script will use a subset of the most differentiated snps (eg: upper 90 or 95% quantile or even 99% depending on the threshild you defined
#s02=$(msub "$DEPENDS"$s14 "$SCRIPTPATH"/11.prepare.si.architecture.sh | tail -n 1 | awk '{print $1}') #create a new archietecture to run the job

#Then you will have to run the script run_all.job.sh or this:
#This is necessary to make sure that the SI is the best fitting model according to the defined quantile 
p01=$(msub  "$SCRIPTPATH"/models/model_job_array.1.sh  | tail -n 1 | awk '{print $1}') #si n.hetero
p02=$(msub  "$SCRIPTPATH"/models/model_job_array.2.sh  | tail -n 1 | awk '{print $1}')  #si n.homo
p03=$(msub  "$SCRIPTPATH"/models/model_job_array.3.sh  | tail -n 1 | awk '{print $1}')  #am n.hetero.m.hetero
p04=$(msub  "$SCRIPTPATH"/models/model_job_array.4.sh  | tail -n 1 | awk '{print $1}')  #am n.homo.m.hetero
p05=$(msub  "$SCRIPTPATH"/models/model_job_array.5.sh  | tail -n 1 | awk '{print $1}')  #am n.hetero.m.homo
p06=$(msub  "$SCRIPTPATH"/models/model_job_array.6.sh  | tail -n 1 | awk '{print $1}')  #am n.homo.m.homo
p07=$(msub  "$SCRIPTPATH"/models/model_job_array.7.sh  | tail -n 1 | awk '{print $1}')  #sc.heterom.heteron
p08=$(msub  "$SCRIPTPATH"/models/model_job_array.8.sh  | tail -n 1 | awk '{print $1}')  #sc.heterom.homon
p09=$(msub  "$SCRIPTPATH"/models/model_job_array.9.sh  | tail -n 1 | awk '{print $1}')  #sc.homom.heteron
p10=$(msub  "$SCRIPTPATH"/models/model_job_array.10.sh  | tail -n 1 | awk '{print $1}') #im.heterom.heteron
p11=$(msub  "$SCRIPTPATH"/models/model_job_array.11.sh  | tail -n 1 | awk '{print $1}') #im.heterom.homon
p12=$(msub  "$SCRIPTPATH"/models/model_job_array.12.sh  | tail -n 1 | awk '{print $1}') #im.homom.heteron
p13=$(msub  "$SCRIPTPATH"/models/model_job_array.13.sh  | tail -n 1 | awk '{print $1}') #sc.homom.homon
p14=$(msub  "$SCRIPTPATH"/models/model_job_array.14.sh  | tail -n 1 | awk '{print $1}') #im.homom.homon
#p15=$(msub  "$SCRIPTPATH"/models/model_job_array.15.sh  | tail -n 1 | awk '{print $1}') #bot.sc.heterom.heteron
#p16=$(msub  "$SCRIPTPATH"/models/model_job_array.16.sh  | tail -n 1 | awk '{print $1}') #2P.sc.heterom.heteron
#p17=$(msub  "$SCRIPTPATH"/models/model_job_array.17.sh  | tail -n 1 | awk '{print $1}') #2P.bot.sc.heterom.heteron

#reshape the simul and process to model choice
s01=$(msub -l depend=$p01:$p02:$p03:$p04:$p05:$p06:$p07:$p08:$p09:$p10:$p11:$p12:$p13:$p14 "$SCRIPTPATH"/00.reshape.sh | tail -n 1 | awk '{print $1}')
s02=$(msub -l depend=$s01 "$SCRIPTPATH"/03.model.selection.sh  | tail -n 1 | awk '{print $1}')

#Once the good threshold is find you can estimate parameters under SI:
s03=$(msub -l depend=$s01 "$SCRIPTPATH"/04.paramestim.si.all.sh | tail -n 1 | awk '{print $1}')
#s04=$(msub -l depend=$s01 "$SCRIPTPATH"/04.paramestim.sc.all.sh | tail -n 1 | awk '{print $1}')
#s04=$(msub "$DEPENDS"$s03 "$SCRIPTPATH"/06.goodness_fit.si.sh  | tail -n 1 | awk '{print $1}')

#confirm job submissions
echo "All jobs submitted"
#end of script
