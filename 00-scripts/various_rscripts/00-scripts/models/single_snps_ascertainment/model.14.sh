#!/bin/bash

# Global variables
FOLDER=$1
NREPS=$(grep -B 1 myfifo spinput.txt |head -1)
BPFILE=../../bpfile_10
SPINPUT=../../spinput.txt

# Create folder and move into it
mkdir "$FOLDER" 2>/dev/null
cd "$FOLDER"

# Copy bpfile and spinput.txt
cp "$BPFILE" .
cp "$SPINPUT" .
NLOC=$(grep -v ^$ $SPINPUT | head -1 )
NLOCUS=$(awk '{print NF}' $BPFILE |sed -n '2p' )

# Create fifo
mknod myfifo p

# Launch ms
../../bin/priorgen5.py \
    bpfile="$BPFILE" n1=0 n1=20 n2=0 n2=20 nA=0  nA=20 tau=0 tau=30 bottleneck=N taubottle=0 taubottle=10 alpha1=1 alpha1=5 alpha2=1 alpha2=5 M1=0 M1=40 \
    M2=0 M2=40 shape1=0 shape1=20 shape2=0 shape2=200 model=IM nreps="$NREPS" \
    Nvariation=homo Mvariation=homo symMig=asym parameters=priorfile | \
    ../../bin/msnsam tbs $(( $NREPS * $NLOCUS )) -s 1 -I 2 \
        tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 2 1 -eN \
        tbs tbs | ../../bin/ascertainement_bias_4.py   $NLOC  $NLOCUS  > myfifo   &
../../bin/mscalc < myfifo #>/dev/null
rm bpfile_* error.txt  myfifo seedms  spinput.txt  spoutput.txt
