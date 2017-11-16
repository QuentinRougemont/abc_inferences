#!/bin/bash

# Global variables
FOLDER=$1
NREPS=10
BPFILE=../../bpfile
SPINPUT=../../spinput.txt
#NLOCUS=5612

# Create folder and move into it
mkdir "$FOLDER" 2>/dev/null
cd "$FOLDER"

# Copy bpfile and spinput.txt
cp "$BPFILE" .
cp "$SPINPUT" .

# Create fifo
mknod myfifo p

# Launch ms
../../bin/mscalc < myfifo #>/dev/null &
../../bin/priorgen5.py \
    bpfile="$BPFILE" n1=0 n1=20 n2=0 n2=20 nA=0 nA=20 tau=0 tau=30 bottleneck=Y \
    taubottle=0 taubottle=10 alpha1=1 alpha1=5 alpha2=1 alpha2=5 M1=0 M1=40 \
    M2=0 M2=40 shape1=0 shape1=20 shape2=0 shape2=200 model=PSC nreps="$NREPS" \
    Nvariation=hetero Mvariation=hetero symMig=asym parameters=priorfile | \
    ../../bin/msnsam tbs $(( $NREPS * $NLOC )) -t 0.08 -r 0.08 80 -I 2 \
        tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -eM tbs 0 -en tbs 1 tbs \
        -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -eM tbs 0 -ej tbs 2 1 -eN tbs tbs > myfifo
