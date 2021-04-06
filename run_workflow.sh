#!/bin/bash

# parameters
DSID=$1
MZP=$2
MDH=$3
MDM=$4
GQ=$5
GX=$6

# generate signal + compute cross-section
cd generate_signal
bash run.sh $DSID $MZP $MDH $MDM $GQ $GX
cd -

# obtain product of acceptance and efficiency
# from truth-level selection in SimpleAnalysis
cd simple_analysis
bash run.sh $DSID
cd -

# compare with cross-section limits for EtMiss + h(bb)
cd evaluate_limits
bash run.sh $DSID
cd -
