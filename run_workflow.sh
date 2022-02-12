#!/bin/bash

# parameters
DSID=$1
MZP=$2
MDH=$3
MDM=$4
GQ=$5
GX=$6
EVENTS=10000
EVENTS=1000

# generate signal + compute cross-section
cd generate_signal
bash run.sh $DSID $MZP $MDH $MDM $GQ $GX $EVENTS
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

# create json file with signal event yields per bin for pyhf
cd evaluate_yields
bash run.sh $DSID $MZP $MDH $MDM $GQ $GX $EVENTS
cd -

# clean up to save some space
rm -rf generate_signal/workdir_evgen/${DSID}/
rm -rf generate_signal/workdir_derivation/${DSID}/
rm data/${DSID}/DAOD_TRUTH1.test_DMSbb_${DSID}.root
