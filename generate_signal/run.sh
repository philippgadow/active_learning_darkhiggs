#!/bin/bash

# generate signal point and associated TRUTH1 derivation in one go for testing
DSID=100000
EVNT_FILE="test_DMSbb_${DSID}.EVNT.root"
DERIVATION_FILE="test_DMSbb_${DSID}.root"

# signal model parameters
# (eventually to be set via a RECAST workflow)
MZP=500
MDH=50
MDM=100
GQ=0.25
GX=1.0

# event generation
env -i
source setup_evgen.sh
python scripts/makeJobOption.py \
  --mzp $MZP \
  --mdh $MDH \
  --mdm $MDM \
  --gq $GQ \
  --gx $GX \
  --dsid $DSID \
  --outputDir workdir_evgen/
cd workdir_evgen
export ATHENA_PROC_NUMBER=8
Gen_tf.py \
    --ecmEnergy=13000. \
    --firstEvent=1 \
    --maxEvents=-1 \
    --jobConfig=100000 \
    --outputEVNTFile=$EVNT_FILE
cd -

# launch derivation
env -i
source setup_derivation.sh
cd workdir_derivation
Reco_tf.py \
    --inputEVNTFile ../workdir_evgen/$EVNT_FILE \
    --outputDAODFile $DERIVATION_FILE \
    --reductionConf TRUTH1
cd -
