#!/bin/bash

##############
# parameters #
##############

DSID=$1
MZP=$2
MDH=$3
MDM=$4
GQ=$5
GX=$6

# definitions
DATA_DIR=$PWD/../data/${DSID}

source setup.sh
python scripts/getYields.py ${DATA_DIR}/histograms.root \
    --acceptances ${DATA_DIR}/acceptances.txt \
    --cross_section ${DATA_DIR}/xsec.json \
    --signalname monoSbb_mzp${MZP}_mdh${MDH}_mdm${MDM}_gq${GQ}_gx${GX} \
    --mzp ${MZP} \
    --mdh ${MDH} \
    --mdm ${MDM} \
    --gq ${GQ} \
    --gx ${GX} \
    --output_file ${DATA_DIR}/patch.json
deactivate

cat ${DATA_DIR}/patch.json
