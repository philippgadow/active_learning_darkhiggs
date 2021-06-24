#!/bin/bash

##############
# parameters #
##############

DSID=$1

# definitions
DATA_DIR=$PWD/../data/${DSID}

python scripts/estimateLimit.py \
    --acceptances ${DATA_DIR}/acceptances.txt \
    --cross_section ${DATA_DIR}/xsec.json \
    --output_file ${DATA_DIR}/sensitivity.txt
