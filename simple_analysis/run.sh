#!/bin/bash

##############
# parameters #
##############
# # DSID - just for bookkeeping
# if [[ -z "${DSID}" ]]; then
#   DSID=100000
# fi

DSID=$1

# define data directory and analysis
ANALYSIS="METHbb2018"
INPUT_DIR="$PWD/../data/${DSID}"
INPUT_FILE="${INPUT_DIR}/DAOD_TRUTH1.test_DMSbb_${DSID}.root"
RUN_DIR=SimpleAnalysis/run/$DSID

# setup simple analysis
source setup.sh

# run over file
mkdir -p $RUN_DIR
cd $RUN_DIR
simpleAnalysis -s layout=run2 -a ${ANALYSIS} ${INPUT_FILE}
cd -

# copy file with acceptances to output directory
cp $RUN_DIR/${ANALYSIS}.txt ${INPUT_DIR}/acceptances.txt
cp $RUN_DIR/${ANALYSIS}.root ${INPUT_DIR}/histograms.root
cat ${INPUT_DIR}/acceptances.txt

# clean up
rm -rf ${RUN_DIR}
