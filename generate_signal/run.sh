#!/bin/bash
##############
# parameters #
##############
set -x

DSID=$1
MZP=$2
MDH=$3
MDM=$4
GQ=$5
GX=$6
EVENTS=$7

# generate signal point and associated TRUTH1 derivation in one go for testing
EVNT_FILE="test_DMSbb_${DSID}.EVNT.root"
DERIVATION_FILE="test_DMSbb_${DSID}.root"

# define directories
export EVGEN_WORKDIR="workdir_evgen/${DSID}"
export DERIVATION_WORKDIR="workdir_derivation/${DSID}"
OUTPUT_DIR="../data/${DSID}"
mkdir -p ${OUTPUT_DIR} ${EVGEN_WORKDIR} ${DERIVATION_WORKDIR}

# event generation
source setup_evgen.sh
python scripts/makeJobOption.py \
  --mzp $MZP \
  --mdh $MDH \
  --mdm $MDM \
  --gq $GQ \
  --gx $GX \
  --dsid $DSID \
  --outputDir $EVGEN_WORKDIR/
cd $EVGEN_WORKDIR
Gen_tf.py \
    --ecmEnergy=13000. \
    --firstEvent=1 \
    --maxEvents=$EVENTS \
    --jobConfig=$DSID \
    --outputEVNTFile=$EVNT_FILE
cd -

# launch derivation
source setup_derivation.sh
cd $DERIVATION_WORKDIR
Reco_tf.py \
    --inputEVNTFile ../../$EVGEN_WORKDIR/$EVNT_FILE \
    --outputDAODFile $DERIVATION_FILE \
    --reductionConf TRUTH1
cd -

# extract cross-section
python scripts/extractRunInfo.py -i "$EVGEN_WORKDIR/log.generate" -o "${OUTPUT_DIR}/xsec.json"

# copy derivation file to output directory
cp $DERIVATION_WORKDIR/DAOD_TRUTH1.test_DMSbb_${DSID}.root ${OUTPUT_DIR}/

# copy log file to output directory
cp $EVGEN_WORKDIR/log.generate ${OUTPUT_DIR}/

# clean up
rm -rf ${EVGEN_WORKDIR}/
rm -rf ${DERIVATION_WORKDIR}/
