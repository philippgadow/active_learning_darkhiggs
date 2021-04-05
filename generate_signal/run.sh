#!/bin/bash
##############
# parameters #
##############
# # DSID - just for bookkeeping
# if [[ -z "${DSID}" ]]; then
#   DSID=100000
# fi

# # signal model parameters
# # Z' boson mass
# if [[ -z "${MZP}" ]]; then
#   MZP=2500
# fi

# # dark Higgs boson mass
# if [[ -z "${MDH}" ]]; then
#   MDH=150
# fi

# # dark matter particle mass
# if [[ -z "${MDM}" ]]; then
#   MDM=200
# fi

# # coupling of Z' boson to quarks
# if [[ -z "${GQ}" ]]; then
#   GQ=0.25
# fi

# # dark sector coupling
# if [[ -z "${GX}" ]]; then
#   GX=1.0
# fi

DSID=$1
MZP=$2
MDH=$3
MDM=$4
GQ=$5
GX=$6

# generate signal point and associated TRUTH1 derivation in one go for testing
EVNT_FILE="test_DMSbb_${DSID}.EVNT.root"
DERIVATION_FILE="test_DMSbb_${DSID}.root"

# define directories
export EVGEN_WORKDIR="workdir_evgen/${DSID}"
export DERIVATION_WORKDIR="workdir_derivation/${DSID}"
OUTPUT_DIR="../data/${DSID}"
mkdir -p ${OUTPUT_DIR} ${EVGEN_DIR} ${DERIVATION_WORKDIR}

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
    --maxEvents=-1 \
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
cp $DERIVATION_WORKDIR/DAOD_TRUTH1.test_DMSbb_100000.root ${OUTPUT_DIR}/
