#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

MYRELEASE="AthDerivation,21.2.118.0,here"
if [[ -z "${DERIVATION_WORKDIR}" ]]; then
  DERIVATION_WORKDIR="workdir_derivation"
fi

rm -rf $DERIVATION_WORKDIR
mkdir $DERIVATION_WORKDIR
cd $DERIVATION_WORKDIR
asetup ${MYRELEASE}
cd -
