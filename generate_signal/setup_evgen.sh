#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

MYRELEASE="AthGeneration,21.6.61,here"
if [[ -z "${EVGEN_WORKDIR}" ]]; then
  EVGEN_WORKDIR="workdir_evgen"
fi


rm -rf $EVGEN_WORKDIR
mkdir $EVGEN_WORKDIR
cd $EVGEN_WORKDIR
asetup ${MYRELEASE}
cd -
