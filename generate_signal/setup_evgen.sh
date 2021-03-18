#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

MYRELEASE="AthGeneration,21.6.61,here"
WORKDIR="workdir_evgen"

rm -rf $WORKDIR
mkdir $WORKDIR
cd $WORKDIR
asetup ${MYRELEASE}
cd -
