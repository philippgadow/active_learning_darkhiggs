#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

MYRELEASE="AthDerivation,21.2.65.0,here"
WORKDIR="workdir_derivation"

rm -rf $WORKDIR
mkdir $WORKDIR
cd $WORKDIR
asetup ${MYRELEASE}
cd -
