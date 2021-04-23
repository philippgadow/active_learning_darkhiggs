#!/bin/bash
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
asetup AnalysisBase,21.2.169

cd SimpleAnalysis/build
source x86_64-centos7-gcc8-opt/setup.sh
cd ../..
