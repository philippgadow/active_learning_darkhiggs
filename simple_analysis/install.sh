#!/bin/bash
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup git
asetup AnalysisBase,21.2.169

if [ ! -d "SimpleAnalysis" ] 
then
    echo "SimpleAnalysis DOES NOT exists. Cloning it" 
    git clone --recursive https://:@gitlab.cern.ch:8443/atlas-phys-susy-wg/SimpleAnalysis.git
fi

cd SimpleAnalysis/

#Compile
mkdir -p build run
cd build
#With b-jet weighting
cmake -DDO_TRUTHTAGGING=ON ..
make

source x86_64-centos7-gcc8-opt/setup.sh
cd ../..
