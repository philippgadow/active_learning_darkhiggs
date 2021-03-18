export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup git
asetup AnalysisBase,21.2.100
git clone --recursive https://:@gitlab.cern.ch:8443/atlas-phys-susy-wg/SimpleAnalysis.git
cd SimpleAnalysis/

#Compile
mkdir build; mkdir run
cd build
#With b-jet weighting
cmake -DDO_TRUTHTAGGING=ON ..

make
source x86_64-centos7-gcc8-opt/setup.sh
cd ../run
