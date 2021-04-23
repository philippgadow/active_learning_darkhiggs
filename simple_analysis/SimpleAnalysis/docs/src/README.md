# Simplified ATLAS SUSY analysis framework

Holds a collections of SUSY analyses. These can be run over samples in different
input formats:

 * DAOD_TRUTH      (TRUTH1 and TRUTH3 tested)
 * xAOD            (either truth-level and or reco-level - the latter with some constraints)
 * slimmed ntuples (Reduced ntuples produced from above input)

It provides the analysis acceptance per control and signal region as well as
optionally histograms or ntuples with event level objects. Smearing of the truth-level objects
can optionally be done using the latest upgrade performance smearing functions.

## Docs

Documentation is found at [https://simpleanalysis.docs.cern.ch](https://simpleanalysis.docs.cern.ch).

## Compiling

It should compile out-of-the-box on top of any recent release on top of CentOS7.

```bash
setupATLAS
lsetup git
asetup AnalysisBase,21.2.100
git clone --recursive https://:@gitlab.cern.ch:8443/atlas-phys-susy-wg/SimpleAnalysis.git
cd SimpleAnalysis/

#Compile
mkdir build; mkdir run
cd build
#Without b-jet weighting:
cmake ..
#With b-jet weighting
cmake -DDO_TRUTHTAGGING=ON ..

make
source x86_64-centos7-gcc8-opt/setup.sh
cd ../run
```

### Submission on the grid

Submission to the CERN grid is fairly straightforward:

```bash
lsetup panda
mkdir submit
cd submit
ln -s ../build/x86_64-centos7-gcc8-opt
prun --osMatching --exec 'source x86_64-centos7-gcc8-opt/setup.sh;simpleAnalysis -a ZeroLepton2015,ThreeBjets2015 %IN' --outputs='*.txt,*.root' --extFile \*.root --athenaTag 21.2.100,AnalysisBase --maxFileSize 40000000 --noCompile --followLinks --inDS <inputDS> --outDS user.<userName>.simple.v1
```

## Adding a new analysis

Simply make new C++ routine in `SimpleAnalysisCodes/src/MyAnalysisName.cxx` with the structure:

```bash
#include "SimpleAnalysisFramework/AnalysisClass.h"
DefineAnalysis(MyAnalysisName)
void MyAnalysisName::Init() {}
void MyAnalysis::ProcessEvent(AnalysisEvent *event)
```

Recompilation will automatically make the analysis available to run. See
`ExampleAnalysis.cxx` and the other supplied analysis examples on
how to specify signal regions, kinematic selections and filling
histograms/ntuples.

## Submission of new SUSY analyses

It is suggested to use analysis names that ends in the year of the latest used
data, i.e. "2015" or "2016", and to add "Conf" in front if it is a preliminary
analysis. Examples: `ZeroLepton2016.cxx` or `StopOneLeptonConf2018.cxx`

To submit either fork the main repository in gitlab, push changes to the private
fork and submit a merge request or simply email the relevant code to @aagaard

For more help, please contact @aagaard.

## Developing

### clang-format

Code style needs to adhere to the [`.clang-format`](./clang-format). This can be done by running `clang-format` locally, or using the docker image

```
docker run -it --rm -v "$(pwd)":/workdir -w /workdir kratsg/clang-format -r . -i
```

### documentation

First, clone the repository. If you want to contribute, just edit the respective file(s) and open a Merge Request. An OpenShift instance will then deal with building and deploying a staging version of your changes for you to preview. If everything looks good, it will be merged in and auto deploy via OpenShift.

To set up a development server at [http://localhost:8000](http://localhost:8000) for the public documentation:

```
docker run -it --rm -v $PWD:/opt/app-root/src \
                    -p 8000:8000 \
                    gitlab-registry.cern.ch/authoring/documentation/s2i-mkdocs-container \
                    bash -c 'cd docs/public; pip install -r requirements.txt; mkdocs serve --dev-addr 0.0.0.0:8000'
```

or for the internal documentation:

```
docker run -it --rm -v $PWD:/opt/app-root/src \
                    -p 8000:8000 \
                    gitlab-registry.cern.ch/authoring/documentation/s2i-mkdocs-container \
                    bash -c 'cd docs/internal; pip install -r requirements.txt; mkdocs serve --dev-addr 0.0.0.0:8000'
```


What's really nice about this method is that it will auto-detect changes and rebuild the site for you as you go.

#### Reference

For reference on how to make these docs, visit [mkdocs-material](https://squidfunk.github.io/mkdocs-material/) which is the theme used and [how-to.docs @ CERN](https://how-to.docs.cern.ch/) for how they're hosted at CERN from GitLab.
