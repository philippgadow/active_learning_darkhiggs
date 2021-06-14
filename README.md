## Truth level implementation of a search for dark matter and dark-Higgs bosons in the MET + s(bb) final state

### Executive summary

This project provides code to perform sensitivity estimates for a search in the MET + s(bb) final state using the ATLAS detector at the Large Hadron Collider.

The code in this project provides the means to

1. generate the dark-Higgs boson signal
2. run a truth-level implementation on the signal generation output to provide acceptance and efficiency, as well as dark Higgs candidate mass distributions
3. estimate sensitivity using generic limits on the MET + h(bb) cross-section and provide parametrised input files for a sensitivity estimate based on a pyhf RECAST

### Setup

Clone the project and make sure that your environment has access to [`/cvmfs`](https://cvmfs.readthedocs.io/en/stable/) to be able to run the signal generation step. This should be the case on most computing nodes provided by research institutions, such as CERN's `lxplus.cern.ch`.

Then, build the parts of the code which require compilation:

```
cd simple_analysis
source install.sh
```

Set up and test the python virtualenv:
```
cd evaluate_yields
source setup.sh
# check if pyhf and uproot have been installed
pip freeze
# exit virtual environment
deactivate
```


### Usage

You can run the workflow for a single point in parameter space of the signal model:
```
# DSID: 100022
# Mass of Z' boson: 2000 GeV
# Mass of dark Higgs boson: 130 GeV
# Mass of dark matter particle: 200 GeV
# Coupling gq: 0.25
# Coupling gx: 1.0

bash run_workflow.sh 100022 2000 130 200 0.25 1.0
```

You can also submit a list of jobs to a HTCondor batch system:
```
python submit_batch.py
```


### Description of components

The individual components of the project are described in the following:

#### Signal generation

The dark Higgs boson and dark matter signal is generated using the ATLAS infrastructure and job options.
The relevant parts are located in the directory `generate_signal`.

Two steps are needed to produce the final TRUTH DxAOD which is used as input by the SimpleAnalysis step.
First, the signal is generated for the specific mode configuration which the user has provided as the input.
The top level job option is provided in `generate_signal/assets/MadGraphControl_MadGraphPythia8_N31LO_A14N23LO_DMSbb_CKKWL.py`. As a result of running the first step, an `EVNT` file is generated. 
Second, the EVNT file is processed with the ATLAS derivation framework to provide a TRUTH DxAOD file which can be processed further.

The user can run both steps using the steering script `generate_signal/run.sh`

The parameters which are provided to this script are listed in the following together with the command to run the signal generation.

```
# DSID: dataset identification number, a unique six-digit number used for bookkeeping
# MZP: mass of Z' boson in GeV
# MDH: mass of dark Higgs boson in GeV
# MDM: mass of dark matter particle in GeV
# GQ: coupling of Z' boson and quarks (default: 0.25)
# GX: coupling of Z' boson and dark sector (default: 1.0)
# EVENTS: number of events to be generated (default: 10000)
bash run.sh <DSID> <MZP> <MDH> <MDM> <GQ> <GX> <EVENTS>
```

The script in turn calls `setup_evgen.sh` and `setup_derivation.sh` to run the two steps, respectively. The operations are carried out in `workdir_evgen` and `workdir_derivation`, respectively.
The resulting output TRUTH DxAOD file is located in the `data/` directory (at top level).

#### SimpleAnalysis truth-level implementation of a MET + h(bb) search

The ATLAS search for dark matter produced in association with a Higgs boson decaying to b-quarks ([ATLAS-CONF-2021-006](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2021-006/)) is implemented using truth-level objects with the [SimpleAnalysis](https://simpleanalysis.docs.cern.ch) framework.
The SimpleAnalysis framework takes care of the most relevant steps of the search implementation such as object definition, cuts and even smearing of truth level objects using ATLAS smearing functions.
The relevant parts are located in the directory `simple_analysis`.
The truth analysis itself is implemented in `simple_analysis/SimpleAnalysis/SimpleAnalysisCodes/src/ANA-EXOT-2018-046.cxx`.

Note that you need to compile the code first before you can run it:
```
cd simple_analysis
source install.sh
```

To run over a previously produced signal DxAOD with a certain DSID, execute:
```
cd simple_analysis
bash run.sh <DSID>
```

As a result, you will find the additional output files in the `data/` directory.
These consist of:

- `acceptances.txt`: text file with unweighted events, weighted events, product of acceptance and efficiency, and uncertainty for each signal region (one SR corresponds to one column)
- `histograms.root`: ROOT file with cutflow histograms, monitoring distributions such as b-tag multiplicity, MET, Higgs candidate invariant mass, as well as Higgs candidate mass distributions in respective signal regions

#### Sensitivity estimate

The `acceptances.txt` file from the preceeding step is used to compute the sensitivity of the search using generic limits on the MET + h(bb) final state.

The associated code is provided in `evaluate_limits`. The model-independent limits which are used as input together with the cross-section, the product of acceptance and efficiency provided by SimpleAnalysis, are located in `evaluate_limits/data/limits_EXOT-2018-46.csv`. The script itself is located in `evaluate_limits/scripts/estimateLimit.py`.

You can run this step individually for a DSID which has been both generated and processed with SimpleAnalysis using
```
cd evaluate_limits
bash run.sh <DSID>
```

The output is located in `data/`, as it is the case in the other steps.

#### Inputs for a pyhf RECAST

Finally, you can also provide inputs to a RECAST of the MET + h(bb) search using [`pyhf`](https://github.com/scikit-hep/pyhf).

The input required by the `pyhf` workflow (which is not part of this project) are so-called JSON patches.
These are provided by this step. The associated code is located in `evaluate_yields`.
The JSON template, which is filled with the event yields estimate by the SimpleAnalysis implementation per bin in the dark Higgs candidate mass distribution, is located in `evaluate_yields/data/patch_template.json`.
The script is located in `evaluate_yields/scripts/getYields.py`. As it depends on [`uproot`](https://uproot.readthedocs.io/en/latest/) to avoid the pain and misery typically associated when working with vanilla ROOT, it is necessary to set up a virtual environment and install `uproot` locally:

```
cd evaluate_yields
source setup.sh
# DSID: dataset identification number, a unique six-digit number used for bookkeeping
# MZP: mass of Z' boson in GeV
# MDH: mass of dark Higgs boson in GeV
# MDM: mass of dark matter particle in GeV
# GQ: coupling of Z' boson and quarks (default: 0.25)
# GX: coupling of Z' boson and dark sector (default: 1.0)
bash run.sh <DSID> <MZP> <MDH> <MDM> <GQ> <GX>
```

The output is located in `data/`, as it is the case in the other steps.

### Literature

1. Dark Higgs boson theory paper: https://arxiv.org/abs/1701.08780 