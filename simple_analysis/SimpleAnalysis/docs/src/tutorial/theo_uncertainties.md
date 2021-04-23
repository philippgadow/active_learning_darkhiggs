---
title: Theory uncertainties
---

!!! abstract "This is just a general prescription"
    This page is not a full-blown and detailled tutorial on how to compute theoretical uncertainties (as this would be too involved for this tutorial). It rather provides a general prescription that is intended to get you started with computing the theoretical signal uncertainties in your analysis.

For many analyses, SimpleAnalysis is an integral part of how they compute signal theory uncertainties. SUSY analyses usually need to consider three different kinds of systematic variations:

- Factorization and renormalization scale,

- merging scale,

- radiation, i.e. parton shower tuning or radiation uncertainty.

In this tutorial, we will not go through a full-blown example of how to compute the theory uncertainties for a given signal point, since this involves quite a bit more than only SimpleAnalysis. We will, however, introduce the general prescription and give a few starting points to get you started for computing the signal theory uncertainties in your analysis.

In general, the prescription is quite simple: For any given variation, one creates a new job option, then generates a new EVNT sample, followed by the creation of a TRUTH3 derivation using that newly created EVNT sample. The TRUTH3 derivations can then be processed with the SimpleAnalysis implementation of the analysis. The difference in the yields in region between the nominal sample and the systematic variation can then be used to derive an uncertainty for that signal point in that region.


## Event generation

Let's pick for example the merging scale up and down variation. We first need to create job options that actually vary this parameter in upward and downward direction. This is handled by giving the job options the right name. The control files will pick up the variation from the job option's name and run the correct variation.

For the signal point with the job option `MC15.376310.MGPy8EG_A14N23LO_SS_onestepCC_1000_900_800.py`, we can do event generation with merging scale up and down variation by using the following job options (with the same content):
```sh
MC15.376310.MGPy8EG_A14N23LO_SS_onestepCC_1000_900_800_qcup.py
MC15.376310.MGPy8EG_A14N23LO_SS_onestepCC_1000_900_800_qcdw.py
```

You will of course need to setup a software release for running a `pathena --trf` and command, e.g.
```sh
setupATLAS
asetup MCProd,19.2.5.36.3,here
lsetup panda
```

!!! tip "Which release to use?"
    Look up the release number you need for the e-tag you want to produce using `GetTfCommand.py --AMI=<amitag>`

Then we use the usual `pathena` command to submit a `Generate_tf` job to the grid for the nominal and the systematic variations.
```sh
pathena --trf="Generate_tf.py --ecmEnergy=13000 --runNumber=<runNumber> --firstEvent=1 --maxEvents=<maxEvents> --randomSeed=%RNDM:333 --jobConfig=<jobOption> --outputEVNTFile=%OUT.EVNT.root" --outDS=user.<username>.<gentf_outputname>
```


## TRUTH3 derivation
Once your event generation jobs are done, download the output and use it to submit derivation jobs. First setup a release again:
```sh
setupATLAS
asetup 21.2,AthDerivation,21.2.44.0
lsetup panda
```

!!! tip "Which release to use?"
    Look up the release number you need for the e-tag you want to produce using `GetTfCommand.py --AMI=<amitag>`


And then send off the derivation jobs:
```sh
pathena --trf="Reco_tf.py --maxEvents <maxEvents> --inputEVNTFile %IN --outputDAODFile %OUT.root --reductionConf TRUTH3"  --inDS=user.<username>.<gentf_outputname>_EXT0 --outDS=user.<username>.<truth3_outputname>
```

## SimpleAnalysis

The generated TRUTH3 datasets can now be put through your SimpleAnalysis analysis implementation. For each of TRUTH3 datasets (nominal and all systematic variations) for each signal point, go ahead and do the usual:
```sh
simpleAnalysis -n -a <your_analysis> <input_TRUTH3>
```

!!! tip "Validation of SimpleAnalysis output"
    It may be a good idea to create a bunch of ntuple branches for your analysis variables using `ntupVar()`. This allows you to check the distributions in case there is need to. You may find that you need to loosen some requirements in your regions of interest in order not to run into to much statistical fluctuations when computing the uncertainties. In that case, looking at the distributions helps to loosen only requirements that do not change the shapes too much.

This will give you the yields in each region for nominal and all systematic variations you processed. Computing an uncertainty value &alpha; for a given systematic variation for a given signal point in a specific region can then be as simple as:

$$
\alpha = \frac{N_\mathrm{sys}-N_\mathrm{nominal}}{N_\mathrm{nominal}}
$$

where $N_\mathrm{sys}$ and $N_\mathrm{nominal}$ are the systematic and nominal yields in that region, respectively.
