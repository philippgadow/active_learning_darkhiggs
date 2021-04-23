---
title: Cross-sections
---

You might be wondering how to scale SimpleAnalaysis output to a given cross-section and/or luminosity. There are two main ways of doing this, depending on the output format you want to use downstream.

!!! warning "Cross-sections in SimpleAnalysis"
    Please refrain from putting cross-section treatment into your SimpleAnalysis implementation file itself. This would unnecessarily clutter a file that is thought to be your analysis selection only.

## Acceptance numbers

If you use acceptance numbers and you want to know, given a certain integrated luminosity, how many events this equals to, you will need to compute:

$$
N_{\mathrm{truth}} = a\cdot\sigma\cdot\epsilon_\mathrm{filter}\cdot\mathcal{L}
$$

where *a* is the acceptance from SimpleAnalysis, $\sigma$ is the cross-section of the process, $\epsilon_\mathrm{filter}$ the MC filter efficiency and $\mathcal{L}$ the integrated luminosity.

You can get the cross-section from e.g. the PMG cross-section database in `/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PMGTools/PMGxsecDB_mc16.txt`.

### C++

!!! info "Setting up a software release"
    This requires you to have a software release setup, which is often something you don't want to do, especially if it is only for the sake of getting the right cross-sections.

If you are using C++ downstream of SimpleAnalysis, just setup the `SUSYTools` cross-section tool, which relies on the PMG cross-section database and offers a number of convenience methods you can use:
```cpp
SUSY::CrossSectionDB *xsecDB = new SUSY::CrossSectionDB(PathResolverFindCalibFile("dev/PMGTools/PMGxsecDB_mc16.txt"), true);
double xsec_times_eff = xsecDB->xsectTimesEff(396681);
```

### Python

!!! info "Pandas"
    The below solution uses `pandas` for storing the central cross-section database in memory. Given the size of the file, it would also still be possible to just create an ordinary python dictionary. If you want to stick with the `pandas` solution, you might have to run `python3 -m pip install pandas --user`.

Often times, you'll find yourself using python (it's personal preference at the end of the day) for most of the analysis-related things downstream of ntuple productions. We can do something similar as in the C++ case. First, let's create a directory and a couple of files we need:
```sh
# in the tutorial base dir, set up env variable
export TUTORIAL_DIR=$(pwd)

mkdir $TUTORIAL_DIR/cross_section
cd $TUTORIAL_DIR/cross_section
touch CrossSectionDB.py
```

As implementation for the cross-section DB, you could for example use something like this:
```python
#!/usr/bin/env python

import os
import pandas as pd


class CrossSectionDB:
    """
    A class for getting cross-sections (and filter efficiencies) from the central PMG xsec database.
    """

    db = pd.DataFrame()

    def __init__(
        self,
        dirname="/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/PMGTools/",
        filename="PMGxsecDB_mc16.txt",
        columns=[
            "dataset_number",
            "physics_short",
            "crossSection",
            "genFiltEff",
            "kFactor",
            "relUncertUP",
            "relUncertDOWN",
            "generator_name",
            "etag",
        ],
        separator="\t\t",
    ):

        self.filename = filename
        self.dirname = dirname

        self.db = pd.read_csv(
            os.path.join(dirname, filename),
            sep=separator,
            skiprows=1,
            header=None,
            names=columns,
            engine='python',
        )

    def getMatch(self, field, dsid, etag=None):
        dsid = int(dsid)
        try:
            if etag:
                mask = (self.db["dataset_number"] == dsid) & (self.db["etag"] == etag)
            else:
                mask = (self.db["dataset_number"] == dsid)
            df = self.db.loc[mask, :]
        except Exception as e: # work on python 3.x
            print('Failed to get value: '+ str(e))
            return None

        if len(df.index) > 1:
            print(
                "More than one row with DSID "
                + str(dsid)
                + " found! Picking first one ..."
            )

        return df[field].values[0]

    def xsec(self, dsid, etag=None):
        return self.getMatch("crossSection", dsid, etag)

    def efficiency(self, dsid, etag=None):
        return self.getMatch("genFiltEff", dsid, etag)

    def kFactor(self, dsid, etag=None):
        return self.getMatch("kFactor", dsid, etag)

    def xsecTimesEff(self, dsid, etag=None):
        try:
            return self.xsec(dsid, etag) * self.efficiency(dsid, etag)
        except:
            return None

    def xsecTimesEffTimeskFac(self, dsid, etag=None):
        try:
            return self.xsec(dsid, etag) * self.efficiency(dsid, etag) * self.kFactor(dsid, etag)
        except:
            return None

```

Let's put this into the newly created `CrossSectionDB.py` file.

Next, create another python file:
```sh
touch doStuff.py
```

We'll use that file to actually do something with the SimpleAnalysis output. Starting with getting the SimpleAnalysis output from the previous job we ran (make sure you have `$TUTORIAL_DIR` set), then just scale with the cross-section, filter efficiency and luminosity
```python
import os
from CrossSectionDB import CrossSectionDB

dsid = 396776
lumi_pb = 139000

# xsecDB instance and xsec for this DSID
xsecDB = CrossSectionDB()
xsecTimesEff = xsecDB.xsecTimesEff(dsid=dsid)

# open SimpleAnalysis output and print events. Make sure you have $TUTORIAL_DIR set as env
input = os.path.join(os.getenv("TUTORIAL_DIR"),"run","MyAnalysisName.txt")
with open(input) as f:
    # skip header
    for _ in range(2):
        f.next()
    # print number of events for each ROI
    for line in f:
        region = line.split(",")[0]
        acceptance = float(line.split(",")[2])
        events = acceptance*xsecTimesEff*lumi_pb

        print(region,events)
```

Now you can simply run with
```sh
python doStuff.py
```

## Ntuples

If you are using the ntuples and not the text files containing acceptance values, you can use the `-x [--xsReweight]` option:
```sh
-x [ --xsReweight ] arg     XS reweight nEvents to LUMI '<LUMI
                            (ifb)>,<nEvents>[,<XS DB file>]'
```
to re-weight `nEvents` to `Lumi`:

```sh
simpleAnalysis -n -a MyAnalysisName -x 6000,139000 $TUTORIAL_DIR/inputs/DAOD_TRUTH3/C1N2_Wh_hbb_700p0_0p0_lep_DAOD_TRUTH3.root
```
This will grab the cross-section and filter efficiency from the PMG cross-section database and re-weight the MC event weight with it. All your histograms should now be correctly normalised to your desired luminosity and cross-section. In the ntuple, the `eventWeight` can be used for normalisation.

![TBrowser output](images/weighted_hist.png)
